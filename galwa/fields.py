from typing import List, Dict, Tuple, Optional
import numpy as np

from galwa.elements import FiniteFieldElement
from galwa.utils import valid_repr, refactor_polynom_terms, pad_element, find_all_dividers_of_field_size


class FiniteField(object):
    """
    FiniteField class represents an extension field of the form :math:`F_{p^n}` where :math:`p` is a prime number and n is the degree of \
    polynomial :math:`f(x)` used for the extension.

    .. math::
        l = F_{p^n} = F_p[x] / f(x)

    Example:

    >>> from galwa import FiniteField
    >>> import numpy as np
    >>> p = 2
    >>> f = np.array([1, 1, 1])  # x^2 + x + 1
    >>> F = FiniteField(p, f)
    >>> F.elements
    [FiniteFieldElement(0, f(x)= 1 + x + x², p=2), FiniteFieldElement(1, f(x)= 1 + x + x², p=2), FiniteFieldElement(x, f(x)= 1 + x + x², p=2), FiniteFieldElement(1 + x, f(x)= 1 + x + x², p=2)]

    """

    def __init__(self, p, f: np.ndarray, representation: Optional[str] = "polynomial"):
        """
        Initialize the FiniteField class.

        Args:
            p (int): the prime number for the field.
            f (np.ndarray): the irreducible polynomial f(x) used for the extension.
            representation Optional[str]: the representation of the elements in the field - polynomial, vector, or matrix, default is polynomial.

        Raises:
            AssertionError: if the polynomial f(x) is reducible over F_p or the representation is invalid
        """
        assert representation in ["polynomial", "vector", "matrix"], "Invalid representation"
        self._repr = representation
        self.p = p
        self.f = f % p
        assert self._check_that_f_is_irreducible(), f"f(x) is reducible over F{p}"
        self.poly_f = np.polynomial.Polynomial(self.f)
        self.x_powers = self._calculate_illegal_powers()
        self._elements = self._create_elements()
        self._generators = None
        self.field_size_dividers = self._find_all_dividers_of_field_size()

    def __eq__(self, other) -> bool:
        """
        Equality builtin function for FiniteField.

        Args:
            other (FiniteField): the other field to compare to.

        Returns:
            bool: True if the fields are equal, False otherwise.
        """
        return self.p == other.p and np.array_equal(self.f, other.f)

    def __repr__(self) -> str:
        """
        Representation builtin function for FiniteField
        Prints the representation of the field - polynomial, vector

        Returns:
            str: the representation of the field
        """
        if self._repr == "polynomial":
            poly_repr = refactor_polynom_terms(f"{self.poly_f}")
            repr_str = f"{self.__class__.__name__}(p= {self.p}, f(x)= {poly_repr})"
            return repr_str
        else:
            return f"{self.__class__.__name__}(p= {self.p}, f(x)= {self.f})"

    @property
    def representation(self) -> str:
        """
        Get the representation of the elements in the field - polynomial, vector, or matrix

        Returns:
            str: the representation of the elements in the field

        Example:

        >>> from galwa import FiniteField
        >>> import numpy as np
        >>> p = 2
        >>> f = np.array([1, 1, 1])
        >>> F = FiniteField(p, f)
        >>> F.representation
        'polynomial'
        """
        return self._repr

    @representation.setter
    @valid_repr
    def representation(self, value) -> None:
        """
        Sets the representation of the elements in the field - polynomial, vector, or matrix. \
        So once printing the elements, they will be printed as polynomials, vectors, or matrices.

        Args:
            value (str): the representation to set. (polynomial, vector, matrix)

        Returns:
            None
        """
        self._repr = value
        for element in self._elements.values():
            element.representation = value

    def elements_as_vectors(self) -> None:
        """
        Changes the representation of the elements to vectors.

        Returns:
            None

        Example:

        >>> from galwa import FiniteField
        >>> import numpy as np
        >>> p = 2
        >>> f = np.array([1, 1, 1])
        >>> F = FiniteField(p, f)
        >>> F.elements_as_vectors()
        >>> F.elements
        [FiniteFieldElement([0 0], f(x) = [1 1 1] p=2), FiniteFieldElement([1 0], f(x) = [1 1 1] p=2), FiniteFieldElement([0 1], f(x) = [1 1 1] p=2), FiniteFieldElement([1 1], f(x) = [1 1 1] p=2)]
        """
        self.representation = "vector"

    def elements_as_matrices(self) -> None:
        """
        Changes the representation of the elements to matrices.

        Returns:
            None

        Example:

        >>> from galwa import FiniteField
        >>> import numpy as np
        >>> p = 2
        >>> f = np.array([1, 1, 1])
        >>> F = FiniteField(p, f)
        >>> F.elements_as_matrices()
        >>> F.elements

        [FiniteFieldElement(None, f(x) = [1 1 1] p=2),
         FiniteFieldElement([[1 0][0 1]], f(x) = [1 1 1] p=2),
         FiniteFieldElement([[0 1][1 1]], f(x) = [1 1 1] p=2),
         FiniteFieldElement([[1 1][1 0]], f(x) = [1 1 1] p=2)]
        """
        self.representation = "matrix"

    def elements_as_polynomials(self) -> None:
        """
        Changes the representation of the elements to polynomials.

        Returns:
            None

        Example:

        >>> from galwa import FiniteField
        >>> import numpy as np
        >>> p = 2
        >>> f = np.array([1, 1, 1])
        >>> F = FiniteField(p, f)
        >>> F.elements_as_polynomials()
        >>> F.elements
        [FiniteFieldElement(0, f(x)= 1 + x + x², p=2), FiniteFieldElement(1, f(x)= 1 + x + x², p=2), FiniteFieldElement(x, f(x)= 1 + x + x², p=2), FiniteFieldElement(1 + x, f(x)= 1 + x + x², p=2)]
        """
        self.representation = "polynomial"
    
    def _find_all_dividers_of_field_size(self):
        """
        Find all the dividers of the size of the extension field.

        Returns:
            List[int]: sorted array with all the dividers of the extension field
        """
        return find_all_dividers_of_field_size(self)

    def _check_that_f_is_irreducible(self) -> bool:
        """
        Checks if the polynomial :math:`f(x)` is irreducible over the field :math:`F_{p}` for degree 2 or 3 \
        For bigger degrees we assume that the polynomial is irreducible.
        A polynomial :math:`f(x)` of degree 2 or 3 is irreducible iff it has no roots in the field :math:`F_{p}`.

        Returns:
            bool: True if the polynomial is irreducible, False otherwise.
        """
        if len(self.f) > 4:
            return True
        for i in range(self.p):
            possible_root = sum([self.f[j] * (i ** j) for j in range(len(self.f))]) % self.p
            if possible_root == 0:
                return False
        return True

    @property
    def generators(self) -> List[FiniteFieldElement]:
        """
        Property to get all the generators in the field

        .. note::
            At the first call of this property the generators are calculated and stored in the generators attribute \
            The reason for this is that the calculation of the generators is an expensive operation.
            First call might be slow depend on the p value, but once the generators are calculated they are stored \
             and can be accessed quickly.

        Returns:
            List[FiniteFieldElement]: a list of all the generators in the field

        Example:

        >>> from galwa import FiniteField
        >>> import numpy as np
        >>> p = 2
        >>> f = np.array([1, 1, 1])
        >>> F = FiniteField(p, f)
        >>> F.generators
        [FiniteFieldElement(x, f(x)= 1 + x + x², p=2), FiniteFieldElement(1 + x, f(x)= 1 + x + x², p=2)]
        """
        if self._generators is None:
            self._generators = self._find_generators()
        return self._generators

    @property
    def order(self) -> int:
        """
        Property to get the order of the field, which is :math:`p^n` where :math:`n` is the degree of the polynomial \
         :math:`f(x)`.

        Returns:
            int: the order of the field

        Example:

        >>> from galwa import FiniteField
        >>> import numpy as np
        >>> p = 2
        >>> f = np.array([1, 1, 1])
        >>> F = FiniteField(p, f)
        >>> F.order
        4
        """
        return self.p ** (len(self.f) - 1)

    def _find_generators(self):
        """
        Finds all the generators in the field.
        A generator element is an element whose order is equal to the order of the multiplicative group
        In other words, a generator element is an element whose powers generate all the elements in the field.

        .. math::
           <g> = \{g^0, g^1, g^2, ..., g^{p^n-2}\} = l^{x} = l - \{0\}


        Returns:
            List[FiniteFieldElement]: a list of all the generators in the field
        """
        generators = []
        for element in self.elements:
            if element.is_generator():  # -1 because p^r-1 is the order of the multiplicative group
                generators.append(element)
        return generators

    def _create_elements(self) -> Dict[Tuple[np.ndarray], FiniteFieldElement]:
        """
        Create all possible elements in the field.

        The elements are of the form :math:`a_{0} + a_{1}x + a_{2}x^2 + ... + a_{n-1}x^{n-1}` where :math:`n` is\
         the degree of :math:`f(x)`.

        So we create a permutations of all possible coefficients for the elements in the field and use them\
         to create the elements.

        Returns:
            Dict[Tuple[np.ndarray], FiniteFieldElement]: a dictionary of all elements in the field where the key \
            is the vector representation of the element.

        """
        n = len(self.f) - 1
        coefficients = np.arange(self.p)
        coefficient_combinations = np.meshgrid(*([coefficients] * n), indexing='ij')
        coefficient_permutations = np.stack(coefficient_combinations, axis=-1).reshape(-1, n)
        elements = {tuple(coefficient[::-1]): FiniteFieldElement(coefficient[::-1], self, representation=self._repr) for
                    coefficient in
                    coefficient_permutations}
        return elements

    @property
    def elements(self) -> List[FiniteFieldElement]:
        """
        Returns all elements in the field as a list.

        Returns:
            List[FiniteFieldElement]: a list of all the elements in the field.

        Example:

        >>> from galwa import FiniteField
        >>> import numpy as np
        >>> p = 2
        >>> f = np.array([1, 1, 1])
        >>> F = FiniteField(p, f)
        >>> F.elements
        [FiniteFieldElement(0, f(x)= 1 + x + x², p=2), FiniteFieldElement(1, f(x)= 1 + x + x², p=2), FiniteFieldElement(x, f(x)= 1 + x + x², p=2), FiniteFieldElement(1 + x, f(x)= 1 + x + x², p=2)]
        """
        return list(self._elements.values())

    def get_element(self, a: np.ndarray) -> FiniteFieldElement:
        """
        Given a vector representation of an element in the field, returns the element.

        Args:
            a (np.ndarray): the vector representation of the element.

        Returns:
            FiniteFieldElement: the element in the field, None if the element is not in the field

        Raises:
            AssertionError: if the element is not in the field.

        Example:

        >>> from galwa import FiniteField
        >>> import numpy as np
        >>> p = 2
        >>> f = np.array([1, 1, 1])
        >>> F = FiniteField(p, f)
        >>> F.get_element(np.array([1, 1]))
        FiniteFieldElement(1 + x, f(x)= 1 + x + x², p=2)
        """
        assert FiniteFieldElement.check_that_element_is_in_field(a, self), "Element is not in the field"
        a = pad_element(a % self.p, self.f)
        return self._elements.get(tuple(a), None)

    def _calculate_illegal_powers(self) -> List[np.ndarray]:
        """
        Calculate the representation of :math:`x^i` for  :math:`n<=i<= 2n-2` as :math:`x^0, x^1, x^2, ..., x^{n-1}` terms.

        Methodology: since an element in the extension can have a degree of at most :math:`n-1`, and the highest degree\
        basis vector is :math:`n-1`, all "illegal" :math:`x` powers can be in\
         the range of :math:`n` to :math:`2(n-1) = 2n-2`.

        So we calculate each :math:`x^i` using an induction:

         :math:`x^i = x^{i-1} * x`

        Where :math:`x^{i-1}` is the previous :math:`x^i` and :math:`x` is the basis vector.

        So we start from :math:`x^n` and calculate :math:`x^(n+1)` using :math:`x^n` and :math:`x`, and so on until\
         we reach :math:`x^(2n-2)`.

        The multiplication of :math:`x^i` by :math:`x` can be seen as a shift of :math:`x^i` to the right,\
         so we shift :math:`x^i` to the right. Then we split the vector to a "valid" :math:`x^i` and\
        "illegal" :math:`x^i`, the valid :math:`x^i` are the :math:`x^i` which have a degree smaller then :math:`n`,
        and the illegal :math:`x^i` are the :math:`x^i` which have a degree bigger or equal to :math:`n`.

        The illegal :math:`x^i` will be represented as a sum of the valid :math:`x^i`.
        So we will be left with:


            some_valid_vector + illgeal_degree_coef * converted_illegal_vector

        Returns:
            List[np.ndarray]: a list of the representation of x^i for i in range n to 2n-2
        """
        n = len(self.f)
        r = n - 1
        x_max_pow_value = -self.f[:-1].copy()
        x_pow_list = [-self.f[:-1]]
        for i in range(r + 2, 2 * r):
            x_i = np.concatenate(([0], x_pow_list[-1]))
            x_i = x_i[-1] * x_max_pow_value + x_i[:-1]
            x_pow_list.append(x_i)
        return x_pow_list