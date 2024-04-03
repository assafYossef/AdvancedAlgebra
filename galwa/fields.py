from typing import List, Dict, Tuple
import numpy as np
from galwa.elements import FiniteFieldElement
from galwa.utils import valid_repr, refactor_polynom_terms, pad_element


class FiniteField(object):
    """
    FiniteField class represents an extension field of the form F_p^n where p is a prime number and n is the degree of \
    polynomial f(x) used for the extension.
    l = F_p^n[x]/<f(x)>
    """

    def __init__(self, p, f: np.ndarray, representation: str = "polynomial"):
        """
        Initialize the FiniteField class

        Args:
            p (int): the prime number for the field
            f (np.ndarray): the irreducible polynomial f(x) used for the extension
            representation (str): the representation of the elements in the field - polynomial, vector, matrix

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

    def __eq__(self, other) -> bool:
        return self.p == other.p and np.array_equal(self.f, other.f)

    def __repr__(self) -> str:
        """
        Representation builtin function for FiniteField
        Print the representation of the field - polynomial, vector

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
        """
        return self._repr

    @representation.setter
    @valid_repr
    def representation(self, value) -> None:
        """
        Sets the representation of the elements in the field - polynomial, vector, or matrix
        So once printing the elements, they will be printed as polynomials, vectors, or matrices

        Args:
            value (str): the representation to set (polynomial, vector, matrix)

        Returns:
            None
        """
        self._repr = value
        for element in self._elements.values():
            element.representation = value

    def elements_as_vectors(self) -> None:
        """
        Change the representation of the elements to vector
        So once printing the elements, they will be printed as vectors

        Returns:
            None
        """
        self.representation = "vector"

    def elements_as_matrices(self) -> None:
        """
        Change the representation of the elements to matrix
        So once printing the elements, they will be printed as matrices

        Returns:
            None
        """
        self.representation = "matrix"

    def elements_as_polynomials(self) -> None:
        """
        Change the representation of the elements to polynomial
        So once printing the elements, they will be printed as polynomials

        Returns:
            None
        """
        self.representation = "polynomial"

    def _check_that_f_is_irreducible(self) -> bool:
        """
        Check if the polynomial f(x) is irreducible over the field F_p for degree 2 or 3 \
        For bigger degrees we assume that the polynomial is irreducible.
        A polynomial f(x) of degree 2 or 3 is irreducible iff it has no roots in the field F_p

        Returns:
            bool: True if the polynomial is irreducible, False otherwise
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

        Remark: At first call the generators are calculated and stored in the generators attribute \
        since calculating the generators is expensive so only calculate upon request

        Returns:
            List[FiniteFieldElement]: a list of all the generators in the field
        """
        if self._generators is None:
            self._generators = self._find_generators()
        return self._generators

    @property
    def order(self) -> int:
        """
        Property to get the order of the field, which is p^n

        Returns:
            int: the order of the field
        """
        return self.p ** (len(self.f) - 1)

    def _find_generators(self):
        """
        Find all the generators in the field
        A generator element is an element whose order is equal to the order of the multiplicative group
        In other words, a generator element is an element whose powers generate all the elements in the field.

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
        Create all possible elements in the field
        The elements are of the form a_0 + a_1x + a_2x^2 + ... + a_{n-1}x^{n-1} where n is the degree of f(x)
        So we create a permutation of all possible coefficients for the elements in the field and use them to create the elements

        Returns:
            Dict[Tuple[np.ndarray], FiniteFieldElement]: a dictionary of all the elements in the field where the key is the vector representation of the element

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
        Returns all the elements in the field as a list

        Returns:
            List[FiniteFieldElement]: a list of all the elements in the field
        """
        return list(self._elements.values())

    def get_element(self, a: np.ndarray) -> FiniteFieldElement:
        """
        Given a vector representation of an element in the field, return the element

        Args:
            a (np.ndarray): the vector representation of the element

        Returns:
            FiniteFieldElement: the element in the field, None if the element is not in the field
        """
        assert FiniteFieldElement.check_that_element_is_in_field(a, self), "Element is not in the field"
        a = pad_element(a % self.p, self.f)
        return self._elements.get(tuple(a), None)

    def _calculate_illegal_powers(self) -> List[np.ndarray]:
        """
        Calculate the representation of x^i for i in range n to 2n-2 as x^0, x^1, x^2, ..., x^(n-1) terms.
        Methodology: since an element in the extension can have a degree of at most n-1, and the highest degree basis \
        vector is n-1, so all "illegal" x powers can be in the range of n to 2n-2.
        We calculate each x_i using a induction - x_i = x_i-1 * x, where x_i-1 is the previous x_i and x is the basis vector.
        so we start from x^n and calculate x^(n+1) using x^n and x, and so on until we reach x^(2n-2)
        The multiplication of x^i by x can seen as a shift of x^i to the right, so we shift x^i to the right.
        Then we split the vector to "valid" x^i and "illegal" x^i, the valid x^i are the x^i which have a degree < n,
        and the illegal x^i are the x^i which have a degree >= n. The illgeal x^i will be represented as a sum of the valid x^i
        so we will left with some_valid_vector + illgeal_degree_coef* converted_illegal_vector

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