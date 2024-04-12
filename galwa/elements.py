from typing import Union, Optional
import numpy as np

from galwa.utils import valid_repr, refactor_polynom_terms, pad_element, zero_element_check, same_field, invert_matrix
from galwa.utils import xgcd, same_prime_field


class PrimeFieldElement:
    """
    PrimeFieldElement class represents an element in a prime field :math:`F_{p}` where p is a prime number.
    The element represent the value a mod p

    Example:

    >>> from galwa import PrimeFieldElement
    >>> a = PrimeFieldElement(3, 5)
    >>> a
    PrimeFieldElement(value= 3,prime= 5)
    >>> print(a)
    3 mod 5
    >>> b = PrimeFieldElement(4, 5)
    >>> a + b
    PrimeFieldElement(value= 2,prime= 5)
    >>> a - b
    PrimeFieldElement(value= 4,prime= 5)
    >>> a * b
    PrimeFieldElement(value= 2,prime= 5)
    >>> a / b
    PrimeFieldElement(value= 2,prime= 5)
    >>> a ** 2
    PrimeFieldElement(value= 4,prime= 5)
    >>> a**-1
    PrimeFieldElement(value= 2,prime= 5)
    >>> a.inverse
    PrimeFieldElement(value= 2,prime= 5)
    >>> a == b
    False
    """
    def __init__(self, a: int, p: int):
        """
        Initialize the PrimeFieldElement class.

        Args:
            a (int): the element in the prime field.
            p (int): the prime number for the prime field.
        """
        self.a = a % p
        self.p = p
        self.i = self._find_inverse()

    def __str__(self) -> str:
        """
        String builtin function for PrimeFieldElement.

        Returns:
            str: the string representation of the element.
        """
        return f"{self.a} mod {self.p}"

    def __repr__(self) -> str:
        return f"PrimeFieldElement(value= {self.a},prime= {self.p})"

    @same_prime_field
    def __eq__(self, other):
        """
        Check if two prime field elements are equal.

        Args:
            other (PrimeFieldElement): the element to compare to.

        Returns:
            bool: True if the elements are equal, False otherwise.

        Raises:
            AssertionError: if the elements are not in the same prime field.
        """
        return self.a == other.a and self.p == other.p

    @same_prime_field
    def __add__(self, other):
        """
        Addition of two prime field elements.

        Args:
            other (PrimeFieldElement): the element to add.

        Returns:
            PrimeFieldElement: the result of the addition.

        Raises:
            AssertionError: if the elements are not in the same prime field.
        """
        return self.__class__((self.a + other.a) % self.p, self.p)

    @same_prime_field
    def __sub__(self, other):
        """
        Subtraction of two prime field elements.

        Args:
            other (PrimeFieldElement): the element to subtract.

        Returns:
            PrimeFieldElement: the result of the subtraction.

        Raises:
            AssertionError: if the elements are not in the same prime field.
        """
        return self.__class__((self.a - other.a) % self.p, self.p)

    @same_prime_field
    def __mul__(self, other):
        """
        Multiplication of two prime field elements.

        Args:
            other (PrimeFieldElement): the element to multiply by.

        Returns:
            PrimeFieldElement: the result of the multiplication.

        Raises:
            AssertionError: if the elements are not in the same prime field.
        """
        return self.__class__((self.a * other.a) % self.p, self.p)

    @same_prime_field
    def __truediv__(self, other):
        """
        Division of two prime field elements.

        Args:
            other (PrimeFieldElement): the element to divide by.

        Returns:
            PrimeFieldElement: the result of the division.

        Raises:
            AssertionError: if the elements are not in the same prime field.
        """
        if other.a == 0:
            raise ValueError("Division by zero")
        return self.__class__((self.a * other.i) % self.p, self.p)

    def __pow__(self, n):
        """
        Exponentiation of a prime field element.

        .. note::
          For this class we are NOT using exponentiation by squaring, since the built-in pow function is faster.

        Args:
            n (int): the exponent.

        Returns:
            PrimeFieldElement: the result of the exponentiation.
        """
        if n == 0:
            return self.__class__(1, self.p)
        if n < 0:
            return self.__class__(pow(self.i, -n, self.p), self.p)
        return self.__class__(pow(self.a, n, self.p), self.p)

    @property
    def inverse(self):
        """
        Returns the inverse of the element. (for multiplicative group)

        The inverse is being defined as the element b such that :math:`a*b = 1 mod p`.

        Returns:
            PrimeFieldElement: the inverse of the element.

        Example:

            >>> from galwa import PrimeFieldElement
            >>> a = PrimeFieldElement(3, 5)
            >>> a.inverse
            PrimeFieldElement(value= 2,prime= 5)
        """
        if self.a != 0:
            return PrimeFieldElement(self.i, self.p)
        else:
            return None

    def _find_inverse(self):
        """
        Finds the inverse of the element. (for multiplicative group)


        The inverse of an element a is the element b such that :math:`a*b = 1 mod p`, to find b we are using \
         the extended Euclidean algorithm.

        Returns:
            int: the inverse of the element.

        Raises:
            ValueError: if the element does not have an inverse.
        """
        if self.a == 0:
            return 0
        d, s, t = xgcd(self.a, self.p)
        if d != 1:
            raise ValueError("Element does not have an inverse")
        return s % self.p


class FiniteFieldElement(object):
    """
    FiniteFieldElement class represents an element in an extension :math:`F_{p^n}` where :math:`p` is a prime number \
    and :math:`n` is the degree of polynomial :math:`f(x)` used for the extension.
    In other words, a class to represent an element :math:`a \in l = F_{p^n}[x]/<f(x)>`

    Example:

    >>> from galwa import FiniteField, FiniteFieldElement
    >>> import numpy as np
    >>> f = np.array([1, 1, 0, 1])
    >>> p = 2
    >>> field = FiniteField(p, f)
    >>> a = FiniteFieldElement(np.array([1, 0, 1]), field)
    >>> a
    FiniteFieldElement(1 + x², f(x)= 1 + x + x³, p=2)
    """

    def __init__(self, a: np.ndarray, field: "FiniteField", representation: Optional[str] = "polynomial"):
        """
        Initialize the FiniteFieldElement class.

        Args:
            a (np.ndarray): the element in the field.
            field (FiniteField): the field the element belongs to.
            representation Optional[str]: the representation of the element (polynomial, vector, matrix), default is polynomial.

        Raises:
            AssertionError: if the representation is invalid or the element is not in the field or wrong type.
        """
        assert isinstance(a, np.ndarray), "Element must be a numpy array"
        assert representation in ["polynomial", "vector", "matrix"], "Invalid representation"
        assert self.check_that_element_is_in_field(a, field), "Element is not in the field"
        self.field = field
        self.p = field.p
        self.a = pad_element(a % field.p, field.f)
        self.gln_a = self.embed_in_gln()
        self.poly_a = np.polynomial.Polynomial(self.a)
        self.ord = None
        self._repr = representation

    def __eq__(self, other):
        """
        Check if two finite field elements are equal.

        Args:
            other (FiniteFieldElement): the element to compare to.

        Returns:
            bool: True if the elements are equal, False otherwise.
        """
        return np.array_equal(self.a, other.a) and self.field == other.field

    def __str__(self):
        """
        String builtin function for FiniteFieldElement.

        Returns:
            str: the string representation of the element by the representation type (polynomial, vector, matrix)
        """
        if self._repr == "polynomial":
            return refactor_polynom_terms(f"{self.poly_a}")
        elif self._repr == "vector":
            return self.a.__str__()
        else:
            return self.gln_a.__str__()

    def __repr__(self) -> str:
        """
        Representation builtin function for FiniteFieldElement.

        Returns:
            str: the representation of the element by the representation type (polynomial, vector, matrix)
        """
        if self._repr == "polynomial":
            a_poly_repr = refactor_polynom_terms(f"{self.poly_a}")
            f_poly_repr = refactor_polynom_terms(f"{self.field.poly_f}")
            return f"{self.__class__.__name__}({a_poly_repr}, f(x)= {f_poly_repr}, p={self.field.p})"
        elif self._repr == "vector":
            return f"{self.__class__.__name__}({self.a}, f(x) = {self.field.f} p={self.field.p})"
        else:
            return f"{self.__class__.__name__}({self.gln_a}, f(x) = {self.field.f} p={self.field.p})"

    @property
    def representation(self):
        """
        Get the representation of the element (polynomial, vector, matrix)

        Returns:
            str: the representation of the element.

        Example:

            >>> from galwa import FiniteField, FiniteFieldElement
            >>> import numpy as np
            >>> f = np.array([1, 1, 0, 1])
            >>> p = 2
            >>> field = FiniteField(p, f)
            >>> a = FiniteFieldElement(np.array([1, 0, 1]), field)
            >>> a.representation
            'polynomial'
        """
        return self._repr

    @representation.setter
    @valid_repr
    def representation(self, value) -> None:
        """
        Set the representation of the element.
        Valid representations are: "polynomial", "vector", "matrix".

        Args:
            value (str): the representation to set.

        Returns:
            None

        Raises:
            ValueError: if the representation is invalid.
        """
        self._repr = value

    @property
    def dimension(self) -> int:
        """
        Get the dimension of the element. (vector dimension)

        Returns:
            int: the dimension of the element
        """
        return len(self.a)

    @same_field
    def __add__(self, other) -> "FiniteFieldElement":
        """
        Calculate the addition of the element by another element, the addition is simple vector addition modulo p

        .. note::
         We dont use PrimeFieldElement for the addition since it is much easier to make vector addition and take \n
          the result modulo p. The result of addition of two integer will always be integer.

        Args:
            other (FiniteFieldElement): the element to add.

        Returns:
            FiniteFieldElement: the result of the addition.

        Raises:
            AssertionError: if the elements are not in the same field.

        Example:

            >>> from galwa import FiniteField, FiniteFieldElement
            >>> import numpy as np
            >>> f = np.array([1, 1, 0, 1])
            >>> p = 2
            >>> field = FiniteField(p, f)
            >>> a = FiniteFieldElement(np.array([1, 0, 1]), field)
            >>> b = FiniteFieldElement(np.array([1, 1, 0]), field)
            >>> a + b
            FiniteFieldElement(x + x², f(x)= 1 + x + x³, p=2)
        """
        return self.__class__((self.a + other.a) % self.field.p, self.field, representation=self._repr)

    @same_field
    def __sub__(self, other) -> "FiniteFieldElement":
        """
        Calculate the subtraction of the element by another element, the subtraction is simple vector subtraction modulo p

        .. note::
            We dont use PrimeFieldElement for the subtraction since it is much easier to make vector subtraction and take \n
            the result modulo p. The result of subtraction of two integer will always be integer.

        Args:
            other (FiniteFieldElement): the element to subtract by.

        Returns:
            FiniteFieldElement: the result of the subtraction.

        Raises:
            AssertionError: if the elements are not in the same field.

        Example:

            >>> from galwa import FiniteField, FiniteFieldElement
            >>> import numpy as np
            >>> f = np.array([1, 1, 0, 1])
            >>> p = 2
            >>> field = FiniteField(p, f)
            >>> a = FiniteFieldElement(np.array([1, 0, 1]), field)
            >>> b = FiniteFieldElement(np.array([1, 1, 0]), field)
            >>> a - b
            FiniteFieldElement(x + x², f(x)= 1 + x + x³, p=2)
        """
        return self.__class__((self.a - other.a) % self.field.p, self.field, representation=self._repr)

    @same_field
    @zero_element_check
    def __mul__(self, other) -> "FiniteFieldElement":
        """
        Calculate the multiplication of the element by another element
        First we check that either the element or the other element is not 0, since 0 is not part of the \
        multiplicative group so the multiplication is not defined.
        Then we calculate the multiplication of the element and the other element by multiplying the gln_a of the element \
        and the gln_a of the other element and take the result modulo p, the first column vector is the element vector.

        .. note::
         We dont use PrimeFieldElement for the multiplication since it is much easier to make matrix multiplication and take \n
            the result modulo p. The result of multiplication of two integer will always be integer.

        Args:
            other (FiniteFieldElement): the element to multiply by.

        Returns:
            FiniteFieldElement: the result of the multiplication.

        Raises:
            AssertionError: if the elements are not in the same field. or if the element is 0

        Example:

            >>> from galwa import FiniteField, FiniteFieldElement
            >>> import numpy as np
            >>> f = np.array([1, 1, 0, 1])
            >>> p = 2
            >>> field = FiniteField(p, f)
            >>> a = FiniteFieldElement(np.array([1, 0, 1]), field)
            >>> b = FiniteFieldElement(np.array([1, 1, 0]), field)
            >>> a * b
            FiniteFieldElement(x², f(x)= 1 + x + x³, p=2)
        """
        return self.__class__(np.matmul(self.gln_a, other.gln_a)[:, 0] % self.field.p, self.field, representation=self._repr)

    @same_field
    @zero_element_check
    def __truediv__(self, other) -> "FiniteFieldElement":
        """
        Calculate the division of the element by another element.
        First we check that either the element or the other element is not 0, since 0 is not part of the \
        multiplicative group so division is not defined.
        Then we calculate the inverse of the other element and multiply it by the element, the inverse is used by \
        inverting the gln_a of the other element and multiply it by the gln_a of the element.
        the first column vector is the element vector.

        .. note::
        This is the only operation that we use PrimeFieldElement, since regular matrix inverse can return float numbers \
        We must use PrimeFieldElement to perform the division between two elements.

        Args:
            other (FiniteFieldElement): the element to divide by.

        Returns:
            FiniteFieldElement: the result of the division.

        Raises:
            AssertionError: if the elements are not in the same field. or if the element is 0.

        Example:

            >>> from galwa import FiniteField, FiniteFieldElement
            >>> import numpy as np
            >>> f = np.array([1, 1, 0, 1])
            >>> p = 2
            >>> field = FiniteField(p, f)
            >>> a = FiniteFieldElement(np.array([1, 0, 1]), field)
            >>> b = FiniteFieldElement(np.array([1, 1, 0]), field)
            >>> a / b
            FiniteFieldElement(1 + x, f(x)= 1 + x + x³, p=2)
        """
        inverse = other._get_inverse()
        return self * inverse

    def __pow__(self, power) -> "FiniteFieldElement":
        """
        Built-in function for exponentiation of the element by a power
        For 0 power, the result is the identity element 1.

        If its negative we first get the inverse of the element \
        and then calculate the exponentiation by squaring of the inverse element.

        For positive power, we calculate the exponentiation by squaring of the element.

        Args:
            power (int): the power to exponentiate the element by.

        Returns:
            FiniteFieldElement: the result of the exponentiation.

        Example:

                >>> from galwa import FiniteField, FiniteFieldElement
                >>> import numpy as np
                >>> f = np.array([1, 1, 0, 1])
                >>> p = 2
                >>> field = FiniteField(p, f)
                >>> a = FiniteFieldElement(np.array([1, 0, 1]), field)
                >>> a ** 2
                FiniteFieldElement(1 + x + x², f(x)= 1 + x + x³, p=2)
        """
        if power == 0:
            return self.__class__(np.array([1] + [0] * (self.dimension - 1)), self.field, representation=self._repr)
        if power > 0:
            res = self._exponentiation_by_squaring(self, power)
            return res
        else:
            inverse = self._get_inverse()
            res = self._exponentiation_by_squaring(inverse, -power)
            return res

    def __hash__(self):
        """
        Hash builtin function, so we hash an element in a set or a dictionary.

        Returns:
            hash: the hash of the element
        """
        return hash(tuple(self.a))
 
    
    @staticmethod
    def _exponentiation_by_squaring(a, n) -> "FiniteFieldElement":
        """
        This function calculate the exponentiation of the element :math:`a^n` using the exponentiation by squaring algorithm.

        Args:
            a (FiniteFieldElement): the element to exponentiate.
            n (int): the power to exponentiate the element by.

        Returns:
            FiniteFieldElement: the result of the exponentiation.
        """

        if n == 0:
            return 1
        if n == 1:
            if a.is_identity_of_multiplication():
                return a
            return a
        res = FiniteFieldElement._exponentiation_by_squaring(a, n // 2)
        if n % 2 == 0:
            res = res * res
        else:
            res = a * res * res

        return res

    def _get_inverse(self) -> "FiniteFieldElement":
        """
        Get the inverse of the element (for multiplicative group)
        The inverse of the element is the element that when multiplied by the element gives the identity element 1

        .. note::
            This is the only operation that we use PrimeFieldElement, since regular matrix inverse can return float numbers. \
            We must use PrimeFieldElement to perform the division between two elements.

        Returns:
            FiniteFieldElement: the inverse of the element
        """
        gln_a_list = [[PrimeFieldElement(value, self.field.p) for value in row] for row in self.gln_a]
        inverse_a = np.array(invert_matrix(gln_a_list))[:,0]
        inverse_a = np.array([element.a for element in inverse_a])
        return self.__class__(inverse_a, self.field, representation=self._repr)

    def as_vector(self) -> None:
        """
        Change sthe representation of the element to vector.
        So once printing the element, it will be printed as a vector.

        Returns:
            None

        Example:

            >>> from galwa import FiniteField, FiniteFieldElement
            >>> import numpy as np
            >>> f = np.array([1, 1, 0, 1])
            >>> p = 2
            >>> field = FiniteField(p, f)
            >>> a = FiniteFieldElement(np.array([1, 0, 1]), field)
            >>> a.as_vector()
            >>> a
            FiniteFieldElement([1 0 1], f(x) = [1 1 0 1] p=2)
            >>> print(a)
            [1 0 1]
        """
        self.representation = "vector"

    def as_matrix(self) -> None:
        """
        Changes the representation of the element to matrix.
        So once printing the element, it will be printed as a matrix.

        Returns:
            None

        Example:

            >>> from galwa import FiniteField, FiniteFieldElement
            >>> import numpy as np
            >>> f = np.array([1, 1, 0, 1])
            >>> p = 2
            >>> field = FiniteField(p, f)
            >>> a = FiniteFieldElement(np.array([1, 0, 1]), field)
            >>> a.as_matrix()
            >>> a
            FiniteFieldElement([[1 0 1]
                                [1 0 0]
                                [0 1 0]], f(x) = [1 1 0 1] p=2)
            >>> print(a)
            [[1 0 1]
            [1 0 0]
            [0 1 0]]
        """
        self.representation = "matrix"

    def as_polynomial(self) -> None:
        """
        Changes the representation of the element to polynomial.
        So once printing the element, it will be printed as a polynomial.

        Returns:
            None

        Example:

            >>> from galwa import FiniteField, FiniteFieldElement
            >>> import numpy as np
            >>> f = np.array([1, 1, 0, 1])
            >>> p = 2
            >>> field = FiniteField(p, f)
            >>> a = FiniteFieldElement(np.array([1, 0, 1]), field)
            >>> a.as_polynomial()
            >>> a
            FiniteFieldElement(1 + x², f(x)= 1 + x + x³, p=2)
            >>> print(a)
            1 + x²
        """
        self.representation = "polynomial"

    @property
    def order(self) -> Union[int, None]:
        """
        Returns the multiplicative order of the element.
        The order is define as

        .. math::
            O(a) = min(n > 0 : a^n = 1 mod p, n \in N)

        .. note::
            At first call, the order is calculated and stored in the ord attribute.
            The reason for that is that the calculation of the order can be an expensive operation.
            First call might take some time depending on the value of p, but once calculated the order can be accessed \
            quickly.

        Returns:
            int: the multiplicative order of the element. None if the element is 0
        """
        if self.ord is None:
            self.ord = self.multiplicative_order()
        return self.ord

    def multiplicative_order(self) -> Union[int, None]:
        """
        Calculates the multiplicative order of the element in the field.

        The multiplicative order of an element a in a finite field is the smallest positive integer n such that\
         :math:`a^n = 1`
        For the 0 element, the order is not defined since 0 is not part of the multiplicative group.

        Methodology:

        From lagrange theorem we know that for :math:`H` as subgroup of :math:`G` then the order of any element in :math:`G` divides the order of :math:`G`.\

        That is true for all :math:`g \in G , |<g>| | |G|`

        In our case the multiplicative group of :math:`F_{p}^x = F_{p} - \{0\}` is a subgroup of the multiplicative group :math:`l^x`.

        So for  all :math:`a \in l^x, O(a) | O(l^x)`.

        So first, we calculate all divisors of the order of the multiplicative group of the field, the complexity of this operation is :math:`O(\sqrt{n})`\ where :math:`n` is the order of the multiplicative group.

        The divisors array will be sorted, we will start from the smallest and calculate :math:`a^{d}` for each divisor :math:`d` and check if the result is the identity element.

        Calculating :math:`a^{d}` can be done using exponentiation by squaring algorithm, the complexity of this operation is :math:`O(\log(d))`.

        So the complexity in the best case will be :math:`O(\log(d)` and in the worst case :math:`O(k\log(d))` where :math:`k` is the number of divisors.


        Returns:
            int: the multiplicative order of the element, None if the element is 0.

        Example:

            >>> from galwa import FiniteField, FiniteFieldElement
            >>> import numpy as np
            >>> f = np.array([1, 1, 0, 1])
            >>> p = 2
            >>> field = FiniteField(p, f)
            >>> a = FiniteFieldElement(np.array([1, 0, 1]), field)
            >>> a.multiplicative_order()
            7
        """
        if self.ord is not None:
            return self.ord
        if self.gln_a is None or np.all(self.a == 0):
            return None
        for divider in self.field.field_size_dividers:
            res = self._exponentiation_by_squaring(self, divider)

            if res.is_identity_of_multiplication():
                return divider

        return None

    def is_generator(self):
        """
        Check if the element is a generator in the field.

        A generator element is an element whose order is equal to the order of the multiplicative group.
        In other words, a generator element is an element whose powers generate all the elements in the field.

        Returns:
            bool: True if the element is a generator, False otherwise

        Example:
            >>> from galwa import FiniteField, FiniteFieldElement
            >>> import numpy as np
            >>> f = np.array([1, 1, 0, 1])
            >>> p = 2
            >>> field = FiniteField(p, f)
            >>> a = FiniteFieldElement(np.array([1, 0, 1]), field)
            >>> a.is_generator()
            True
        """
        return self.order == (self.field.order - 1)

    def is_identity_of_multiplication(self):
        """
        Checks if the element is the identity element of the multiplication operation.
        The identity element of the multiplication operation is the element 1.

        Returns:
            bool: True if the element is the identity element, False otherwise

        Example:
            >>> from galwa import FiniteField, FiniteFieldElement
            >>> import numpy as np
            >>> f = np.array([1, 1, 0, 1])
            >>> p = 2
            >>> field = FiniteField(p, f)
            >>> a = FiniteFieldElement(np.array([1, 0, 1]), field)
            >>> a.is_identity_of_multiplication()
            False
        """
        return np.array_equal(self.a, np.array([1] + [0] * (self.dimension - 1)))

    @staticmethod
    def check_that_element_is_in_field(element: np.ndarray, field: "FiniteField") -> bool:
        """
        Check if the element is in the field with the given irreducible polynomial :math:`f(x)`

        A valid element must be of the form :math:`a = a_{0} + a_{1}x + a_{2}x^2 + ... + a_{n-1}x^{n-1}` where\
         :math:`n` is the degree of :math:`f(x)`.
        In other words, the degree of the element must be less than the degree of :math:`f(x)`

        Args:
            element (np.ndarray): the element to check.
            field (FiniteField): the field to check the element in.

        Returns:
            bool: True if the element is in the field, False otherwise.

        Example:
            >>> from galwa import FiniteField, FiniteFieldElement
            >>> import numpy as np
            >>> f = np.array([1, 1, 0, 1])
            >>> p = 2
            >>> field = FiniteField(p, f)
            >>> a = np.array([1, 0, 1])
            >>> FiniteFieldElement.check_that_element_is_in_field(a, field)

        """
        return len(field.f) > len(element)

    def embed_in_gln(self) -> Union[np.ndarray, None]:
        """
        Embed the element in :math:`GL_{n}(p)` where :math:`n` is the degree of the irreducible polynom :math:`f(x)`.

        Methodology: for an element a in the field, we calculate the multiplication of :math:`a*x^i` for :math:`0<=i<=(n-1)`
        and store the results in the columns of a matrix.

        If the multiplication has a degree higher than :math:`n-1`, we use previously calculated vectors\
        which represent :math:`x^i` for  :math:`n <= i <= 2n-2`.

        The maximum degree of the multiplication is at most 2n-2 since the highest degree in an element is n-1 and the highest degree basis
        vector is also n-1. So the highest degree in the multiplication is :math:`2(n-1) = 2n-2`

        .. note::
            For the zero element, there is no representation in :math:`GL_{n}(p)` since the zero element is not part of the multiplicative group.
            In this case we return None

        Returns:
            np.ndarray: the matrix representation of the element in :math:`GL_{n}(p)`

        Example:
            >>> from galwa import FiniteField, FiniteFieldElement
            >>> import numpy as np
            >>> f = np.array([1, 1, 0, 1])
            >>> p = 2
            >>> field = FiniteField(p, f)
            >>> a = FiniteFieldElement(np.array([1, 0, 1]), field)
            >>> a.embed_in_gln()
            array([[1, 1, 0],
                   [0, 0, 1],
                   [1, 0, 0]])
        """
        if np.all(self.a == 0):
            return None
        n = len(self.field.f)
        r = n - 1
        x_pow_list = self.field.x_powers
        emmbed_matrix = np.zeros((r, r), dtype=int)
        for i in range(r):
            basis_vector = np.eye(1, r, i).flatten()
            a_mul_basis = np.convolve(self.a, basis_vector, mode='full')
            a_mul_basis_valid = a_mul_basis[:r].copy()
            for j in range(r, 2 * r - 1):
                x_j_co = a_mul_basis[j]
                x_j = x_pow_list[j - r]
                v = x_j_co * x_j
                a_mul_basis_valid += v
            emmbed_matrix[:, i] = a_mul_basis_valid
        return emmbed_matrix % self.field.p