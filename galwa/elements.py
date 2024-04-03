from typing import Union
import numpy as np
from galwa.utils import valid_repr, refactor_polynom_terms, pad_element, zero_element_check, same_field, invert_matrix
from galwa.utils import xgcd, same_prime_field


class PrimeFieldElement:
    """
    PrimeFieldElement class represents an element in a prime field of the form F_p where p is a prime number.
    The element represent the value a mod p
    """
    def __init__(self, a: int, p: int):
        """
        Initialize the PrimeFieldElement class

        Args:
            a (int): the element in the prime field
            p (int): the prime number for the prime field
        """
        self.a = a % p
        self.p = p
        self.i = self._find_inverse()

    def __str__(self):
        """
        String builtin function for PrimeFieldElement

        Returns:
            str: the string representation of the element

        Example:
        PrimeFieldElement(3, 5)
        """
        return f"{self.a} mod {self.p}"

    def __repr__(self):
        return f"PrimeFieldElement(value= {self.a},prime= {self.p})"

    @same_prime_field
    def __eq__(self, other):
        return self.a == other.a and self.p == other.p and self.inverse == other.inverse

    @same_prime_field
    def __add__(self, other):
        """
        Addition of two prime field elements
        Args:
            other (PrimeFieldElement): the element to add
        Returns:
            PrimeFieldElement: the result of the addition
        """
        return self.__class__((self.a + other.a) % self.p, self.p)

    @same_prime_field
    def __sub__(self, other):
        """
        Subtraction of two prime field elements
        Args:
            other (PrimeFieldElement): the element to subtract
        Returns:
            PrimeFieldElement: the result of the subtraction
        """
        return self.__class__((self.a - other.a) % self.p, self.p)

    @same_prime_field
    def __mul__(self, other):
        """
        Multiplication of two prime field elements
        Args:
            other (PrimeFieldElement): the element to multiply by
        Returns:
            PrimeFieldElement: the result of the multiplication
        """
        return self.__class__((self.a * other.a) % self.p, self.p)

    @same_prime_field
    def __truediv__(self, other):
        """
        Division of two prime field elements
        Args:
            other (PrimeFieldElement): the element to divide by
        Returns:
            PrimeFieldElement: the result of the division
        """
        if other.a == 0:
            raise ValueError("Division by zero")
        return self.__class__((self.a * other.i) % self.p, self.p)

    def __pow__(self, n):
        """
        Exponentiation of a prime field element
        Args:
            n (int): the exponent
        Returns:
            PrimeFieldElement: the result of the exponentiation
        """
        if n == 0:
            return self.__class__(1, self.p)
        if n < 0:
            return self.__class__(pow(self.i, -n, self.p), self.p)
        return self.__class__(pow(self.a, n, self.p), self.p)

    @property
    def inverse(self):
        """
        Returns the unit of the element
        Returns:
            int: the unit of the element
        """
        return PrimeFieldElement(self.i, self.p)

    def _find_inverse(self):
        """
        find the unit of the element
        The unit of an element a is the element b such that a*b = 1 mod p, to find b we use the extended euclidean algorithm
        Returns:
            int: the unit of the element

        Raises:
            ValueError: if the element does not have an inverse
        """
        if self.a == 0:
            return 0
        d, s, t = xgcd(self.a, self.p)
        if d != 1:
            raise ValueError("Element does not have an inverse")
        return s % self.p


class FiniteFieldElement(object):
    """
    FiniteFieldElement class represents an element in an extension field of the form F_p^n where p is a prime number \
     and n is the degree of polynomial f(x) used for the extension.
     element that belongs to l = F_p^n[x]/<f(x)>
    """

    def __init__(self, a: np.ndarray, field, representation: str = "polynomial"):
        """
        Initialize the FiniteFieldElement class

        Args:
            a (np.ndarray): the element in the field
            field (FiniteField): the field the element belongs to
            representation (str): the representation of the element - polynomial, vector, matrix

        Raises:
            AssertionError: if the representation is invalid or the element is not in the field
        """
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
        return np.array_equal(self.a, other.a) and self.field == other.field

    def __str__(self):
        """
        String builtin function for FiniteFieldElement
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
        Representation builtin function for FiniteFieldElement
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
            str: the representation of the element
        """
        return self._repr

    @representation.setter
    @valid_repr
    def representation(self, value) -> None:
        """
        Set the representation of the element
        Valid representations are: polynomial, vector, matrix
        Args:
            value (str): the representation to set

        Returns:
            None
        """
        self._repr = value

    @property
    def dimension(self) -> int:
        """
        Get the dimension of the element

        Returns:
            int: the dimension of the element
        """
        return len(self.a)

    @same_field
    def __add__(self, other):
        """
        Calculate the addition of the element by another element, the addition is simple vector addition modulo p
        Args:
            other (FiniteFieldElement): the element to add

        Returns:
            FiniteFieldElement: the result of the addition
        """
        return self.__class__((self.a + other.a) % self.field.p, self.field, representation=self._repr)

    @same_field
    def __sub__(self, other):
        """
        Calculate the subtraction of the element by another element, the subtraction is simple vector subtraction modulo p
        Args:
            other (FiniteFieldElement): the element to subtract by

        Returns:
            FiniteFieldElement: the result of the subtraction
        """
        return self.__class__((self.a - other.a) % self.field.p, self.field, representation=self._repr)

    @same_field
    @zero_element_check
    def __mul__(self, other):
        """
        Calculate the multiplication of the element by another element
        First we check that either the element or the other element is not 0, since 0 is not part of the \
        multiplicative group so multiplication is not defined.
        Then we calculate the multiplication of the element and the other element by multiplying the gln_a of the element \
        and the gln_a of the other element and take the result modulo p, the first column vector is the element vector.

        Args:
            other (FiniteFieldElement): the element to multiply by

        Returns:
            FiniteFieldElement: the result of the multiplication
        """
        return self.__class__(np.matmul(self.gln_a, other.gln_a)[:, 0] % self.field.p, self.field, representation=self._repr)

    @same_field
    @zero_element_check
    def __truediv__(self, other):
        """
        Calculate the division of the element by another element
        First we check that either the element or the other element is not 0, since 0 is not part of the \
        multiplicative group so division is not defined.
        Then we calculate the inverse of the other element and multiply it by the element, the inverse is used by \
        inverting the gln_a of the other element and multiply it by the gln_a of the element and take the result modulo p.
        the first column vector is the element vector.
        Args:
            other (FiniteFieldElement): the element to divide by

        Returns:
            FiniteFieldElement: the result of the division
        """
        inverse = other._get_inverse()
        return self * inverse

    def __pow__(self, power):
        """
        Built-in function for exponentiation of the element by a power
        For 0 power, the result is the identity element 1, if its negative we first get the inverse of the element \
        and then calculate the exponentiation by squaring of the inverse element.
        For positive power, we calculate the exponentiation by squaring of the element.
        :param power:
        :return:
        """
        if power == 0:
            return self.__class__(np.array([1] + [0] * (self.dimension - 1)), self.field, representation=self._repr)
        if power > 0:
            res, element_order = self._exponentiation_by_squaring_with_order(self, power)
            if self.ord is None:
                self.ord = element_order
            return res
        else:
            inverse = self._get_inverse()
            res, element_order = self._exponentiation_by_squaring_with_order(inverse, -power)
            if self.ord is None:
                self.ord = element_order
            return res

    def __hash__(self):
        return hash(tuple(self.a))

    @staticmethod
    def _exponentiation_by_squaring_with_order(a, n):
        """
        This function have 2 purposes:
        first calculate the exponentiation of the element a by n using the exponentiation by squaring algorithm.
        secondly, we use the fact that if a^n = 1 then the order of a divides n, so we can try and calculate the order \
         of the element when calculating the exponentiation by squaring.
         This is true from lagrange theorem, that for H a subgroup of G, the order of any element in G divides \
            the order of the group G. For all g in G, |<g>| | |G|
        In our case the multiplicative group of Fp* = Fp - {0} is a subgroup of the multiplicative group l* where l \
        is the extended field l = Fp[x]/<f(x)>. So for any element a in l*, O(a) | O(l*)

        Args:
            a (FiniteFieldElement): the element to exponentiate
            n (int): the power to exponentiate the element by

        Returns:
            FiniteFieldElement: the result of the exponentiation
            int: the order of the element, None if the order couldn't be found for the specific n
        """
        if n == 0:
            return 1, None
        if n == 1:
            if a.is_identity_of_multiplication():
                return a, 1
            return a, None
        res, element_order = FiniteFieldElement._exponentiation_by_squaring_with_order(a, n // 2)
        if n % 2 == 0:
            res = res * res
        else:
            res = a * res * res
        if res.is_identity_of_multiplication() and element_order is None:
            element_order = n
        return res, element_order

    def _get_inverse(self):
        """
        Get the inverse of the element
        The inverse of the element is the element that when multiplied by the element gives the identity element 1

        Returns:
            FiniteFieldElement: the inverse of the element
        """
        gln_a_list = [[PrimeFieldElement(value, self.field.p) for value in row] for row in self.gln_a]
        inverse_a = np.array(invert_matrix(gln_a_list))[:,0]
        inverse_a = np.array([element.a for element in inverse_a])
        return self.__class__(inverse_a, self.field, representation=self._repr)

    def as_vector(self) -> None:
        """
        Change the representation of the element to vector
        So once printing the element, it will be printed as a vector

        Returns:
            None
        """
        self.representation = "vector"

    def as_matrix(self) -> None:
        """
        Change the representation of the element to matrix
        So once printing the element, it will be printed as a matrix

        Returns:
            None
        """
        self.representation = "matrix"

    def as_polynomial(self) -> None:
        """
        Change the representation of the element to polynomial
        So once printing the element, it will be printed as a polynomial

        Returns:
            None
        """
        self.representation = "polynomial"

    @property
    def order(self) -> int:
        """
        Calculate the multiplicative order of the element in the field

        Remark: At first call the order is calculated and stored in the ord attribute, since calculating the order is expensive
        so only calculate upon request

        Returns:
            int: the multiplicative order of the element
        """
        if self.ord is None:
            self.ord = self.multiplicative_order()
        return self.ord

    def multiplicative_order(self) -> Union[int, None]:
        """
        Calculate the multiplicative order of the element in the field
        The multiplicative order of an element a in a finite field is the smallest positive integer n such that a^n = 1
        For the 0 element, the order is not defined since 0 is not part of the multiplicative group

        Returns:
            int: the multiplicative order of the element, None if the element is 0
        """
        if self.ord is not None:
            return self.ord
        if self.gln_a is None or np.all(self.a == 0):
            raise ValueError("0 element has no order, its not part of the multiplicative group")
        _, element_order = self._exponentiation_by_squaring_with_order(self, self.field.order - 1)
        return element_order

    def is_generator(self):
        """
        Check if the element is a generator in the field
        A generator element is an element whose order is equal to the order of the multiplicative group
        In other words, a generator element is an element whose powers generate all the elements in the field.

        Returns:
            bool: True if the element is a generator, False otherwise
        """
        return self.order == (self.field.order - 1)

    def is_identity_of_multiplication(self):
        """
        Check if the element is the identity element of the multiplication operation
        The identity element of the multiplication operation is the element 1

        Returns:
            bool: True if the element is the identity element, False otherwise
        """
        return np.array_equal(self.a, np.array([1] + [0] * (self.dimension - 1)))

    @staticmethod
    def check_that_element_is_in_field(element, field) -> bool:
        """
        Check if the element is in the field with the given irreducible polynomial f(x)
        A valid element must be of the form a = a_0 + a_1x + a_2x^2 + ... + a_{n-1}x^{n-1} where n is the degree of f(x)
        In other words, the degree of the element must be less than the degree of f(x)

        Args:
            element (np.ndarray): the element to check
            field (FiniteField): the field to check the element in

        Returns:
            bool: True if the element is in the field, False otherwise
        """
        return len(field.f) > len(element)

    def embed_in_gln(self) -> Union[np.ndarray, None]:
        """
        Embed the element in GL(n, p) where n is the degree of the extension field.
        Methodology: for an element a in the field, we calculate the multiplication of a by x^i for i in range(n-1)
        and put the results in the columns of a matrix. If the multiplication has a degree higher than n-1, we
        use previously calculated vectors which represent x^i  n <= i < 2n-2
        the max degree of the multiplication is at most 2n-2 since the highest degree in a is n-1 and the highest degree basis
        vector is n-1. So the highest degree in the multiplication is 2n-2

        Remark: for the 0 element, there is no representation in GL(n, p) since the 0 element is not part of the multiplicative group

        Returns:
            np.ndarray: the matrix representation of the element in GL(n, p)
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