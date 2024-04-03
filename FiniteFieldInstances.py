from typing import List, Dict, Tuple, Union
import numpy as np
from utils import valid_repr, refactor_polynom_terms, pad_element, zero_element_check, same_field, invert_matrix
from PrimeFieldElementClass import PrimeFieldElement

class FiniteField(object):
    """
    FiniteField class represents an extension field of the form F_p^n where p is a prime number and n is the degree of \
    polynomial f(x) used for the extension.
    l = F_p^n[x]/<f(x)>
    """

    def __init__(self, p, f: np.ndarray, representation: str = "polynomial"):
        assert representation in ["polynomial", "vector", "matrix"], "Invalid representation"
        self._repr = representation
        self.p = p
        self.f = f % p
        assert self._check_that_f_is_irreducible(), f"f(x) is reducible over F{p}"
        self.poly_f = np.polynomial.Polynomial(self.f)
        self.x_powers = self._calculate_illegal_powers()
        self._elements = self._create_elements()
        self._generators = None

    def __eq__(self, other):
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
    def generators(self):
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

    def _create_elements(self):
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
    def elements(self):
        """
        Returns all the elements in the field as a list

        Returns:
            List[FiniteFieldElement]: a list of all the elements in the field
        """
        return list(self._elements.values())

    def get_element(self, a: np.ndarray):
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


class FiniteFieldElement(object):
    """
    FiniteFieldElement class represents an element in an extension field of the form F_p^n where p is a prime number \
     and n is the degree of polynomial f(x) used for the extension.
     element that belongs to l = F_p^n[x]/<f(x)>
    """

    def __init__(self, a: np.ndarray, field, representation: str = "polynomial"):
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
        if self.gln_a is None:
            return None
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
