import numpy as np
import re


def bsgs(generator, element, group_order):
    """
    Baby-step Giant-step algorithm to solve the discrete logarithm problem.

    .. math::
        g^x = h

    Args:
        generator (FiniteFieldElement): g in :math:`g^x = h`
        element (FiniteFieldElement): h in :math:`g^x = h`
        group_order (int): group order , to initialize the table.

    Returns:
        int: x such that g^x = h mod p

    Raises:
        ValueError: if the discrete logarithm is not found

    Example:
        >>> from galwa import FiniteFieldElement, FiniteField
        >>> from galwa.utils import bsgs
        >>> import numpy as np
        >>> f = np.array([2, 0, 0, 2, 1])
        >>> p = 3
        >>> F = FiniteField(p, f)
        >>> g = FiniteFieldElement(np.array([1, 1, 0, 0]), F)
        >>> h = g ** 10
        >>> h
        FiniteFieldElement(2, f(x)= 2 + 2·x³ + x⁴, p=3)
        >>> order = F.order - 1
        >>> bsgs(g, h, order)
        10
    """
    m = int(np.ceil(np.sqrt(group_order)))
    table = {}
    # Baby-step
    for j in range(m):
        table[pow(generator, j)] = j
    # Giant-step
    c = pow(generator, -m)
    for i in range(m):
        y = (element * pow(c, i))
        if y in table:
            return i * m + table[y]
    raise ValueError("Discrete logarithm not found")


def xgcd(a, b):
    """
    Extended Euclidean algorithm to find the greatest common divisor and the coefficients of Bezout's identity.

    Args:
        a (int): first number
        b (int): second number

    Returns:
        tuple: gcd, s, t such that :math:`gcd(a, b) = s * a + t * b`

    Example:

        >>> from galwa.utils import xgcd
        >>> xgcd(240, 46)
        (2, -9, 47)
    """
    if b == 0:
        if a < 0:
            return -a, -1, 0
        return a, 1, 0
    else:
        q, r = a // b, a % b
        d, s, t = xgcd(b, r)
        return d, t, s - q * t


def same_prime_field(func):
    """
    Helper function to check if the operands belong to the same prime field, or if one of them is an integer,\
     convert it to a prime field element.
    """

    def wrapper(self, other):
        if isinstance(other, int):
            other = self.__class__(other, self.p)
            return func(self, other)
        if self.p != other.p:
            raise ValueError("Operands must belong to the same prime field")
        return func(self, other)

    return wrapper


def same_field(func):
    """
    Helper function to check if the operands belong to the same field.
    """

    def wrapper(self, other):
        if self.field != other.field:
            raise ValueError("Operands must belong to the same field")
        return func(self, other)

    return wrapper


def zero_element_check(func):
    """
    Helper function to check if an element is zero for the multiplicative group, so multiplication or division cannot be done.
    """

    def wrapper(self, other):
        if np.all(self.a == 0) or np.all(other.a == 0):
            raise ValueError("Zero is not part of the multiplicative group, so it cannot be used in this operation")
        return func(self, other)

    return wrapper


def valid_repr(func):
    """
    Helper function to check if the representation is valid.

    A valid representation is either "polynomial", "vector" or "matrix"
    """

    def wrapper(self, value):
        assert value in ["polynomial", "vector", "matrix"], "Invalid representation"
        return func(self, value)

    return wrapper


def refactor_polynom_terms(poly_repr: str) -> str:
    """
    Helper function to refactor the polynomial representation.

    The polynomial representation of numpy returns float coefficients, as well it doesnt drop zero coefficients.

    This function will refactor the polynomial representation to a more human readable form

    Example:
    poly_repr = :math:`1.0 x^2 + 0.0 x^1 + 3.0 x^0`
    refactor_polynom_terms(poly_repr) -> :math:`x^2 + 3`

    Args:
        poly_repr (str): polynomial representation.

    Returns:
        str: refactored polynomial representation.
    """
    poly_repr = f"{poly_repr}".replace(".0", "")
    new_str = ""
    split_str = re.split(r'\s*([-+])\s*', poly_repr)

    for term in split_str:
        term = term.strip()
        if term[0] != "0":
            if term[0] == "1":
                if len(term) > 1:
                    new_str += f"{term[2:]} "
                else:
                    new_str += f"{term} "
            else:
                new_str += f"{term} "
        else:
            if len(new_str) > 0:
                new_str = new_str[:-2]
    new_str = new_str.strip("+- ")
    return new_str if new_str else "0"


def pad_element(element: np.array, f: np.array) -> np.array:
    """
    Helper function to pad the element with zeros to match the size of the irreducible polynomial.

    Args:
        element (np.array): element to pad.
        f (np.array): irreducible polynomial to match the size.

    Returns:
        np.array: padded element (if needed)
    """
    return np.pad(element, pad_width=(0, len(f) - len(element) - 1), mode='constant')


def determinant(matrix: list) -> int:
    """
    Calculates the determinant of a square matrix recursively.

    Args:
        matrix (list): a square matrix represented as a list of lists.

    Returns:
        the determinant of the matrix.
    """
    n = len(matrix)
    if n == 1:
        return matrix[0][0]
    elif n == 2:
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    else:
        det = 0
        for j in range(n):
            minor = [row[:j] + row[j + 1:] for row in matrix[1:]]
            det = matrix[0][j] * determinant(minor) * ((-1) ** j) + det
        return det


def _transpose(matrix: list) -> list:
    """
    Transposes a matrix.

    Args:
        matrix (list): a matrix represented as a list of lists.

    Returns:
        list: the transposed matrix.
    """
    return [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]


def _cofactor_matrix(matrix: list) -> list:
    """
    Calculates the cofactor matrix of a square matrix.

    Args:
        matrix (list): a square matrix represented as a list of lists.

    Returns:
        list: the cofactor matrix.
    """
    n = len(matrix)
    cofactors = []
    for i in range(n):
        cofactor_row = []
        for j in range(n):
            minor = [row[:j] + row[j + 1:] for row in (matrix[:i] + matrix[i + 1:])]
            cofactor_row.append(determinant(minor) * ((-1) ** (i + j)))
        cofactors.append(cofactor_row)
    return cofactors


def _adjugate_matrix(matrix: list) -> list:
    """
    Calculates the adjugate matrix of a square matrix.

    Args:
        matrix (list): a square matrix represented as a list of lists.

    Returns:
        list: the adjugate matrix.
    """
    return _transpose(_cofactor_matrix(matrix))


def invert_matrix(matrix: list) -> list:
    """
    Inverts a square matrix using the determinant and adjugate.

    Args:
        matrix (list): a square matrix represented as a list of lists

    Returns:
        list: the inverse of the matrix
    """
    det = determinant(matrix)
    if det == 0:
        print("Matrix is singular. Inverse does not exist.")
        return None
    adj = _adjugate_matrix(matrix)
    n = len(matrix)
    inverse = [[adj[i][j] / det for j in range(n)] for i in range(n)]
    return inverse


def find_all_dividers_of_number(num: int) -> list:
    """
    Find all the dividers of the size of the given number.

    Args:
        num (int): number to find the dividers of.

    Returns:
        List[int]: sorted array with all the dividers of the given number

    Example:

        >>> find_all_dividers_of_number(12)
        [1, 2, 3, 4, 6, 12]

    """

    divisors = set()
    for i in range(1, int(np.sqrt(num)) + 1):
        if num % i == 0:
            divisors.add(i)
            divisors.add(num // i)

    return sorted(list(divisors))
