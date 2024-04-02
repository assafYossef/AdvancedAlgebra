import numpy as np


def xgcd(a,b):
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
    Helper function to check if the operands belong to the same prime field
    This is used for PrimeFieldElement class, which we implemented but didnt use it in FiniteFieldElement or FiniteField.
    """
    def wrapper(self, other):
        if self.p != other.p:
            raise ValueError("Operands must belong to the same prime field")
        return func(self, other)
    return wrapper


def zero_element_check(func):
    """
    Helper function to check if an element is zero for the multiplicative group, so multiplication or division cannot be done
    """
    def wrapper(self, other):
        if np.all(self.a == 0) or np.all(other.a == 0):
            raise ValueError("Zero is not part of the multiplicative group, so it cannot be used in this operation")
        return func(self, other)
    return wrapper


def valid_repr(func):
    """
    Helper function to check if the representation is valid
    A valid representation is either "polynomial", "vector" or "matrix"
    """
    def wrapper(self, value):
        assert value in ["polynomial", "vector", "matrix"], "Invalid representation"
        return func(self, value)
    return wrapper


def refactor_polynom_terms(poly_repr: str) -> str:
    """
    Helper function to refactor the polynomial representation
    The polynomial representation of numpy returns float coefficients, as well it doesnt drop zero coefficients
    This function will refactor the polynomial representation to a more human readable form

    Example:
    poly_repr = "1.0 x^2 + 0.0 x^1 + 3.0 x^0"
    refactor_polynom_terms(poly_repr) -> "x^2 + 3"

    Args:
        poly_repr (str): polynomial representation

    Returns:
        str: refactored polynomial representation
    """
    poly_repr = f"{poly_repr}".replace(".0", "")
    new_str = ""
    split_str = poly_repr.split(" ")
    for term in split_str:
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
                # remove the last term, it will be a arthimatic operator and space
                new_str = new_str[:-2]
    # remove space from end and beginning
    new_str = new_str.strip("+- ")
    # if new_str is empty return 0
    return new_str if new_str else "0"


def pad_element(element, f):
    """
    Helper function to pad the element with zeros to match the size of f
    Args:
        element (np.array): element to pad
        f (np.array): irreducible polynomial to match the size

    Returns:
        np.array: padded element (if needed)
    """
    return np.pad(element, pad_width=(0, len(f) - len(element) -1 ), mode='constant')