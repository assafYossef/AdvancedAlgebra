import numpy as np

from FiniteFieldClass import FiniteField

#TODO build GLn matrice - > open EmbededToGln.txt to find the steps to do it | Should we do it here or in FiniteFieldClass?

class FiniteFieldElement:
    def __init__(self, a: np.ndarray):
        self.a = a

    # TODO implement pretty printing for polynomials
    def __repr__(self):
        pass

    # TODO implement pretty printing for polynomials (should ovveride this function or the __repr__?)
    def __str__(self):
        pass

    # TODO sum coefficients and at the end do modulu p
    def __add__(self, other):
        pass

    # TODO sub coefficients and at the end do modulu p
    def __sub__(self, other):
        pass

    # TODO convert to GLn matrice
    # TODO multiply both matrices and convert back to polynomial
    def __mul__(self, other):
        pass

    # TODO convert to GLn matrice
    # TODO find the opposite matrice of other
    # TODO multiply them and convert the result back to polynomial
    def __truediv__(self, other):        
        pass

    # TODO find the order of the current element in the extended field
    @property
    def order(self):
        pass

    # check that a is a valid element in the extenede field
    # TODO check each element in a if in Fp
    # TODO check if a length is less the degree of f(x)
    def _check_that_a_is_in_field(self, field:FiniteField):
        # return np.all(self.a < field.p) and np.all(self.a >= 0) and len(self.a) <= len(field.f)x     - 1
        pass