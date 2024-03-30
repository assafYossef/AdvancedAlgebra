import numpy as np

from PrimeFieldElementClass import PrimeFieldElement

from FiniteFieldClass import FiniteField
from Polynom import Polynom

#TODO build GLn matrice - > open EmbededToGln.txt to find the steps to do it | Should we do it here or in FiniteFieldClass?

class FiniteFieldElement:
    def __init__(self, a: np.ndarray, l: FiniteField):
        # convert the array to include PrimeFieldElements
        polynom = np.ndarray((len(a),), dtype=PrimeFieldElement)
        for index in range(len(a)):
            polynom[index] = PrimeFieldElement(a=a[index], p=l.p)
        
        self.a = Polynom(poly=polynom)
        self.l = l
        
        # check a is in l
        if not self._check_that_a_is_in_field():
            raise Exception(f"{self.a} does not belong to {l}")
        self.variable_order = -1

    # TODO implement pretty printing for polynomials
    def __repr__(self):
        pass

    # TODO implement pretty printing for polynomials (should ovveride this function or the __repr__?)
    def __str__(self):
        return self.a.str_rep

    # TODO sum coefficients and at the end do modulu p
    def __add__(self, other: "FiniteFieldElement") -> "FiniteFieldElement":
        return self.__class__(self.a + other.a, self.p)

    # TODO sub coefficients and at the end do modulu p
    def __sub__(self, other: "FiniteFieldElement") -> "FiniteFieldElement":
        return self.__class__(self.a - other.a, self.p)

    # TODO convert to GLn matrice
    # TODO multiply both matrices and convert back to polynomial
    def __mul__(self, other: "FiniteFieldElement"):
        pass

    # TODO convert to GLn matrice
    # TODO find the opposite matrice of other
    # TODO multiply them and convert the result back to polynomial
    def __truediv__(self, other: "FiniteFieldElement"):        
        pass

    # TODO find the order of the current element in the extended field
    @property
    def order(self):
        pass

    # check that a is a valid element in the exteneded field
    def _check_that_a_is_in_field(self):
        belong_to_l = True
        for element in self.a:
            if element.p != self.l.p:
                belong_to_l = False
                break
        
        return belong_to_l and len(self.a) < len(self.l.f)
        


# extended_field = FiniteField(p=2, f=np.array([1, 1, 0, 0, 1]))
# a = FiniteFieldElement(a=np.array([1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]), l=extended_field)