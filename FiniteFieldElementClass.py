import numpy as np

from PrimeFieldElementClass import PrimeFieldElement

from FiniteFieldClass import FiniteField

#TODO build GLn matrice - > open EmbededToGln.txt to find the steps to do it | Should we do it here or in FiniteFieldClass?

class FiniteFieldElement:
    def __init__(self, a: np.ndarray, p: int):
        # convert the array to include PrimeFieldElements
        polynom = np.ndarray((len(a),), dtype=PrimeFieldElement)
        for index in range(len(a)):
            polynom[index] = PrimeFieldElement(a=a[index], p=p)
        
        self.a = polynom
        self.p = p
        self.variable_order = -1

    # TODO implement pretty printing for polynomials
    def __repr__(self):
        pass

    # TODO implement pretty printing for polynomials (should ovveride this function or the __repr__?)
    def __str__(self):
        poly = []
        for pow_i in range(len(self.a) - 1, -1, -1):
            if pow_i == 0:
                poly.append(str(self.a[pow_i].a))
            elif self.a[pow_i].a == 0:
                continue
            elif self.a[pow_i].a == 1:
                if pow_i == 1:
                    poly.append("X")
                else:
                    poly.append(f"X^{pow_i}")
            else:
                if pow_i == 1:
                    poly.append(f"{self.a[pow_i].a}*X")
                else:
                    poly.append(f"{self.a[pow_i].a}*X^{pow_i}")

        return " + ".join(poly)

    # TODO sum coefficients and at the end do modulu p
    def __add__(self, other: "FiniteFieldElement") -> "FiniteFieldElement":
        return self.__class__(self.a + other, self.p)

    # TODO sub coefficients and at the end do modulu p
    def __sub__(self, other: "FiniteFieldElement") -> "FiniteFieldElement":
        return self.__class__(self.a - other, self.p)

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

    # check that a is a valid element in the extenede field
    def _check_that_a_is_in_field(self, field: FiniteField):
        belong_to_l = True
        for element in self.a:
            if element.p != field.p:
                belong_to_l = False
                break
        
        return belong_to_l and len(self.a) < len(field.f)
        
