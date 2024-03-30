from itertools import product
import numpy as np

from PrimeFieldElementClass import PrimeFieldElement
from Polynom import Polynom

class FiniteField:
    def __init__(self, p: int, f: np.ndarray):
        self.p = p

        # convert the array to include PrimeFieldElements
        polynom = np.ndarray((len(f),), dtype=PrimeFieldElement)
        for index in range(len(f)):
            polynom[index] = PrimeFieldElement(a=f[index], p=p)

        # validate polynom over Fp
        if not self._check_that_f_is_valid_over_p(polynom):
            raise Exception("polynom is not valid over Fp")
        elif not self._check_that_f_is_irreducible(polynom):
            raise Exception("polynom is reducible over Fp")
        else:
            # polynom is valid over Fp
            self.f = Polynom(poly=polynom)
        
        self._field_elements = []
        self._field_primitives = []
        self._field_generators = []
        self._field_general_element = []

        for pow_i in range(len(self.f) - 1):
            self._field_general_element.append(f"a{pow_i}")
        
        self._field_general_element = Polynom(poly=self._field_general_element)

    def _check_that_f_is_irreducible(self, polynom: np.ndarray) -> bool:
        # if the polynom degree is greater than 3, we need to belive them it's true
        if self.p > 3:
            return True

        # else p is 2 or 3
        Fp = np.ndarray((self.p - 1,), dtype=PrimeFieldElement)

        for i in range(1, self.p):
            Fp[i - 1] = PrimeFieldElement(i, self.p)

        irreducible = True
        for element in Fp:
            result = PrimeFieldElement(0, self.p)
            for pow_i in range(len(polynom)):
                result += polynom[pow_i] * (element ** pow_i)

            if result.a == 0:
                irreducible = False
                break

        return irreducible

    def _check_that_f_is_valid_over_p(self, polynom: np.ndarray):
        for element in polynom:
            if abs(element.a) >= self.p:
                return False

        return True
    
    def __format__(self, __format_spec: str) -> str:
        return f"F{self.p}[X]/<{self.irreducible_poly}>"

    @property
    def irreducible_poly(self):
        return self.f.str_rep

    # find generators over extended field
    # TODO run on all extended field elements
    # TODO check if the order of each element is equal to (p^r) - 1 | where r is the degree of f(x)
    @property
    def generators(self):
        pass

    # represent the genral element that belongs to this field
    @property
    def general_element(self):
       return self._field_general_element.str_rep
    
    # find all the elements in the extended field
    @property
    def elements(self):
        if len(self._field_elements) == 0:            
            n = self.p ** len(self.f)
            items = list(product(range(0, self.p), repeat=len(self._field_general_element)))
            self._field_elements = np.ndarray((len(items,),), dtype=Polynom)

            for index in range(len(items)):
                self._field_elements[index] = Polynom(items[index])
    
        return self._field_elements


    # find all the primitives in the extended field
    # TODO - How do we find primitives? - their order is the size of the field minus 1
    @property
    def primitives(self):
        pass

    # find GLn matrice
    # TODO find the

# p = 2
# extended_field = FiniteField(p=p, f=np.array([1, 1, 0, 0, 1]))
# print(extended_field.irreducible_poly)
# print(extended_field.general_element)
# for element in extended_field.elements:
#     print(element)

# print(len(extended_field.elements))