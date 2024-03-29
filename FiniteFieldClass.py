import numpy as np

from PrimeFieldElementClass import PrimeFieldElement

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
            self.f = polynom
        
        self.field_elements = []
        self.field_primitives = []
        self.field_generators = []

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

    @property
    def irreducible_poly(self):
        poly = []
        for pow_i in range(len(self.f) - 1, -1, -1):
            if pow_i == 0:
                poly.append(str(self.f[pow_i].a))
            elif self.f[pow_i].a == 0:
                continue
            elif self.f[pow_i].a == 1:
                if pow_i == 1:
                    poly.append("X")
                else:
                    poly.append(f"X^{pow_i}")
            else:
                if pow_i == 1:
                    poly.append(f"{self.f[pow_i].a}*X")
                else:
                    poly.append(f"{self.f[pow_i].a}*X^{pow_i}")

        return " + ".join(poly)

    # find generators over extended field
    # TODO run on all extended field elements
    # TODO check if the order of each element is equal to (p^r) - 1 | where r is the degree of f(x)
    @property
    def generators(self):
        pass

    # find all the elements in the extended field
    # TODO find a solution for f(x) = 0 mod(f(x))
    # TODO calculate its pow from 0 to r-1 modulu p | where r is the degree of f(x)
    @property
    def elements(self):
        pass

    # find all the primitives in the extended field
    # TODO - How do we find primitives? 
    @property
    def primitives(self):
        pass

    # find GLn matrice
    # TODO find the

p = 47
extended_field = FiniteField(p=p, f=np.array([5, 40, 8, 0, 1]))
print(extended_field.irreducible_poly)