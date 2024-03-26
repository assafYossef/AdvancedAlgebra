import numpy as np

class FiniteField:
    def __init__(self, p, f: np.ndarray):
        self.p = p
        self.f = f

    # TODO check f is irreducible over p. (we can only check that if the degree of f is 2 or 3 else we need to belive them)
    def _check_that_f_is_irreducible(self):
        pass

    # TODO check that f coefficients are valid over p
    def _check_that_f_is_valid_over_p(self):
        pass

    # find generators over extended field
    # TODO run on all extended field elements
    # TODO check if the order of each element is equal to (p^r) - 1 | where r is the degree of f(x)
    def _find_generator(self):
        pass

    # find all the elements in the extended field
    # TODO find a solution for f(x) = 0 mod(f(x))
    # TODO calculate its pow from 0 to r-1 modulu p | where r is the degree of f(x)
    def _find_elements(self):
        pass

    # find all the primitives in the extended field
    # TODO - How do we find primitives? 
    def _find_primitives(self):
        pass