import numpy as np

class Polynom:
    def __init__(self, poly: np.ndarray):
        self.poly = poly
    
    def __add__(self, other):
        return self.poly + other.poly
    
    def __sub__(self, other):
        return self.poly - other.poly
    
    def __str__(self) -> str:
        return self.str_rep
    
    def __format__(self, __format_spec: str) -> str:
        return self.str_rep

    @property
    def str_rep(self):
        str_poly = []
        for pow_i in range(len(self.poly) - 1, -1, -1):
            if pow_i == 0:
                str_poly.append(str(self.poly[pow_i]))
            elif self.poly[pow_i] == 0:
                continue
            elif self.poly[pow_i] == 1:
                if pow_i == 1:
                    str_poly.append("X")
                else:
                    str_poly.append(f"X^{pow_i}")
            else:
                if pow_i == 1:
                    str_poly.append(f"{self.poly[pow_i]}*X")
                else:
                    str_poly.append(f"{self.poly[pow_i]}*X^{pow_i}")
        
        # remove unnecessary 0
        if str_poly[len(str_poly) - 1] == "0" and len(str_poly) > 1:
            str_poly = str_poly[:-1]

        return " + ".join(str_poly)

    def __len__(self):
        return len(self.poly)

    def __iter__(self):
        return self.poly.__iter__()
    