import math

# extended euclidean algorithm
def xgcd(a: int,b: int) -> int:
    if b == 0:
        if a < 0:
            return -a, -1, 0
        return a, 1, 0
    else:
        q, r = a // b, a % b
        d, s, t = xgcd(b, r)
        return d, t, s - q * t

# TODO add edge cases
class PrimeFieldElement:
    def __init__(self, a, p):
        if a >= p or a < 0:
            raise Exception(f"{a} is not in F of {p}")
        
        if a != 0:
            if math.gcd(p, a) != 1:
                raise Exception(f"{a} is not in F of {p}")
        
        self.p = p
        self.a = a
    
    # decorator for validating that operations will be on elements with the same prime p
    def validate_other_p(func):
        def wrapper(*args, **kw):
            if args[0].p != args[1].p:
                raise Exception("Not Equal P!")

            return func(*args, **kw)

        return wrapper
    
    def _exponentional_by_squering(self, base: "PrimeFieldElement", pow_number: int) -> "PrimeFieldElement":
        if pow_number == 0:
             return PrimeFieldElement(1, base.p)
        
        if pow_number % 2 == 0:
            return self._exponentional_by_squering(base * base, pow_number / 2)
        else:
            return base * self._exponentional_by_squering(base * base, math.floor(pow_number / 2))
    
    # overload print function
    def __str__(self):
        return str(self.a)
    
    # overload add function
    @validate_other_p
    def __add__(self, other: "PrimeFieldElement") -> "PrimeFieldElement":
        return self.__class__((self.a + other.a) % self.p, self.p)
    
    # overload subtraction function
    @validate_other_p
    def __sub__(self, other: "PrimeFieldElement") -> "PrimeFieldElement":
        return self.__class__((self.a - other.a) % self.p, self.p)
    
    @validate_other_p
    def __mul__(self, other: "PrimeFieldElement") -> "PrimeFieldElement":
        return self.__class__((self.a * other.a) % self.p, self.p)

    def __pow__(self, pow_number: int) -> "PrimeFieldElement":        
        is_negative = pow_number < 0

        pow_number = abs(pow_number)
        result = self._exponentional_by_squering(self, pow_number)

        if is_negative:
            result = result.inverse

        return result

    @property
    def inverse(self) -> "PrimeFieldElement":
        if self.a == 0:
            raise Exception("{a} doesnt have inverse in Fp")
        
        gcd, s, _ = xgcd(self.a, self.p)

        if gcd != 1:
            raise Exception("a is not foreign for p")
        
        if s < 0:
            s = self.p + s

        return self.__class__(s, self.p)
    
    @validate_other_p
    def __truediv__(self, other: "PrimeFieldElement") -> "PrimeFieldElement":
        return self.__class__((self.a * other.inverse.a) % self.p, self.p)