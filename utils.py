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
    def wrapper(self, other):
        if self.p != other.p:
            raise ValueError("Operands must belong to the same prime field")
        return func(self, other)
    return wrapper