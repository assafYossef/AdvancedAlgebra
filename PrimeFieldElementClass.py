# TODO decide if to use this function or a library
def xgcd(a,b):
    if b == 0:
        if a < 0:
            return -a, -1, 0
        return a, 1, 0
    else:
        q, r = a // b, a % b
        d, s, t = xgcd(b, r)
        return d, t, s - q * t


# TODO decorator for checking that comparison between two different prime field elements have the same prime
# TODO add edge cases
class PrimeFieldElement:
    def __init__(self, a, p):
        self.a = a % p
        self.p = p
    
    # overload print function
    def __str__(self):
        return str(self.a)
    
    # overload add function
    def __add__(self, other):
        return self.__class__((self.a + other.a) % self.p, self.p)
    
    # overload subtraction function
    def __sub__(self, other):
        return self.__class__((self.a - other.a) % self.p, self.p)
    
    def __mul__(self, other):
        return self.__class__((self.a * other.a) % self.p, self.p)

    # TODO call the extended euclidean algorithm
    # TODO make sure the element has inverse (means its not equal to 0 modulu p)
    @property
    def inverse(self):
        pass
    
    # TODO validate the other has inverse so we wont get an error
    def __truediv__(self, other):
        return self.__class__((self.a * other.inverse) % self.p)
    

    



p = 2

x = PrimeFieldElement(8, p)
y = PrimeFieldElement(9, p)


z = x - y + x

print(z)