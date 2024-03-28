from utils import xgcd, same_prime_field


class PrimeFieldElement:
    def __init__(self, a, p):
        self.a = a % p
        self.p = p
        self.i = self._find_inverse()

    # overload print function
    def __str__(self):
        return f"{self.a} mod {self.p}"

    def __eq__(self, other):
        return self.a == other.a and self.p == other.p and self.inverse == other.inverse
    # overload add function
    @same_prime_field
    def __add__(self, other):
        """
        Addition of two prime field elements
        Args:
            other (PrimeFieldElement): the element to add
        Returns:
            PrimeFieldElement: the result of the addition
        """
        return self.__class__((self.a + other.a) % self.p, self.p)

    @same_prime_field
    def __sub__(self, other):
        """
        Subtraction of two prime field elements
        Args:
            other (PrimeFieldElement): the element to subtract
        Returns:
            PrimeFieldElement: the result of the subtraction
        """
        return self.__class__((self.a - other.a) % self.p, self.p)

    @same_prime_field
    def __mul__(self, other):
        """
        Multiplication of two prime field elements
        Args:
            other (PrimeFieldElement): the element to multiply by
        Returns:
            PrimeFieldElement: the result of the multiplication
        """
        return self.__class__((self.a * other.a) % self.p, self.p)

    @same_prime_field
    def __truediv__(self, other):
        """
        Division of two prime field elements
        Args:
            other (PrimeFieldElement): the element to divide by
        Returns:
            PrimeFieldElement: the result of the division
        """
        return self.__class__((self.a * other.inverse) % self.p, self.p)

    @property
    def inverse(self):
        """
        Returns the unit of the element
        Returns:
            int: the unit of the element
        """
        return self.i

    def _find_inverse(self):
        """
        find the unit of the element
        The unit of an element a is the element b such that a*b = 1 mod p, to find b we use the extended euclidean algorithm
        Returns:
            int: the unit of the element
        """
        d, s, t = xgcd(self.a, self.p)
        if d != 1:
            raise ValueError("Element does not have an unit")
        return s % self.p
    


    
