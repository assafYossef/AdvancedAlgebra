from utils import xgcd, same_prime_field


class PrimeFieldElement:
    def __init__(self, a, p):
        self.a = a % p
        self.p = p
        self.i = self._find_inverse()

    def __str__(self):
        return f"{self.a} mod {self.p}"

    def __repr__(self):
        return f"PrimeFieldElement(value= {self.a},prime= {self.p})"

    @same_prime_field
    def __eq__(self, other):
        return self.a == other.a and self.p == other.p and self.inverse == other.inverse

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
        if other.a == 0:
            raise ValueError("Division by zero")
        return self.__class__((self.a * other.i) % self.p, self.p)

    def __pow__(self, n):
        """
        Exponentiation of a prime field element
        Args:
            n (int): the exponent
        Returns:
            PrimeFieldElement: the result of the exponentiation
        """
        if n == 0:
            return self.__class__(1, self.p)
        if n < 0:
            return self.__class__(pow(self.i, -n, self.p), self.p)
        return self.__class__(pow(self.a, n, self.p), self.p)

    @property
    def inverse(self):
        """
        Returns the unit of the element
        Returns:
            int: the unit of the element
        """
        return PrimeFieldElement(self.i, self.p)

    def _find_inverse(self):
        """
        find the unit of the element
        The unit of an element a is the element b such that a*b = 1 mod p, to find b we use the extended euclidean algorithm
        Returns:
            int: the unit of the element

        Raises:
            ValueError: if the element does not have an inverse
        """
        if self.a == 0:
            return 0
        d, s, t = xgcd(self.a, self.p)
        if d != 1:
            raise ValueError("Element does not have an inverse")
        return s % self.p
    


    
