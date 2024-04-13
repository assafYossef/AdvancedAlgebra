import galois
from galwa.elements import PrimeFieldElement


class PrimeFieldTests:
    def __init__(self, p):
        self.p = p

    def test_sum_members_in_field(self):
        """
        Test the sum of members in the field

        Raises:
            AssertionError: If the sum of the members is not the expected one.
        """
        p = self.p

        gf = galois.GF(p)

        zero = PrimeFieldElement(0, 11)
        x = PrimeFieldElement(3, p)
        y = PrimeFieldElement(5, p)
        z = PrimeFieldElement(7, p)
        w = PrimeFieldElement(9, p)
        n = PrimeFieldElement(10, p)

        assert x + y == PrimeFieldElement(int(gf(x.a) + gf(y.a)), p)
        assert x + y + z == PrimeFieldElement(int(gf(x.a) + gf(y.a) + gf(z.a)), p)
        assert w + n == PrimeFieldElement(int(gf(w.a) + gf(n.a)), p)
        assert x + y + z + w + n == PrimeFieldElement(int(gf(x.a) + gf(y.a) + gf(z.a) + gf(w.a) + gf(n.a)), p)
        assert zero + x == x

    def test_sub_members_in_field(self):
        """
        Test the subtraction of members in the field

        Raises:
            AssertionError: If the subtraction of the members is not the expected one.
        """
        p = self.p

        gf = galois.GF(p)

        zero = PrimeFieldElement(0, 11)
        x = PrimeFieldElement(3, p)
        y = PrimeFieldElement(5, p)
        z = PrimeFieldElement(7, p)
        w = PrimeFieldElement(9, p)
        n = PrimeFieldElement(10, p)

        assert x - y == PrimeFieldElement(int(gf(x.a) - gf(y.a)), p)
        assert y - x == PrimeFieldElement(int(gf(y.a) - gf(x.a)), p)
        assert z - w - n == PrimeFieldElement(int(gf(z.a) - gf(w.a) - gf(n.a)), p)
        assert zero + x == x

    @staticmethod
    def test_sum_members_different_fields():
        """
        Test the sum of members in different fields, which should raise a ValueError

        Raises:
            AssertionError: If the sum of the members doesnt throw a ValueError
        """
        p = 11
        p2 = 13

        a = PrimeFieldElement(3, p)
        b = PrimeFieldElement(3, p2)

        try:
            c = a + b
        except Exception as e:
            assert type(e) == ValueError
            assert e.args[0] == "Operands must belong to the same prime field"

    @staticmethod
    def test_sub_members_different_fields():
        """
        Test the subtraction of members in different fields, which should raise a ValueError

        Raises:
            AssertionError: If the subtraction of the members doesnt throw a ValueError
        """
        p = 11
        p2 = 13

        a = PrimeFieldElement(3, p)
        b = PrimeFieldElement(3, p2)

        try:
            c = a - b
        except Exception as e:
            assert type(e) == ValueError
            assert e.args[0] == "Operands must belong to the same prime field"

    def test_mul_members_in_field(self):
        """
        Test the multiplication of members in the field

        Raises:
            AssertionError: If the multiplication of the members is not the expected one.
        """
        p = self.p

        gf = galois.GF(p)

        zero = PrimeFieldElement(0, 11)
        x = PrimeFieldElement(3, p)
        y = PrimeFieldElement(5, p)
        z = PrimeFieldElement(7, p)
        w = PrimeFieldElement(9, p)
        n = PrimeFieldElement(10, p)

        assert x * y == PrimeFieldElement(int(gf(x.a) * gf(y.a)), p)
        assert y * x == PrimeFieldElement(int(gf(y.a) * gf(x.a)), p)
        assert z * w * n == PrimeFieldElement(int(gf(z.a) * gf(w.a) * gf(n.a)), p)
        assert zero * x == zero

    @staticmethod
    def test_mul_members_different_fields():
        """
        Test the multiplication of members in different fields, which should raise a ValueError

        Raises:
            AssertionError: If the multiplication of the members doesnt throw a ValueError
        """
        p = 11
        p2 = 13

        a = PrimeFieldElement(3, p)
        b = PrimeFieldElement(3, p2)

        try:
            c = a * b
        except Exception as e:
            assert type(e) == ValueError
            assert e.args[0] == "Operands must belong to the same prime field"

    def test_divide_members_in_field(self):
        """
        Test the division of members in the field

        Raises:
            AssertionError: If the division of the members is not the expected one.
        """
        p = self.p

        gf = galois.GF(p)

        zero = PrimeFieldElement(0, 11)
        x = PrimeFieldElement(3, p)
        y = PrimeFieldElement(5, p)
        z = PrimeFieldElement(7, p)
        w = PrimeFieldElement(9, p)
        n = PrimeFieldElement(10, p)

        assert x / y == PrimeFieldElement(int(gf(x.a) / gf(y.a)), p)
        assert y / x == PrimeFieldElement(int(gf(y.a) / gf(x.a)), p)
        assert z / w / n == PrimeFieldElement(int(gf(z.a) / gf(w.a) / gf(n.a)), p)
        assert zero / x == zero

        try:
            c = x / zero
        except Exception as e:
            assert type(e) == ValueError
            assert e.args[0] == "Division by zero"

    @staticmethod
    def test_divide_members_different_fields():
        """
        Test the division of members in different fields, which should raise a ValueError

        Raises:
            AssertionError: If the division of the members doesnt throw a ValueError
        """
        p = 11
        p2 = 13

        a = PrimeFieldElement(3, p)
        b = PrimeFieldElement(3, p2)

        try:
            c = a / b
        except Exception as e:
            assert type(e) == ValueError
            assert e.args[0] == "Operands must belong to the same prime field"

    def test_find_legal_inverse(self):
        """
        Test the inverse of members in the field

        Raises:
            AssertionError: If the inverse of the members is not the expected one.
        """
        p = self.p

        gf = galois.GF(p)

        x = PrimeFieldElement(3, p)
        y = PrimeFieldElement(5, p)
        z = PrimeFieldElement(7, p)
        w = PrimeFieldElement(9, p)
        n = PrimeFieldElement(10, p)

        assert x.inverse == int(gf(x.a) ** -1)
        assert y.inverse == int(gf(y.a) ** -1)
        assert z.inverse == int(gf(z.a) ** -1)
        assert w.inverse == int(gf(w.a) ** -1)
        assert n.inverse == int(gf(n.a) ** -1)

    def test_illegal_inverse(self):
        """
        Test the inverse of zero, which should be None

        Raises:
            AssertionError: If the inverse of zero is not None
        """
        p = self.p

        zero = PrimeFieldElement(0, p)

        assert zero.inverse is None
