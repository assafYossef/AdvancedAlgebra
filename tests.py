from PrimeFieldElementClass import PrimeFieldElement

def test_prime_field_element():
    p = 11

    x = PrimeFieldElement(3, p)
    y = PrimeFieldElement(5, p)
    z = PrimeFieldElement(7, p)
    w = PrimeFieldElement(9, p)
    n = PrimeFieldElement(10, p)

    assert x + y == PrimeFieldElement(8, p)
    assert x - y == PrimeFieldElement(9, p)
    assert x * y == PrimeFieldElement(4, p)
    assert x / y == PrimeFieldElement(5, p)
    assert x.inverse == 4
    assert y.inverse == 9
    assert z.inverse == 8
    assert w.inverse == 5
    assert n.inverse == 10

    print("All tests passed")
test_prime_field_element()