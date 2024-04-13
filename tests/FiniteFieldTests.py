import galois
import numpy as np
import random

from galwa.fields import FiniteField
from galwa.utils import bsgs


def str_rep(poly):
    str_poly = []
    for pow_i in range(len(poly) - 1, -1, -1):
        if pow_i == 0:
            str_poly.append(str(poly[pow_i]))
        elif poly[pow_i] == 0:
            continue
        elif poly[pow_i] == 1:
            if pow_i == 1:
                str_poly.append("X")
            else:
                str_poly.append(f"X^{pow_i}")
        else:
            if pow_i == 1:
                str_poly.append(f"{poly[pow_i]}*X")
            else:
                str_poly.append(f"{poly[pow_i]}*X^{pow_i}")
    
    # remove unnecessary 0
    if str_poly[len(str_poly) - 1] == "0" and len(str_poly) > 1:
        str_poly = str_poly[:-1]

    return " + ".join(str_poly)


class FiniteFieldTests:

    def test_elements_in_field_easy(self):
        p = 2

        poly = [1,1,1]

        gf = galois.GF(p**(len(poly) - 1), repr="poly", irreducible_poly=poly[::-1])
        l = FiniteField(p=p, f=np.array(poly), representation="vector")

        assert len(l.elements) == len(gf.elements)

        convert = []
        # convert the elements to galois elements
        for element in l.elements:
            convert.append(gf(str_rep(element.a)))
        
        for element in gf.elements:
            assert element in convert


    
    def test_elements_in_field_medium(self):
        p = 2

        poly = [1, 1, 0, 1]

        gf = galois.GF(p**(len(poly) - 1), repr="poly", irreducible_poly=poly[::-1])
        l = FiniteField(p=p, f=np.array(poly), representation="vector")

        assert len(l.elements) == len(gf.elements)

        convert = []
        # convert the elements to galois elements
        for element in l.elements:
            convert.append(gf(str_rep(element.a)))
        
        for element in gf.elements:
            assert element in convert

    def test_elements_in_field_hard(self):
        p = 3

        poly = [2, 0, 0, 2, 1]

        gf = galois.GF(p**(len(poly) - 1), repr="poly", irreducible_poly=poly[::-1])
        l = FiniteField(p=p, f=np.array(poly), representation="vector")

        assert len(l.elements) == len(gf.elements)
        
        convert = []
        # convert the elements to galois elements
        for element in l.elements:
            convert.append(gf(str_rep(element.a)))
        
        for element in gf.elements:
            assert element in convert

    def test_elements_in_field_hell(self):
        p = 47

        poly = [5, 45, 1]

        gf = galois.GF(p**(len(poly) - 1), repr="poly", irreducible_poly=poly[::-1])
        l = FiniteField(p=p, f=np.array(poly), representation="vector")

        assert len(l.elements) == len(gf.elements)

        convert = []
        # convert the elements to galois elements
        for element in l.elements:
            convert.append(gf(str_rep(element.a)))
        
        for element in gf.elements:
            assert element in convert
    
    def test_sum_elements_in_field(self):
        p = 5

        poly = [2, 4, 1]

        gf = galois.GF(p**(len(poly) - 1), repr="poly", irreducible_poly=poly[::-1])
        l = FiniteField(p=p, f=np.array(poly), representation="vector")
        
        for element1 in l.elements:
            for element2 in l.elements:
                # convert the elements to be galois type
                gf_element1 = gf(str_rep(element1.a))
                gf_element2 = gf(str_rep(element2.a))
                
                # sum the elements and then convert it to galois type
                sum_result = gf(str_rep((element1 + element2).a))

                assert sum_result == gf_element1 + gf_element2
    
    def test_sub_elements_in_field(self):
        p = 5

        poly = [2, 4, 1]

        gf = galois.GF(p**(len(poly) - 1), repr="poly", irreducible_poly=poly[::-1])
        l = FiniteField(p=p, f=np.array(poly), representation="vector")
        
        for element1 in l.elements:
            for element2 in l.elements:
                # convert the elements to be galois type
                gf_element1 = gf(str_rep(element1.a))
                gf_element2 = gf(str_rep(element2.a))
                
                # sub the elements and then convert it to galois type
                sum_result = gf(str_rep((element1 - element2).a))

                assert sum_result == gf_element1 - gf_element2
    
    def test_mul_elements_in_field(self):
        p = 5

        poly = [2, 4, 1]

        gf = galois.GF(p**(len(poly) - 1), repr="poly", irreducible_poly=poly[::-1])
        l = FiniteField(p=p, f=np.array(poly), representation="vector")
        
        for element1 in l.elements:
            for element2 in l.elements:
                # convert the elements to be galois type
                gf_element1 = gf(str_rep(element1.a))
                gf_element2 = gf(str_rep(element2.a))

                # check if one of the elements is 0, so we skip it
                if gf_element1 == gf(0):
                    break

                if gf_element2 == gf(0):
                    continue
                
                # multiply the elements and then convert it to galois type
                sum_result = gf(str_rep((element1 * element2).a))

                assert sum_result == gf_element1 * gf_element2
    

    def test_divide_elements_in_field(self):
        p = 5

        poly = [2, 4, 1]

        gf = galois.GF(p**(len(poly) - 1), repr="poly", irreducible_poly=poly[::-1])
        l = FiniteField(p=p, f=np.array(poly), representation="vector")
        
        for element1 in l.elements:
            for element2 in l.elements:
                # convert the elements to be galois type
                gf_element1 = gf(str_rep(element1.a))
                gf_element2 = gf(str_rep(element2.a))

                # check if one of the elements is 0, so we skip it
                if gf_element1 == gf(0):
                    break

                if gf_element2 == gf(0):
                    try:
                        x = element1 / element2
                    except Exception as e:
                        assert type(e) == ValueError
                        assert e.args[0] == "Zero is not part of the multiplicative group, so it cannot be used in this operation"
                        continue

                gf_element1.multiplicative_order()
                
                # divide the elements and then convert it to galois type
                sum_result = gf(str_rep((element1 / element2).a))

                assert sum_result == gf_element1 / gf_element2


    def test_inverse(self):
        p = 97

        poly = [5, 96, 1]

        gf = galois.GF(p ** (len(poly) - 1), repr="poly", irreducible_poly=poly[::-1])
        l = FiniteField(p=p, f=np.array(poly))

        assert len(l.elements) == len(gf.elements)

        for element in l.elements:
            gf_element = gf(str_rep(element.a))

            # skip on the zero and element
            if gf_element == gf(0):
                try:
                    x = element ** -1
                except Exception as e:
                    assert type(e) == ValueError
                    assert e.args[0] == "Zero element doesn't have inverse"
                    continue

            assert gf_element ** -1 == gf(str_rep((element ** -1).a))
    

    def test_pow(self):
        p = 47

        poly = [5, 45, 1]

        gf = galois.GF(p ** (len(poly) - 1), repr="poly", irreducible_poly=poly[::-1])
        l = FiniteField(p=p, f=np.array(poly))

        for power in range(-5, 5):
            element = l.elements[random.randint(0, len(l.elements) - 1)]
            gf_element = gf(str_rep(element.a))

            # if we chose 0, we pick another one
            while gf_element == gf(0):
                element = l.elements[random.randint(0, len(l.elements) - 1)]
                gf_element = gf(str_rep(element.a))

            assert gf_element ** power == gf(str_rep((element ** power).a))


    def test_order_of_elements(self):
        p = 7

        poly = [4, 0, 6, 1]

        gf = galois.GF(p ** (len(poly) - 1), repr="poly", irreducible_poly=poly[::-1])
        l = FiniteField(p=p, f=np.array(poly))

        assert len(l.elements) == len(gf.elements)        

        for element in l.elements:
            gf_element = gf(str_rep(element.a))            

            # skip on the zero element
            if gf_element == gf(0):
                continue
            
            assert gf_element.multiplicative_order() == element.order

    
    def test_generators(self):
        p = 7

        poly = [3, 6, 1]

        gf = galois.GF(p ** (len(poly) - 1), repr="poly", irreducible_poly=poly[::-1])
        l = FiniteField(p=p, f=np.array(poly))

        generators = l.generators

        assert len(generators) == len(gf.primitive_elements)

        for generator in generators:
            gf_element = gf(str_rep(generator.a))

            assert gf_element in gf.primitive_elements
    
    def test_bsgs(self):
        p = 7

        poly = [3, 6, 1]

        gf = galois.GF(p ** (len(poly) - 1), repr="poly", irreducible_poly=poly[::-1])
        l = FiniteField(p=p, f=np.array(poly))

        generators = l.generators

        for element in l.elements:
            gf_element = gf(str_rep(element.a))

            if gf_element == gf(0):
                continue

            g = generators[random.randint(0, len(generators) - 1)]

            x = bsgs(generator=g, element=element, group_order=g.order)

            assert g ** x == element
        
        g = l.elements[random.randint(0, len(l.elements) - 1)]
        g_gf = gf(str_rep(g.a))

        while g_gf == gf(0):
            g = l.elements[random.randint(0, len(l.elements) - 1)]
            g_gf = gf(str_rep(g.a))
            
        for element in l.elements:
            gf_element = gf(str_rep(element.a))

            if gf_element == gf(0):
                continue

            try:
                x = bsgs(generator=g, element=element, group_order=g.order)
                assert g ** x == element
            except ValueError:
                x = gf_element.log(g_gf)
                assert g_gf ** x != gf_element

    def test_operations_on_elemetns_different_fields(self):
        p1 = 7
        poly1 = [3, 6, 1]

        l1 = FiniteField(p=p1, f=np.array(poly1))
        
        p2 = 5
        poly2 = [2, 4, 1]

        l2 = FiniteField(p=p2, f=np.array(poly2))

        for element1 in l1.elements:
            for element2 in l2.elements:
                
                # dont do check on zero element
                if np.all(element2.a == 0):
                    continue

                # try sum
                try:
                    x = element1 + element2
                except Exception as e:
                    assert type(e) == ValueError
                    assert e.args[0] == "Operands must belong to the same field"
                
                # try sub
                try:
                    x = element1 - element2
                except Exception as e:
                    assert type(e) == ValueError
                    assert e.args[0] == "Operands must belong to the same field"
                
                # try mul
                try:
                    x = element1 * element2
                except Exception as e:
                    assert type(e) == ValueError
                    assert e.args[0] == "Operands must belong to the same field"
                
                # try truediv
                try:
                    x = element1 / element2
                except Exception as e:
                    assert type(e) == ValueError
                    assert e.args[0] == "Operands must belong to the same field"
                
                break


