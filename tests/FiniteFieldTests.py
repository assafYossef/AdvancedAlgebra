import galois
import numpy as np

from galwa.fields import FiniteField


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
                    continue
                
                # divide the elements and then convert it to galois type
                sum_result = gf(str_rep((element1 / element2).a))

                assert sum_result == gf_element1 / gf_element2
    

    def test_inverse(self):
        p = 97

        poly = [5, 96, 1]

        


        
    
    # def test_order_of_elements(self):
    #     p = 7

    #     poly = [4, 0, 6, 1]

    #     gf = galois.GF(p ** (len(poly) - 1), repr="poly", irreducible_poly=poly[::-1])
    #     l = FiniteField(p=p, f=np.array(poly))

    #     assert len(l.elements) == len(gf.elements)

    #     for element in l.elements:
    #         gf_element = gf(str_rep(element.a))

    #         # skip on the zero element
    #         if gf_element == gf(0):
    #             continue

    #         try:
    #             assert gf_element.multiplicative_order() == element.order
    #         except AssertionError as e:
    #             print(e)
    #             print(f"Galois Result: {gf_element.multiplicative_order()}") # 6
    #             print(f"Our Result: {element.order}")  # 42
    #             print(f"Polynom Representaion: {str_rep(element.a)}") # 3
    #             print(f"Our Representions: {str(element)}")
    #             print(f"Vector Representaion: {element.a}")
    #             print()


    # def test(self):
    #     p = 7

    #     poly = [4, 0, 6, 1]

    #     gf = galois.GF(p**(len(poly) - 1), repr="poly", irreducible_poly=poly[::-1])
    #     l = FiniteField(p=p, f=np.array(poly), representation="vector")

    #     for element in l.elements:
    #         if str_rep(element.a) == str_rep([3, 0, 0]):
    #             print("here")
    #             gf_element = gf(str_rep(element.a))
    #             print(element.order)
                
                # counter = 1
                # temp_element = element
                # while gf_element != gf(1):
                #     temp_element = element * temp_element
                #     gf_element = gf(str_rep(temp_element.a))
                #     counter += 1
                
                # print(counter)
                # print(element.multiplicative_order())



    # def test_gln(self):
    #     p = 2

    #     poly = [1,1,1]

    #     gf = galois.GF(p**2, repr="poly", irreducible_poly=poly)
        
    #     l = FiniteField(p=2, f=np.array(poly), representation="vector")

    #     assert len(l.elements) == len(gf.elements)
    #     # convert the elements to galois elements
    #     for element in l.elements:
    #         convert.append(gf(str_rep(element.a)))
        
        # for element in gf.elements:
        #     if element.characteristic_poly() in convert:
        #         print("yes")
        
        # convert = np.array(convert)
        # assert (convert == gf.elements).all()
