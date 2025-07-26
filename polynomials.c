#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "polynomials.h"


struct cyclotomic_cosets*
generate_cosets_from_n_needed_elements(struct finite_field* field, int n)
{
    int was = 0;
    struct cyclotomic_coset** cosets_from = malloc(sizeof(struct cyclotomic_coset*) * n);
    int num_of_elemets = 0;
    uint64_t p = pow(field->characteristic, field->power) - 1;
    for(int i = 1; i <= n; ++i){
        if((was & (1 << i)) != 0){
            continue;
        }
        cosets_from[num_of_elemets] = malloc(sizeof(struct cyclotomic_coset));
        cosets_from[num_of_elemets]->from_el = i;
        cosets_from[num_of_elemets]->size = 1;
        int k = i;
        was |= (1 << i);
        if(field->characteristic == 2){
            k <<= 1;
        }else{
            k *= field->characteristic;
        }
        while(k != i){
            was |= (1 << k);
            if(field->characteristic == 2){
                k <<= 1;
            }else{
                k *= field->characteristic;
            }
            
            k %= p;

            (cosets_from[num_of_elemets]->size)++;
        }
        ++num_of_elemets;
    }
    struct cyclotomic_cosets* to_return = malloc(sizeof(struct cyclotomic_cosets));
    to_return->cosets = cosets_from;
    to_return->number_of_cosets = num_of_elemets;
    return to_return;
}

void
free_cosets_struct(struct cyclotomic_cosets* to_free)
{
    for(int i = 0; i < to_free->number_of_cosets; ++i){
        free(to_free->cosets[i]);
    }
    free(to_free->cosets);
    free(to_free);
}

void
free_extended_polynomial(struct extended_polynomial* to_free)
{
    free(to_free->coefs);
    free(to_free);
}

void
print_extended_polynomial(uint64_t* to_print)
{
    //if(to_print == NULL){
    //    printf("Poly is NULL\n");
    //    return;
    //}
    //printf("Degree: %d\t", to_print->degree);
    //if(to_print->degree == 0) {
    //    printf("\n\n");
    //    return;
    //}
    //for(int i = 0; i < to_print->degree; ++i){
    //    if(to_print->coefs[i] != 0){
    //        printf("%lu * x^%d + ", to_print->coefs[i], i);   
    //    }
    //}
    for(int i = 0; i < MAX_DEGREE; ++i){
        printf("%lu \t", to_print[i]);   
    }
    printf("\n\n");
}

void
copy_extended_polynomial(uint64_t* poly, uint64_t* res)
{
    for(int i = 0; i < MAX_DEGREE; ++i){
        res[i] = poly[i];
    }
}

void
add_two_extended_polynomials(struct finite_field* field, uint64_t* poly1, uint64_t* poly2)
{
    for(int i = 0; i < MAX_DEGREE; ++i){
        poly1[i] = add_in_field(field, poly1[i], poly2[i]);
    }
}

void
multiply_two_extended_polynomials(struct finite_field* field, uint64_t* poly1, uint64_t* poly2)
{
    uint64_t tmp[MAX_DEGREE];
    memcpy(tmp, poly1, sizeof(uint64_t) * MAX_DEGREE);

    memset(poly1, 0, sizeof(uint64_t) * MAX_DEGREE);

    int deg_poly1 = get_degree(tmp), deg_poly2 = get_degree(poly2);

    for (int i = 0; i <= deg_poly1; ++i) {
        uint64_t a = tmp[i];
        if (a == 0) continue;
        for (int j = 0; j <= deg_poly2 && (i + j) < MAX_DEGREE; ++j) {
            uint64_t b = poly2[j];
            if (b == 0) continue;
            uint64_t prod = multiply_in_field(field, a, b);
            poly1[i + j] = add_in_field(field, poly1[i + j], prod);
        }
    }
    
}


void
construct_minimal_polynomial_from_coset(struct finite_field* field, struct cyclotomic_coset* coset, uint64_t res[MAX_DEGREE])
{
    int base_element = coset->from_el;

    uint64_t first_multiplier[MAX_DEGREE];
    uint64_t second_multiplier[MAX_DEGREE];
    memset(first_multiplier, 0, MAX_DEGREE * sizeof(uint64_t));
    uint64_t primitive_power = find_primitive_in_power(field, base_element);
    first_multiplier[1] = 1;
    
    first_multiplier[0] = construct_inverse_element_add(field, primitive_power);

    memset(second_multiplier, 0, MAX_DEGREE * sizeof(uint64_t));

    for(int i = 0; i < coset->size - 1; ++i){
        base_element *= field->characteristic;
        second_multiplier[1] = 1;
        second_multiplier[0] = construct_inverse_element_add(field, find_primitive_in_power(field, base_element));

        
        multiply_two_extended_polynomials(field, first_multiplier, second_multiplier);
        memset(second_multiplier, 0, MAX_DEGREE * sizeof(uint64_t));
    }
    memcpy(res, first_multiplier, sizeof(uint64_t) * MAX_DEGREE);
}

polynomial
extended_polynomial_to_polynomial(struct finite_field* field, uint64_t* poly)
{
    polynomial result = 0;
    for(int i = 0; i < MAX_DEGREE; ++i){
        if(field->characteristic == 2 && poly[i] == 1){
            result += (1UL << i);
        }
    }
    return result;
}

void
construct_generator_polynomial(struct finite_field* field, int number_of_errors, uint64_t res[MAX_DEGREE])
{
    struct cyclotomic_cosets* cosets = generate_cosets_from_n_needed_elements(field, 2 * number_of_errors);

    uint64_t first_multiplier[MAX_DEGREE], second_multiplier[MAX_DEGREE];
    construct_minimal_polynomial_from_coset(field, cosets->cosets[0], first_multiplier);
    for(int i = 1; i < cosets->number_of_cosets; ++i){
        construct_minimal_polynomial_from_coset(field, cosets->cosets[i], second_multiplier);
        multiply_two_extended_polynomials(field, first_multiplier, second_multiplier);
    }

    free_cosets_struct(cosets);

    memcpy(res, first_multiplier, sizeof(uint64_t) * MAX_DEGREE);
}

void
find_inverse_by_add_for_polynomial(struct finite_field* field, uint64_t a[MAX_DEGREE], uint64_t res[MAX_DEGREE])
{

    for (int i = 0; i < MAX_DEGREE; ++i) {
        res[i] = (needed_type)construct_inverse_element_add(field, a[i]);
    }
}

void
difference_of_two_polynomials(struct finite_field* field, uint64_t* a, uint64_t* b, uint64_t res[MAX_DEGREE])
{
    for(int i = 0; i < MAX_DEGREE; ++i){
        res[i] = add_in_field(field, a[i], construct_inverse_element_add(field, b[i]));
    }
}

void
divide_polynomials_with_remainder(struct finite_field* field, uint64_t dividend[],
    uint64_t divisor[], uint64_t result[MAX_DEGREE], uint64_t remainder[MAX_DEGREE])
{
    memset(result, 0, MAX_DEGREE * sizeof(uint64_t));
    memcpy(remainder, dividend, MAX_DEGREE * sizeof(uint64_t));

    int deg_divisor = get_degree(divisor);
    int deg_dividend = get_degree(dividend);

    if(deg_divisor < 0 || deg_dividend < deg_divisor){
        return;
    }
    for(int i = deg_dividend - deg_divisor; i > -1;--i){
        if(remainder[deg_divisor + i] == 0) continue;
        result[i] = 1;
        for(int j = 0; j <= deg_divisor; ++j){
            remainder[j + i] = add_in_field(field, remainder[j + i], divisor[j]);
        }
    }
}


uint64_t
find_value_from_root(struct finite_field* field, uint64_t poly[MAX_DEGREE], int primitive_power, int deg, uint64_t dir)
{
    uint64_t result = 0;
    int d = deg;
    for(int i = d; i >= 0; --i){
        result = multiply_in_field(field, result, dir);
        result = add_in_field(field, result, poly[i]);
    }
    return result;
}

int
equal_polynomials(uint64_t* a, uint64_t* b){
    if(a == NULL || b == NULL){
        return 0;
    }
    for(int i = 0; i < MAX_DEGREE; ++i){
        if(a[i] != b[i]) return 0;
    }

    return 1;
}


void
extended_euclidean_algorithm(struct finite_field* field, uint64_t* el, uint64_t* f, 
    uint64_t s0[MAX_DEGREE], uint64_t t0[MAX_DEGREE], uint64_t r0[MAX_DEGREE])
{
    int k = 0;
    copy_extended_polynomial(el, r0);
    memset(s0, 0, sizeof(uint64_t) * MAX_DEGREE);
    memset(t0, 0, sizeof(uint64_t) * MAX_DEGREE);
    

    uint64_t r1[MAX_DEGREE], s1[MAX_DEGREE], t1[MAX_DEGREE], q[MAX_DEGREE], r[MAX_DEGREE], temp[MAX_DEGREE], multiply_temp[MAX_DEGREE];
    copy_extended_polynomial(f, r1);
    memset(s1, 0, sizeof(uint64_t) * MAX_DEGREE);
    memset(t1, 0, sizeof(uint64_t) * MAX_DEGREE);
    s0[0] = 1;
    t1[0] = 1;

    while(get_degree(r1) >= 0){
        ++k;
        divide_polynomials_with_remainder(field, r0, r1, q, r);
        copy_extended_polynomial(r0, temp);
        copy_extended_polynomial(r1, r0);
        copy_extended_polynomial(q, multiply_temp);

        multiply_two_extended_polynomials(field, multiply_temp, r1);
        difference_of_two_polynomials(field, temp, multiply_temp, r1);

        
        copy_extended_polynomial(s0, temp);
        copy_extended_polynomial(s1, s0);
        copy_extended_polynomial(q, multiply_temp);
        
    

        multiply_two_extended_polynomials(field, multiply_temp, s1);
        difference_of_two_polynomials(field, temp, multiply_temp, s1);
        


        copy_extended_polynomial(t0, temp);
        copy_extended_polynomial(t1, t0);

        multiply_two_extended_polynomials(field, q, t1);


        difference_of_two_polynomials(field, temp, q, t1);
    }
}

uint64_t
construct_inverse_element_multiply(struct finite_field* field, uint64_t el)
{
    uint64_t t0[MAX_DEGREE], s0[MAX_DEGREE], r0[MAX_DEGREE];

    
    uint64_t el1[MAX_DEGREE], el2[MAX_DEGREE];
    construct_polynomial_from_field_element(field, get_primitive_polynomial(field), el1);
    construct_polynomial_from_field_element(field, el, el2);
    extended_euclidean_algorithm(field, el1, el2, s0, t0, r0);

    return get_field_el_from_polynomial(field, t0);
}

void
construct_polynomial_from_field_element(struct finite_field* field, uint64_t el, uint64_t res[MAX_DEGREE])
{
    memset(res, 0, sizeof(uint64_t) * MAX_DEGREE);
    for(int i = 0; el != 0; ++i){
        res[i] = el % field->characteristic;
        if(field->characteristic == 2){
            el >>= 1;
        }else{
            el /= field->characteristic;       
        }
    }
}

uint64_t
get_field_el_from_polynomial(struct finite_field* field, uint64_t poly[MAX_DEGREE])
{
    uint64_t result = 0;
    for(int i = 0; i < MAX_DEGREE; ++i){
        result += pow(field->characteristic, i) * poly[i];
    }
    return result;
}

void
multiply_polynomial_by_element(struct finite_field* field, uint64_t poly[MAX_DEGREE], uint64_t el)
{
    for(int i = 0; i < MAX_DEGREE; ++i){
        poly[i] = multiply_in_field(field, poly[i], el);
    }

}

int
get_degree(uint64_t poly[MAX_DEGREE])
{
    for(int i = MAX_DEGREE - 1; i >= 0; --i){
        if(poly[i] != 0){
            return i;
        }
    }
    return -1;
}
