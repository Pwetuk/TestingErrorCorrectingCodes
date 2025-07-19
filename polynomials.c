#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "polynomials.h"


struct cyclotomic_cosets*
generate_cosets_from_n_needed_elements(struct finite_field* field, int n)
{
    int was = 0;
    struct cyclotomic_coset** cosets_from = malloc(sizeof(struct cyclotomic_coset*) * n);
    int num_of_elemets = 0;
    unsigned long long int p = pow(field->characteristic, field->power) - 1;
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
    //printf("Freeing %lld\n", to_free->coefs);
    free(to_free->coefs);
    //printf("After coefs\n");
    free(to_free);
}

struct extended_polynomial*
construct_extended_polynomial_from_coefs(int* coefs, int n){
    struct extended_polynomial* result = malloc(sizeof(struct extended_polynomial));
    result->coefs = malloc(sizeof(unsigned long long int) * n);
    for(int i = 0; i < n; ++i){
        result->coefs[i] = coefs[i];
    }
    result->degree = n;
    return result;
}

void
print_extended_polynomial(struct extended_polynomial* to_print)
{
    if(to_print == NULL){
        printf("Poly is NULL\n");
        return;
    }
    printf("Degree: %d\t", to_print->degree);
    if(to_print->degree == 0) {
        printf("\n\n");
        return;
    }
    for(int i = 0; i < to_print->degree; ++i){
        if(to_print->coefs[i] != 0){
            printf("%lld * x^%d + ", to_print->coefs[i], i);   
        }
    }
    printf("\n\n");
}

struct extended_polynomial*
copy_extended_polynomial(struct extended_polynomial* poly)
{
    struct extended_polynomial* result = malloc(sizeof(struct extended_polynomial));
    result->degree = poly->degree;
    result->coefs = malloc(sizeof(unsigned long long int) * (result->degree + 1));
    for(int i = 0; i < poly->degree; ++i){
        result->coefs[i] = poly->coefs[i];
    }
    return result;
}

struct extended_polynomial*
add_two_extended_polynomials(struct finite_field* field, struct extended_polynomial* poly1, struct extended_polynomial* poly2)
{
    int total_degree = max(poly1->degree, poly2->degree);
    struct extended_polynomial* result = malloc(sizeof(struct extended_polynomial));
    int last_non_zero = 0; 
    result->coefs = malloc(sizeof(unsigned long long) * total_degree);

    for(int i = 0; i < max(poly1->degree, poly2->degree); ++i){
        if(poly1->degree < i){
            result->coefs[i] = poly2->coefs[i];
            last_non_zero = i + 1;
        }else if(poly2->degree < i){
            result->coefs[i] = poly1->coefs[i];
            last_non_zero = i + 1;
        }else{
            result->coefs[i] = add_in_field(field, poly1->coefs[i], poly2->coefs[i]);
            if(result->coefs[i] != 0){
                last_non_zero = i + 1;
            }
        }
        
    }
    result->degree = last_non_zero;
    return result;
}

struct extended_polynomial*
multiply_two_extended_polynomials(struct finite_field* field, struct extended_polynomial* poly1, struct extended_polynomial* poly2)
{
    //printf("We are now multiplying\n");
    //print_extended_polynomial(poly1);
    //print_extended_polynomial(poly2);
    if(poly1->degree == 0 || poly2->degree == 0) return make_zero_polynomial(1);
    int total_degree = poly1->degree + poly2->degree - 1;
    struct extended_polynomial* result = malloc(sizeof(struct extended_polynomial));
    result->degree = total_degree;
    result->coefs = calloc(total_degree, sizeof(unsigned long long int));
    for(int i = 0; i < poly1->degree; ++i){
        for(int j = 0; j < poly2->degree; ++j){
            result->coefs[i + j] = add_in_field(field, multiply_in_field(field, poly1->coefs[i], poly2->coefs[j]), result->coefs[i + j]);
        }
    }

    //print_extended_polynomial(result);
    //printf("--------------\n\n");
    return result;
}


struct extended_polynomial*
construct_minimal_polynomial_from_coset(struct finite_field* field, struct cyclotomic_coset* coset)
{
    int base_element = coset->from_el;

    struct extended_polynomial* first_multiplier = malloc(sizeof(struct extended_polynomial));
    struct extended_polynomial* second_multiplier = malloc(sizeof(struct extended_polynomial));
    struct extended_polynomial* temp_polynomial;
    first_multiplier->degree = 2;
    first_multiplier->coefs = malloc(sizeof(unsigned long long int) * 2);
    first_multiplier->coefs[1] = 1;
    first_multiplier->coefs[0] = construct_inverse_element_add(field, find_primitive_in_power(field, base_element));

    second_multiplier->degree = 2;
    second_multiplier->coefs = malloc(sizeof(unsigned long long int) * 2);

    for(int i = 0; i < coset->size - 1; ++i){
        base_element *= field->characteristic;
        second_multiplier->coefs[1] = 1;
        second_multiplier->coefs[0] = construct_inverse_element_add(field, find_primitive_in_power(field, base_element));

        temp_polynomial = multiply_two_extended_polynomials(field, first_multiplier, second_multiplier);
        free_extended_polynomial(first_multiplier);
        
        first_multiplier = temp_polynomial;
        
    }

    free_extended_polynomial(second_multiplier);

    return first_multiplier;
}

polynomial
extended_polynomial_to_polynomial(struct finite_field* field, struct extended_polynomial* poly)
{
    polynomial result = 0;
    for(int i = 0; i < poly->degree; ++i){
        if(field->characteristic == 2 && poly->coefs[i] == 1){
            result += (1UL << i);
        }else if(field->characteristic != 2){
            result += poly->coefs[i] * pow(field->characteristic, i);
        }
    }
    return result;
}

struct extended_polynomial*
construct_generator_polynomial(struct finite_field* field, int number_of_errors)
{
    struct cyclotomic_cosets* cosets = generate_cosets_from_n_needed_elements(field, 2 * number_of_errors);

    struct extended_polynomial* first_multiplier, *second_multiplier;
    first_multiplier = construct_minimal_polynomial_from_coset(field, cosets->cosets[0]);
    for(int i = 1; i < cosets->number_of_cosets; ++i){
        second_multiplier = construct_minimal_polynomial_from_coset(field, cosets->cosets[i]);
        first_multiplier = multiply_two_extended_polynomials(field, first_multiplier, second_multiplier);
    }

    free_cosets_struct(cosets);


    free_extended_polynomial(second_multiplier);

    return first_multiplier;
}

struct extended_polynomial*
find_inverse_by_add_for_polynomial(struct finite_field* field, struct extended_polynomial* a)
{
    struct extended_polynomial* result = malloc(sizeof(struct extended_polynomial));

    result->degree = a->degree;
    result->coefs = malloc(sizeof(unsigned long long) * (a->degree + 2));
    if (result->coefs == NULL) {
        printf("Something went horribly wrong\n");
        free(result);
        return NULL;
    }

    for (int i = 0; i < a->degree; ++i) {
        unsigned long long temp = construct_inverse_element_add(field, a->coefs[i]);
        result->coefs[i] = temp;
    }
    return result;
}

struct extended_polynomial*
difference_of_two_polynomials(struct finite_field* field, struct extended_polynomial* a, struct extended_polynomial* b)
{
    //printf("We are now permoning difference\n");
    //print_extended_polynomial(a);
    //print_extended_polynomial(b);
    int i = 0;
    int result_degree = 0;
    struct extended_polynomial* result = malloc(sizeof(struct extended_polynomial));
    result->coefs = malloc(sizeof(unsigned long long) * max(a->degree, b->degree));
    for(; i < min(a->degree, b->degree); ++i){
        result->coefs[i] = add_in_field(field, a->coefs[i], construct_inverse_element_add(field, b->coefs[i]));
        if(result->coefs[i] != 0) result_degree = i + 1;
    }
    if(i < b->degree){
        for(; i < b->degree; ++i){
            result->coefs[i] = construct_inverse_element_add(field, b->coefs[i]);
        }
        result_degree = b->degree;
    }else if(i < a->degree){
        for(; i < a->degree; ++i){
            result->coefs[i] = a->coefs[i];
        }
        result_degree = a->degree;
    }
    result->degree = result_degree;
    //print_extended_polynomial(result);
    //printf("--------------\n\n");
    return result;
}

void
divide_polynomials_with_remainder(struct finite_field* field, struct extended_polynomial* dividing,
     struct extended_polynomial* divider, struct extended_polynomial** result, struct extended_polynomial** remainder)
{
    int deg_a = dividing->degree - 1;
    int deg_b = divider->degree - 1;
    if(deg_a < deg_b){
        *remainder = copy_extended_polynomial(dividing);
        *result = make_zero_polynomial(1);
    }
    struct extended_polynomial* inverse_poly = find_inverse_by_add_for_polynomial(field, divider);
    *result = make_zero_polynomial(deg_a - deg_b + 1);
    (*result)->degree = deg_a - deg_b + 1;
    unsigned long long int new_coef = 0;
    while(deg_a >= deg_b){
        struct extended_polynomial* shifted_inverse = multiply_extended_polynomial_by_x_n(inverse_poly, deg_a - deg_b, 0);
        new_coef = 0;

        while(dividing->degree - 1 >= deg_a){
            if(new_coef > field->characteristic){
                print_extended_polynomial(dividing);
                print_extended_polynomial(divider);
                print_extended_polynomial(shifted_inverse);
                printf("ASdasdasdas\n");
                *remainder = NULL;
                return;
            }
            dividing = add_two_extended_polynomials(field, dividing, shifted_inverse);
            ++new_coef;
        }
        (*result)->coefs[deg_a - deg_b] = new_coef;
        //printf("Result is: ");
        //print_extended_polynomial(*result);
        deg_a = dividing->degree - 1;

        free_extended_polynomial(shifted_inverse);
    }
    free_extended_polynomial(inverse_poly);
    *remainder = copy_extended_polynomial(dividing);
}


unsigned long long
find_value_from_root(struct finite_field* field, struct extended_polynomial* poly, int primitive_power)
{
    //printf(
    //    "INserting: %lld", find_primitive_in_power(field, primitive_power)
    //);
    unsigned long long result = 0;
    for(int i = 0; i < poly->degree; ++i){
        result = add_in_field(field, 
            multiply_in_field(field, poly->coefs[i], find_primitive_in_power(field, primitive_power * i)),
            result);
        //printf("Temp res: %lld\t", result);
    }
    return result;
}


struct extended_polynomial*
multiply_extended_polynomial_by_x_n(struct extended_polynomial* poly, int n, int need_free){
    struct extended_polynomial* result = malloc(sizeof(struct extended_polynomial));
    unsigned long long* new_coefs = calloc(poly->degree + n, sizeof(unsigned long long));
    for(int i = 0; i < poly->degree; ++i){
        new_coefs[i + n] = poly->coefs[i];
    }
    result->coefs = new_coefs;
    result->degree = poly->degree + n;
    //printf("At shift with n=%d: ", n);
    //print_extended_polynomial(result);
    if(need_free == 1) free_extended_polynomial(poly);
    return result;
}

int
equal_polynomials(struct extended_polynomial* a, struct extended_polynomial* b){
    if(a == NULL || b == NULL){
        return 0;
    }
    if(a->degree != b->degree) return 0;
    for(int i = 0; i < a->degree; ++i){
        if(a->coefs[i] != b->coefs[i]) return 0;
    }

    return 1;
}


struct extended_polynomial*
make_zero_polynomial(int degree){
    struct extended_polynomial* result = malloc(sizeof(struct extended_polynomial));
    result->degree = 0;
    result->coefs = calloc(degree, sizeof(unsigned long long));
    return result;
}

void
extended_euclidean_algorithm(struct finite_field* field, struct extended_polynomial* el, struct extended_polynomial* f, 
    struct extended_polynomial** s0, struct extended_polynomial** t0, struct extended_polynomial** r0)
{
    *r0 = copy_extended_polynomial(el);
    *s0 = make_zero_polynomial(1);
    *t0 = make_zero_polynomial(1); 
    
    struct extended_polynomial* r1 = copy_extended_polynomial(f);
    struct extended_polynomial* s1 = make_zero_polynomial(1);
    struct extended_polynomial* t1 = make_zero_polynomial(1);
    struct extended_polynomial* q, *r, *temp, *multiply_temp;
    (*s0)->coefs[0] = 1;
    t1->coefs[0] = 1;

    (*s0)->degree = 1;
    t1->degree = 1;

    while(r1->degree > 0){
        divide_polynomials_with_remainder(field, (*r0), r1, &q, &r);
        free_extended_polynomial(r);

        temp = (*r0);
        (*r0) = r1;
        multiply_temp = multiply_two_extended_polynomials(field, q, r1);
        r1 = difference_of_two_polynomials(field, temp, multiply_temp);
        
        free_extended_polynomial(multiply_temp);
        free_extended_polynomial(temp);


        temp = (*s0);
        (*s0) = s1;
        multiply_temp = multiply_two_extended_polynomials(field, q, s1);
        s1 = difference_of_two_polynomials(field, temp, multiply_temp);
        
        free_extended_polynomial(multiply_temp);
        free_extended_polynomial(temp);

        temp = (*t0);
        (*t0) = t1;
        multiply_temp = multiply_two_extended_polynomials(field, q, t1);
        t1 = difference_of_two_polynomials(field, temp, multiply_temp);
        
        free_extended_polynomial(multiply_temp);
        free_extended_polynomial(temp);
        
        free_extended_polynomial(q);

    }
    free_extended_polynomial(r1);
    free_extended_polynomial(s1);
    free_extended_polynomial(t1);
}

unsigned long long
construct_inverse_element_multiply(struct finite_field* field, unsigned long long el)
{
    struct extended_polynomial* t0, *s0, *r0;

    extended_euclidean_algorithm(field, construct_polynomial_from_field_element(field, get_primitive_polynomial(field)),
        construct_polynomial_from_field_element(field, el), &s0, &t0, &r0);

    free_extended_polynomial(s0);
    free_extended_polynomial(r0);

    return get_field_el_from_polynomial(field, t0, NEED_FREE);
}

struct extended_polynomial*
construct_polynomial_from_field_element(struct finite_field* field, unsigned long long el)
{
    struct extended_polynomial* result = malloc(sizeof(struct extended_polynomial));
    unsigned long long temp_el = el;
    for(int i = 0; el != 0; ++i){
        if(field->characteristic == 2){
            el >>= 1;
        }else{
            el /= field->characteristic;       
        }
        result->degree = i + 1;
    }

    el = temp_el;

    result->coefs = calloc(result->degree + 1, sizeof(unsigned long long));
    for(int i = 0; el != 0; ++i){
        result->coefs[i] = el % field->characteristic;
        if(field->characteristic == 2){
            el >>= 1;
        }else{
            el /= field->characteristic;       
        }
    }

    return result;
}

unsigned long long
get_field_el_from_polynomial(struct finite_field* field, struct extended_polynomial* poly, int need_free)
{
    unsigned long long int result = 0;
    for(int i = 0; i < poly->degree; ++i){
        result += pow(field->characteristic, i) * poly->coefs[i];
    }
    if(need_free == 1) free_extended_polynomial(poly);
    return result;
}

struct extended_polynomial*
multiply_polynomial_by_element(struct finite_field* field, struct extended_polynomial* poly, unsigned long long int el, int need_free)
{
    struct extended_polynomial* result = copy_extended_polynomial(poly);
    int initial = result->degree;

    for(int i = 0; i < initial; ++i){
        result->coefs[i] = multiply_in_field(field, result->coefs[i], el);
        if(result->coefs[i] != 0) result->degree = i + 1;
    }

    if(need_free == NEED_FREE) free_extended_polynomial(poly);
    return result;
}