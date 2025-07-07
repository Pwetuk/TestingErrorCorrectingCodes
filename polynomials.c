#include <stdlib.h>
#include <stdio.h>

#include "polynomials.h"


struct cyclotomic_cosets*
generate_cosets_from_n_needed_elements(struct finite_field* field, int n)
{
    int was = 0;
    struct cyclotomic_coset** cosets_from = malloc(sizeof(struct cyclotomic_coset*) * n);
    int num_of_elemets = 0;
    for(int i = 1; i <= n; ++i){
        if((was & (1 << i)) != 0){
            continue;
        }
        cosets_from[num_of_elemets] = malloc(sizeof(struct cyclotomic_coset));
        cosets_from[num_of_elemets]->from_el = i;
        cosets_from[num_of_elemets]->size = 1;
        int k = i;
        was |= (1 << i);
        k <<= 1;
        printf("%d, ", k);
        while(k != i){
            was |= (1 << k);
            k <<= 1;
            k %= ((1 << field->power) - 1);
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



polynomial
add_two_polynomials(struct finite_field* field, polynomial poly1, polynomial poly2)
{
    polynomial to_return;
    to_return = add_in_field(field, poly1, poly2);
    return to_return;
}


polynomial
multiply_two_polynomials(struct finite_field* field, polynomial poly1, polynomial poly2)
{
    polynomial to_return = 0;
    unsigned long long int backup_poly = poly1;
    int current_power = 0;

    while(backup_poly > 0){
        if(backup_poly % 2 != 0){
            to_return = add_two_polynomials(field, to_return, (poly2 << current_power));
        }
        backup_poly >>= 1;
        current_power += 1;
    }
    return to_return;
}


void
print_polynomial(struct finite_field* field, polynomial poly)
{
    unsigned long long out_poly = poly;
    for(int i = 0; out_poly != 0; ++i){
        if(out_poly % 2 != 0){
            printf("x^%d + ", i);
        }
        out_poly >>= 1;
    }
    printf("\n");
}


void
free_extended_polynomial(struct extended_polynomial* to_free)
{
    free(to_free->coefs);
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
    result->coefs = malloc(sizeof(unsigned long long int) * result->degree);
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
            last_non_zero = i;
        }else if(poly2->degree < i){
            result->coefs[i] = poly1->coefs[i];
            last_non_zero = i;
        }else{
            result->coefs[i] = add_in_field(field, poly1->coefs[i], poly2->coefs[i]);
            if(result->coefs[i] != 0){
                last_non_zero = i;
            }
        }
        
    }
    ++last_non_zero;
    result->degree = last_non_zero;
    return result;
}

struct extended_polynomial*
multiplty_two_extended_polynomials(struct finite_field* field, struct extended_polynomial* poly1, struct extended_polynomial* poly2)
{
    int total_degree = poly1->degree + poly2->degree;
    struct extended_polynomial* result = malloc(sizeof(struct extended_polynomial));
    result->degree = total_degree;
    result->coefs = calloc(total_degree, sizeof(unsigned long long int));
    for(int i = 0; i < poly1->degree; ++i){
        for(int j = 0; j < poly2->degree; ++j){
            result->coefs[i + j] = add_in_field(field, multiply_in_field(field, poly1->coefs[i], poly2->coefs[j]), result->coefs[i + j]);
        }
    }


    return result;
}


polynomial
construct_minimal_polynomial_from_coset(struct finite_field* field, struct cyclotomic_coset* coset)
{
    int base_element = coset->from_el;

    struct extended_polynomial* first_multiplier = malloc(sizeof(struct extended_polynomial));
    struct extended_polynomial* second_multiplier = malloc(sizeof(struct extended_polynomial));
    struct extended_polynomial* temp_polynomial;
    first_multiplier->degree = 2;
    first_multiplier->coefs = malloc(sizeof(unsigned long long int) * 2);
    first_multiplier->coefs[1] = 1;
    first_multiplier->coefs[0] = find_primitive_in_power(field, base_element);

    
    second_multiplier->degree = 2;
    second_multiplier->coefs = malloc(sizeof(unsigned long long int) * 2);


    for(int i = 0; i < coset->size - 1; ++i){
        base_element <<= 1;
        second_multiplier->coefs[1] = 1;
        second_multiplier->coefs[0] = find_primitive_in_power(field, base_element);

        temp_polynomial = multiplty_two_extended_polynomials(field, first_multiplier, second_multiplier);
        free_extended_polynomial(first_multiplier);
        
        first_multiplier = temp_polynomial;
        
    }

    
    polynomial result = extended_polynomial_to_polynomial(field, first_multiplier);
    
    free_extended_polynomial(first_multiplier);
    free_extended_polynomial(second_multiplier);

    return result;
}

polynomial
extended_polynomial_to_polynomial(struct finite_field* field, struct extended_polynomial* poly)
{
    polynomial result = 0;
    for(int i = 0; i < poly->degree; ++i){
        switch (poly->coefs[i])
        {
        case 1:
            result += (1UL << i);
            break;
        case 0:
            break;
        default:
            return -1;
            break;
        }
    }
    return result;
}

polynomial
construct_generator_polynomial(struct finite_field* field, int number_of_errors)
{
    struct cyclotomic_cosets* cosets = generate_cosets_from_n_needed_elements(field, 2 * number_of_errors);

    polynomial first_multiplier, second_multiplier;
    first_multiplier = construct_minimal_polynomial_from_coset(field, cosets->cosets[0]);
    for(int i = 1; i < cosets->number_of_cosets; ++i){
        second_multiplier = construct_minimal_polynomial_from_coset(field, cosets->cosets[i]);
        //printf("Minimal polynomialials: %lld, %lld\n", first_multiplier, second_multiplier);
        first_multiplier = multiply_two_polynomials(field, first_multiplier, second_multiplier);
    }

    free_cosets_struct(cosets);

    return first_multiplier;
}