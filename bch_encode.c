#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "bch_encode.h"


struct extended_polynomial*
encode_bch(struct finite_field* field, struct extended_polynomial* generator, struct extended_polynomial* data, int n, int degree_data)
{
    struct extended_polynomial* result;
    data = multiply_extended_polynomial_by_x_n(data, n - degree_data, 0);
    struct extended_polynomial *quotient = NULL, *remainder = NULL;
    divide_polynomials_with_remainder(field, data, generator, &quotient, &remainder);
    result = difference_of_two_polynomials(field, data, remainder);
    if(quotient == NULL){
        printf("Something went horribly wrong\n");
    }
    free_extended_polynomial(quotient);
    free_extended_polynomial(remainder);

    free_extended_polynomial(data);

    return result;
}

unsigned long long
decode_bch(struct finite_field* field, struct extended_polynomial* message, int number_of_errors, int original_data_length)
{
    long long unsigned int* syndrome = calculate_syndrome(field, message, number_of_errors);

    struct extended_polynomial* locator = construct_locator_polynomial(field, syndrome, number_of_errors);

    int* errors_pos = chien_search(field, locator, number_of_errors);

    fix_errors_in_bch_binary(message, errors_pos, number_of_errors);

    unsigned long long data = 0;
    for(int i = 0; i < original_data_length; ++i){
        data <<= 1;
        data += message->coefs[message->degree - i - 1];
    }

    free(syndrome);
    free(errors_pos);
    free_extended_polynomial(locator);
    return data;
}

unsigned long long*
calculate_syndrome(struct finite_field* field, struct extended_polynomial* encoded_message, int number_of_errors)
{
    unsigned long long *syndrome = malloc(sizeof(unsigned long long) * 2 * number_of_errors);
    for(int i = 0; i < 2 * number_of_errors; ++i){
        syndrome[i] = find_value_from_root(field, encoded_message, i + 1);
    }
    return syndrome;
}



struct extended_polynomial*
construct_locator_polynomial(struct finite_field* field, unsigned long long* syndrome, int number_of_errors){
    int base_coef[2] = {1, 0};
    struct extended_polynomial* locator = construct_extended_polynomial_from_coefs(base_coef, 1);
    struct extended_polynomial* B = construct_extended_polynomial_from_coefs(base_coef, 1);
    struct extended_polynomial* locator_old, *B_copy;
    int L = 0, m = 1, b = 1;
    unsigned long long delta_r;
    for(int r = 0; r < 2 * number_of_errors; ++r){
        
        delta_r = syndrome[r];
        for(int i = 1; i <= L; ++i){
            delta_r = add_in_field(field, delta_r, multiply_in_field(field, syndrome[r - i], locator->coefs[i]));
        }

        
        
        if(delta_r == 0){
            ++m;
        }else{
            locator_old = locator;
            B_copy = multiply_extended_polynomial_by_x_n(B, m, 0);
            unsigned long long coef = multiply_in_field(field, delta_r, construct_inverse_element_multiply(field, b));
            B_copy = multiply_polynomial_by_element(field, B_copy, coef, NEED_FREE);
            locator = difference_of_two_polynomials(field, 
                locator, B_copy);
                free_extended_polynomial(B_copy);
                if(2 * L <= r){
                    L = r + 1 - L;
                    free_extended_polynomial(B);
                    B = copy_extended_polynomial(locator_old);
                    b = delta_r;
                    m = 1;
                }else{
                    ++m;
                }
                free_extended_polynomial(locator_old);
        }
        //printf("deltar_r: %lld\n", delta_r);
        //printf("Locator is: \t");
        //print_extended_polynomial(locator);
        //printf("B is: \t");
        //print_extended_polynomial(B);
    }

    free_extended_polynomial(B);

    return locator;
}

int*
chien_search(struct finite_field* field, struct extended_polynomial* locator, int number_of_errors)
{
    unsigned long long n = pow(field->characteristic, field->power);
    unsigned long long pos;
    int err_pos = 0;
    int *result = malloc(number_of_errors * sizeof(int));
    for(int i = 0; i < number_of_errors; ++i){
        result[i] = -1;
    }
    for(unsigned long long i = 0; i < n - 1; ++i){
        pos = find_value_from_root(field, locator, -i);
        if(pos == 0){
            result[err_pos] = i;
            ++err_pos;
            if(err_pos >= number_of_errors){
                break;
            }
        }
    }
    return result;
}

void
fix_errors_in_bch_binary(struct extended_polynomial* codeword, int* errors_pos, int number_of_errors)
{
    for(int i = 0; i < number_of_errors; ++i){
        if(errors_pos[i] == -1){
            break;
        }
        codeword->coefs[errors_pos[i]] ^= 1;
    }
}