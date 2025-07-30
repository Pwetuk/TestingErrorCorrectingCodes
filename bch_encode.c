#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bch_encode.h"


void
encode_bch(struct bch_code* bch_code_struct, uint64_t data[MAX_DEGREE], uint64_t result[MAX_DEGREE])
{
    uint64_t n = bch_code_struct->n - 1, tmp[MAX_DEGREE];
    memset(result, 0, sizeof(uint64_t) * MAX_DEGREE);

    memset(tmp, 0, MAX_DEGREE * sizeof(uint64_t));

    int deg_data = get_degree(data);
    
    for(int i = 0; i <= deg_data; ++i){
        tmp[n - bch_code_struct->data_length + 1 + i] = data[i];
    }
    uint64_t quotient[MAX_DEGREE], remainder[MAX_DEGREE];
    divide_polynomials_with_remainder(bch_code_struct->field, tmp, bch_code_struct->generator, quotient, remainder);
    
    difference_of_two_polynomials(bch_code_struct->field, tmp, remainder, result);
}   

void
decode_bch(struct bch_code* bch_code_struct, uint64_t message[MAX_DEGREE], uint64_t result[MAX_DEGREE])
{
    uint64_t locator[MAX_DEGREE];
    int errors_pos[MAX_DEGREE];
    uint64_t syndrome[MAX_DEGREE];
    memset(syndrome, 0, MAX_DEGREE * sizeof(uint64_t));
    int is_null = calculate_syndrome(bch_code_struct->field, message, bch_code_struct->number_of_errors, syndrome);



    if(is_null == 1){
        construct_locator_polynomial(bch_code_struct->field, syndrome, bch_code_struct->number_of_errors, locator);
        chien_search(bch_code_struct->field, locator, bch_code_struct->number_of_errors, errors_pos, bch_code_struct->n);
        
        fix_errors_in_bch_binary(message, errors_pos, bch_code_struct->number_of_errors);
    }

    memset(result, 0, MAX_DEGREE * sizeof(uint64_t));
    
    int result_deg = bch_code_struct->data_length;
    int n = bch_code_struct->n - 1;
    int messg_degree = get_degree(message);
    for(int i = 0; i < bch_code_struct->data_length; ++i){
        result[i] = message[n - result_deg + 1 + i];
    }
}

int
calculate_syndrome(struct finite_field* field, uint64_t encoded_message[MAX_DEGREE], int number_of_errors, uint64_t syndrome[MAX_DEGREE])
{
    int number_of_not_null = 0;
    int d = get_degree(encoded_message);
    uint64_t need_to_check = 2;
    for(int i = 0; i < 2 * number_of_errors; ++i){
        syndrome[i] = find_value_from_root(field, encoded_message, d, need_to_check);
        if(syndrome[i] != 0){
            number_of_not_null = 1;
        }
        need_to_check = multiply_in_field(field, need_to_check, 2);
    }
    return number_of_not_null;
}



void
construct_locator_polynomial(struct finite_field* field, uint64_t* syndrome, int number_of_errors, uint64_t locator[MAX_DEGREE]){
    uint64_t B[MAX_DEGREE], locator_old[MAX_DEGREE], B_copy[MAX_DEGREE], temp[MAX_DEGREE];

    memset(B, 0, MAX_DEGREE * sizeof(uint64_t));
    memset(B_copy, 0, MAX_DEGREE * sizeof(uint64_t));

    memset(locator, 0, sizeof(uint64_t) * MAX_DEGREE);

    locator[0] = 1;
    B[0] = 1;
    int L = 0, m = 1, deg_b;
    uint64_t delta_r, b = 1;
    for(int r = 0; r < 2 * number_of_errors; ++r){
        
        delta_r = syndrome[r];
        for(int i = 1; i <= L && i < MAX_DEGREE; ++i){
            delta_r = add_in_field(field, delta_r, multiply_in_field(field, syndrome[r - i], locator[i]));
        }        

        if(delta_r == 0){
            ++m;
        }else{
            memcpy(locator_old, locator, MAX_DEGREE * sizeof(uint64_t));

            deg_b = get_degree(B);

            memset(B_copy, 0, MAX_DEGREE * sizeof(uint64_t));

            for(int i = 0; i <= deg_b; ++i){
                B_copy[i + m] = B[i];
            }
            uint64_t inversed = construct_inverse_element_multiply(field, b);
            uint64_t coef = multiply_in_field(field, delta_r, inversed);
            multiply_polynomial_by_element(field, B_copy, coef);
            difference_of_two_polynomials(field, locator, B_copy, temp);
            memcpy(locator, temp, sizeof(uint64_t) * MAX_DEGREE);

            if(2 * L <= r){
                L = r + 1 - L;
                memcpy(B, locator_old, sizeof(uint64_t) * MAX_DEGREE);
                b = delta_r;
                m = 1;  
            }else{
                ++m;
            }
            
        }
    }
}

void
chien_search(struct finite_field* field, uint64_t* locator, int number_of_errors, int* result, uint64_t n)
{
    uint64_t pos;
    uint64_t need_to_check = 1;
    uint64_t inversed_primitive = construct_inverse_element_multiply(field, 2);
    int err_pos = 0;
    int deg = get_degree(locator);


    memset(result, -1, sizeof(int) * MAX_DEGREE);

    for(uint64_t i = 0; i < n - 1; ++i){

        pos = find_value_from_root(field, locator, deg, need_to_check);
        if(pos == 0){
            result[err_pos] = i;
            ++err_pos;
            if(err_pos >= number_of_errors){
                break;
            }
        }
        need_to_check = multiply_in_field(field, need_to_check, inversed_primitive);
    }
}

void
fix_errors_in_bch_binary(uint64_t* codeword, int* errors_pos, int number_of_errors)
{
    for(int i = 0; i < number_of_errors; ++i){
        if(errors_pos[i] == -1){
            break;
        }
        codeword[errors_pos[i]] ^= 1;
    }
}

struct bch_code*
init_bch(uint64_t p, uint64_t power, int number_of_errors, uint64_t primitive)
{
    struct bch_code* to_return = malloc(sizeof(struct bch_code));

    if(primitive != 0){
        to_return->field = malloc(sizeof(struct finite_field));
        to_return->field->characteristic = p;
        to_return->field->power = power;
        to_return->field->primative_in_power_n = primitive;
    }else{
        to_return->field = construct_gf_p_m(p, power);
    }

    to_return->number_of_errors = number_of_errors;
    construct_generator_polynomial(to_return->field, number_of_errors, to_return->generator);
    to_return->n = pow(to_return->field->characteristic, to_return->field->power);
    to_return->data_length = to_return->n - get_degree(to_return->generator) - 1;
    return to_return;
}

void
free_bch_code(struct bch_code* bch_code_struct)
{
    free(bch_code_struct->field);
    free(bch_code_struct);
}
