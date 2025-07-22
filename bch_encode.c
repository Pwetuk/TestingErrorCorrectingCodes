#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bch_encode.h"


void
encode_bch(struct bch_code* bch_code_struct, uint64_t data[MAX_DEGREE], uint64_t result[MAX_DEGREE])
{
    uint64_t n = 1 << bch_code_struct->field->power;
    memset(result, 0, sizeof(uint64_t) * MAX_DEGREE);

    int deg_data = get_degree(data);
    for(int i = 0; i < deg_data; ++i){
        data[n - bch_code_struct->data_length - 1 + i] = data[i];
        data[i] = 0;
    }

    uint64_t quotient[MAX_DEGREE], remainder[MAX_DEGREE];
    divide_polynomials_with_remainder(bch_code_struct->field, data, bch_code_struct->generator, quotient, remainder);
    difference_of_two_polynomials(bch_code_struct->field, data, remainder, result);
}

void
decode_bch(struct bch_code* bch_code_struct, uint64_t message[MAX_DEGREE], uint64_t result[MAX_DEGREE])
{
    uint64_t locator[MAX_DEGREE];
    int errors_pos[MAX_DEGREE];
    uint64_t syndrome[MAX_DEGREE];
    calculate_syndrome(bch_code_struct->field, message, bch_code_struct->number_of_errors, syndrome);
    printf("Calculated\n");
    print_extended_polynomial(syndrome);

    uint8_t is_null = 0;
    for(int i = 0; i < bch_code_struct->number_of_errors; ++i){
        if(syndrome[i] != 0){
            is_null = 1;
            break;
        }
    }

    if(is_null == 1){
        construct_locator_polynomial(bch_code_struct->field, syndrome, bch_code_struct->number_of_errors, locator);
        print_extended_polynomial(locator);
        chien_search(bch_code_struct->field, locator, bch_code_struct->number_of_errors, errors_pos);
        printf("Errors: ");
        for(int i = 0; i <= bch_code_struct->number_of_errors; ++i){
            printf("%d\t", errors_pos[i]);
        }
        fix_errors_in_bch_binary(message, errors_pos, bch_code_struct->number_of_errors);
    }

    memset(result, 0, MAX_DEGREE * sizeof(uint64_t));
    
    int result_deg = bch_code_struct->data_length;
    int messg_degree = get_degree(message);

    for(int i = 0; i < bch_code_struct->data_length; ++i){
        result[i] = message[messg_degree - result_deg + i];
    }

}

void
calculate_syndrome(struct finite_field* field, uint64_t encoded_message[MAX_DEGREE], int number_of_errors, uint64_t syndrome[MAX_DEGREE])
{
    for(int i = 0; i < 2 * number_of_errors; ++i){
        printf("i: %d\n", i);
        syndrome[i] = find_value_from_root(field, encoded_message, i + 1);
    }
    printf("Calculated 2\n");
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
chien_search(struct finite_field* field, uint64_t* locator, int number_of_errors, int* result)
{

    uint64_t n = 1 << field->power;
    uint64_t pos;
    int err_pos = 0;

    memset(result, -1, sizeof(int) * MAX_DEGREE);

    for(uint64_t i = 0; i < n - 1; ++i){
        pos = find_value_from_root(field, locator, -i);
        if(pos == 0){
            result[err_pos] = i;
            ++err_pos;
            if(err_pos >= number_of_errors){
                break;
            }
        }
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
    to_return->data_length = pow(to_return->field->characteristic, to_return->field->power) - get_degree(to_return->generator);
    return to_return;
}

void
free_bch_code(struct bch_code* bch_code_struct)
{
    free(bch_code_struct->field);
    free(bch_code_struct);
}

/*
uint64_t*
bch_encode_arr(struct bch_code* bch, uint64_t* to_encode_arr, uint64_t to_encode_length, uint64_t *encoded_length)
{
    uint64_t el;
    int ind = 0, threshold = 64, last_power = 0, not_encoded = 1;
    *encoded_length = 0;
    uint64_t size = 64;

    uint64_t *result = calloc(size, sizeof(uint64_t)), *temp;

    struct extended_polynomial* encode_poly = make_zero_polynomial(bch->data_length), *encoded;
    encode_poly->degree = bch->data_length;

    for(uint64_t i = 0; i < to_encode_length; ++i){
        el = to_encode_arr[i];
        for(int k = 0; (k < threshold) || ((ind != 0) && (i + 1 == to_encode_length)); ++k){
            not_encoded = 1;
            encode_poly->coefs[ind] = el % bch->field->characteristic;
            el >>= 1;
            ++ind;
            if(ind >= bch->data_length){
                not_encoded = 0;
                encoded = encode_bch(bch, encode_poly);
                
                for(int j = 0; j < encoded->degree; ++j){
                    result[*encoded_length] += (1ULL << last_power) * encoded->coefs[j];
                    ++last_power;
                    if(last_power >= threshold){
                        last_power = 0;
                        ++(*encoded_length);
                    }
                    if(*encoded_length >= size){
                        temp = calloc(size * 2, sizeof(uint64_t));
                        if(temp == NULL){
                            printf("NOT ALLOC'D\n");
                            return NULL;
                        }
                        memcpy(temp, result, size * sizeof(uint64_t));
                        size *= 2;
                        free(result);
                        result = temp;
                    }
                }
                free_extended_polynomial(encoded);
                ind = 0;
            }
        }
    }
    *encoded_length += 1;
    free_extended_polynomial(encode_poly);
    return result;
}

uint64_t*
bch_decode_arr(struct bch_code* bch, uint64_t* to_decode_arr, uint64_t to_decode_length, uint64_t *decoded_length)
{
    uint64_t el, last_power = 0;
    int threshold = 64;
    int ind = 0;
    *decoded_length = 0;
    uint64_t size = 64;
    int not_decoded = 1;

    uint64_t *result = calloc(size, sizeof(uint64_t)), *temp;

    struct extended_polynomial* decode_poly = make_zero_polynomial(bch->data_length + bch->generator->degree - 1), *decoded;
    decode_poly->degree = bch->data_length + bch->generator->degree - 1;

    for(uint64_t i = 0; i < to_decode_length; ++i){
        el = to_decode_arr[i];
        for(int k = 0; (k < threshold) || ((ind != 0) && (i + 1 == to_decode_length)); ++k){
            not_decoded = 1;
            decode_poly->coefs[ind] = el % bch->field->characteristic;
            el >>= 1;
            ++ind;
            if(ind >= bch->data_length + bch->generator->degree - 1){
                ind = 0;
                decoded = decode_bch(bch, decode_poly);
                for(int j = 0; j < decoded->degree; ++j){
                    not_decoded = 0;
                    result[*decoded_length] += (1ULL << last_power) * decoded->coefs[j];
                    ++last_power;
                    if(last_power >= (uint64_t)threshold){
                        ++(*decoded_length);
                        last_power = 0;
                    }
                    if(*decoded_length >= size){
                        temp = calloc(size * 2, sizeof(uint64_t));
                        if(temp == NULL){
                            printf("NOT ALLOC'D\n");
                            return NULL;
                        }
                        memcpy(temp, result, size * sizeof(uint64_t));
                        size *= 2;
                        free(result);
                        result = temp;
                    }
                }
                free_extended_polynomial(decoded);
            }
        }
    }
    free_extended_polynomial(decode_poly);
    return result;
}
*/