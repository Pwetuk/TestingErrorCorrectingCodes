#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bch_encode.h"


struct extended_polynomial*
encode_bch(struct bch_code* bch_code_struct, struct extended_polynomial* data)
{
    uint64_t n = pow(bch_code_struct->field->characteristic, bch_code_struct->field->power);
    struct extended_polynomial* result, *rresult = make_zero_polynomial(n);
    data = multiply_extended_polynomial_by_x_n(data, n - bch_code_struct->data_length - 1, 0);
    struct extended_polynomial *quotient = NULL, *remainder = NULL;
    divide_polynomials_with_remainder(bch_code_struct->field, data, bch_code_struct->generator, &quotient, &remainder);
    result = difference_of_two_polynomials(bch_code_struct->field, data, remainder);
    if(quotient == NULL){
        printf("Something went horribly wrong\n");
    }
    free_extended_polynomial(quotient);
    free_extended_polynomial(remainder);

    free_extended_polynomial(data);
    rresult->degree = n - 1;
    for(int i = 0; i < result->degree; ++i){
        rresult->coefs[i] = result->coefs[i];
    }
    free_extended_polynomial(result);
    return rresult;
}

struct extended_polynomial*
decode_bch(struct bch_code* bch_code_struct, struct extended_polynomial* message)
{
    struct extended_polynomial* locator;
    int* errors_pos;
    uint64_t* syndrome = calculate_syndrome(bch_code_struct->field, message, bch_code_struct->number_of_errors);


    uint8_t is_null = 0;
    for(int i = 0; i < bch_code_struct->number_of_errors; ++i){
        if(syndrome[i] != 0){
            is_null = 1;
            break;
        }
    }

    if(is_null == 1){
        locator = construct_locator_polynomial(bch_code_struct->field, syndrome, bch_code_struct->number_of_errors);
    
        errors_pos = chien_search(bch_code_struct->field, locator, bch_code_struct->number_of_errors);
    
        fix_errors_in_bch_binary(message, errors_pos, bch_code_struct->number_of_errors);
    }

    struct extended_polynomial* result = make_zero_polynomial(bch_code_struct->data_length);
    result->degree = bch_code_struct->data_length;

    for(int i = 0; i < bch_code_struct->data_length; ++i){
        result->coefs[i] = message->coefs[message->degree - result->degree + i];
    }

    free(syndrome);
    if(is_null == 1){
        free(errors_pos);
        free_extended_polynomial(locator);
    }
    return result;
}

uint64_t*
calculate_syndrome(struct finite_field* field, struct extended_polynomial* encoded_message, int number_of_errors)
{
    uint64_t *syndrome = malloc(sizeof(uint64_t) * 2 * number_of_errors);
    for(int i = 0; i < 2 * number_of_errors; ++i){
        syndrome[i] = find_value_from_root(field, encoded_message, i + 1);
    }
    return syndrome;
}



struct extended_polynomial*
construct_locator_polynomial(struct finite_field* field, uint64_t* syndrome, int number_of_errors){
    needed_type base_coef[2] = {1, 0};
    struct extended_polynomial* locator = construct_extended_polynomial_from_coefs(base_coef, 1);
    struct extended_polynomial* B = construct_extended_polynomial_from_coefs(base_coef, 1);
    struct extended_polynomial* locator_old, *B_copy;
    int L = 0, m = 1;
    uint64_t delta_r, b = 1;
    for(int r = 0; r < 2 * number_of_errors; ++r){
        
        delta_r = syndrome[r];
        for(int i = 1; i <= L && i < locator->degree; ++i){
            delta_r = add_in_field(field, delta_r, multiply_in_field(field, syndrome[r - i], locator->coefs[i]));
        }        

        if(delta_r == 0){
            ++m;
        }else{
            locator_old = locator;
            B_copy = multiply_extended_polynomial_by_x_n(B, m, 0);
            uint64_t coef = multiply_in_field(field, delta_r, construct_inverse_element_multiply(field, b));
            B_copy = multiply_polynomial_by_element(field, B_copy, coef, NEED_FREE);
            locator = difference_of_two_polynomials(field, locator, B_copy);
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
    }

    free_extended_polynomial(B);
    return locator;
}

int*
chien_search(struct finite_field* field, struct extended_polynomial* locator, int number_of_errors)
{
    uint64_t n = pow(field->characteristic, field->power);
    uint64_t pos;
    int err_pos = 0;
    int *result = malloc(number_of_errors * sizeof(int));
    for(int i = 0; i < number_of_errors; ++i){
        result[i] = -1;
    }
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
    to_return->generator = construct_generator_polynomial(to_return->field, number_of_errors);
    to_return->data_length = pow(to_return->field->characteristic, to_return->field->power) - to_return->generator->degree;
    return to_return;
}

void
free_bch_code(struct bch_code* bch_code_struct)
{
    free(bch_code_struct->field);
    free_extended_polynomial(bch_code_struct->generator);
    free(bch_code_struct);
}

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
