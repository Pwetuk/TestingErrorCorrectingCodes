#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "finite_field.h"


uint64_t
add_in_field(struct finite_field* field, uint64_t a, uint64_t b)
{
    if(field->characteristic == 2) return a ^ b;
    uint64_t a_copy = a, b_copy = b, result = 0;
    int i = 0;
    while (a_copy != 0)
    {
        result += pow((double)field->characteristic, (double)i) * ((a_copy + b_copy) % field->characteristic);
        a_copy /= (field->characteristic);
        b_copy /= (field->characteristic);
        ++i;
    }
    result += pow(field->characteristic, i) * ((a_copy + b_copy));
    return result;
}

uint64_t
multiply_in_field(struct finite_field* field, uint64_t a, uint64_t b)
{   
    if(field->characteristic == 2){
        uint64_t result = 0;
        uint64_t max_val = 1 << field->power;
        uint64_t a_copy = a, b_copy = b;
        while (a_copy) {
            if(a_copy % 2 == 1){
                result ^= b_copy;
            }
            a_copy >>= 1;
            b_copy <<= 1;
            if(b_copy >= max_val){
                b_copy ^= max_val;
                b_copy ^= field->primative_in_power_n;
            }
        }
        return result;
    }
    uint64_t result = 0;
    uint64_t max_val = pow(field->characteristic, field->power);
    uint64_t a_copy = a, b_copy = b;
    while(a_copy != 0){
        if(a_copy % field->characteristic != 0){
            for(int i = 0; (uint64_t)i < a_copy % field->characteristic; ++i){
                result = add_in_field(field, result, b_copy);
            }
        }
        a_copy /= field->characteristic;
        b_copy *= field->characteristic;
        while(b_copy >= max_val){
            b_copy -= max_val;
            b_copy = add_in_field(field, b_copy, field->primative_in_power_n);
        }
    }
    return result;
}

uint64_t
find_primitive_in_power(struct finite_field* field, int power)
{
    if(power < 0){
        uint64_t temp = find_primitive_in_power(field, -power);
        return construct_inverse_element_multiply(field, temp);
    }
    if((uint64_t)power < field->power){
        if(field->characteristic == 2){
            return (1 << power);
        }
        return pow(field->characteristic, power);
    }else if(power % 2 == 0){
        uint64_t a = find_primitive_in_power(field, power/ 2);
        return  multiply_in_field(field, a, a);
    }
    return multiply_in_field(field, find_primitive_in_power(field, power - 1), field->characteristic);
}

uint64_t
construct_inverse_element_add(struct finite_field* field, uint64_t element)
{
    if(field->characteristic == 2) return element;
    uint64_t element_copy = element;
    uint64_t result = 0;
    int i = 0;
    while(element_copy != 0){
        if(element_copy % field->characteristic != 0){
            result += (field->characteristic - (element_copy % field->characteristic)) * pow(field->characteristic, i);
        }
        element_copy /= field->characteristic;
        ++i;
    }
    return result;
}

uint64_t
get_primitive_polynomial(struct finite_field* field)
{
    uint64_t result = pow(field->characteristic, field->power);
    result += construct_inverse_element_add(field, field->primative_in_power_n);
    return result;
}

struct finite_field*
construct_gf_p_m(uint64_t characteristic, uint64_t power)
{
    uint64_t possible_polynomial = (uint64_t) pow((double)characteristic, (double)power);
    for(uint64_t i = 1; i < possible_polynomial; ++i){
        if(check_if_irreducible(possible_polynomial + i, characteristic, power) == 1){
            if(check_if_primitive(possible_polynomial + i, characteristic, power) == 1){
                struct finite_field* result = malloc(sizeof(struct finite_field));
                result->characteristic = characteristic;
                result->power = power;
                result->primative_in_power_n = construct_inverse_element_add(result, i);
                return result;
            }
        }
    }
    return NULL;
}

int
check_if_irreducible(uint64_t possible_polynomial, uint64_t p, uint64_t m)
{
    struct finite_field check_field;
    check_field.characteristic = p;
    check_field.power = m;
    check_field.primative_in_power_n = 100;
    uint64_t s0[MAX_DEGREE], t0[MAX_DEGREE], r0[MAX_DEGREE], possible_primitive[MAX_DEGREE], need_check_gcd[MAX_DEGREE];
    bool result = 1;
    construct_polynomial_from_field_element(&check_field, possible_polynomial, possible_primitive);
    for(uint64_t i = 1; i < m; ++i){
        get_polynomial_for_irruducuble(p, i, need_check_gcd);
        extended_euclidean_algorithm(&check_field, possible_primitive, 
        need_check_gcd, s0, t0, r0
        );
        
        
        if(get_degree(r0) != 1 || r0[0] != 1){
            result = 0;
            break;
        }
        
    }
    return result;
}


void
get_polynomial_for_irruducuble(uint64_t p, uint64_t d, uint64_t need_check_gcd[MAX_DEGREE])
{
    memset(need_check_gcd, 0, MAX_DEGREE * sizeof(uint64_t));
    need_check_gcd[1] = 1;
    
    need_check_gcd[1 << d] += 1;
}

int
check_if_primitive(uint64_t possible_polynomial, uint64_t p, uint64_t m)
{
    uint64_t n = (uint64_t)pow((double)p, (double)m);

    int result = 1;

    struct finite_field test_field;
    test_field.characteristic = p;
    test_field.power = m;
    test_field.primative_in_power_n = construct_inverse_element_add(&test_field, possible_polynomial - n);
    int len_primes = 0;
    uint64_t *primes = find_prime_divisors(n - 1, &len_primes);

    for(int i = 0; i < len_primes; ++i){
        if(find_primitive_in_power(&test_field, n/primes[i]) == 1){
            result = 0;
            break;
        }
    }
    free(primes);
    return result;
}

uint64_t*
find_prime_divisors(uint64_t n, int *length)
{
    uint64_t *primes = malloc(sizeof(uint64_t) * 64), *new_primes;
    int size = 64, count = 0;

    if (n % 2 == 0) {
        primes[count] = 2;
        ++count;
        if(count >= size){
            new_primes = malloc(sizeof(uint64_t) * size * 2);
            memcpy(new_primes, primes, size);
            size *= 2;
            free(primes);
            primes = new_primes;
        }
        while (n % 2 == 0) {
            n >>= 1;
        }
    }
    
    uint64_t limit = (uint64_t) sqrt(n);
    for (uint64_t d = 3; d <= limit; d += 2) {
        if (n % d == 0) {
            primes[count] = d;
            if(count >= size){
                new_primes = malloc(sizeof(uint64_t) * size * 2);
                memcpy(new_primes, primes, size);
                size *= 2;
                free(primes);
                primes = new_primes;
            }
            ++count;
            while (n % d == 0) {
                n /= d;
            }
            limit = (uint64_t) sqrt(n);
        }
    }

    if (n > 1) {
        primes[count] = n;
        ++count;

    }
    *length = count;
    return primes;
}