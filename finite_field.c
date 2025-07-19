#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "finite_field.h"


unsigned long long
add_in_field(struct finite_field* field, unsigned long long a, unsigned long long b)
{
    if(field->characteristic == 2) return a ^ b;
    unsigned long long int a_copy = a, b_copy = b, result = 0;
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

unsigned long long
multiply_in_field(struct finite_field* field, unsigned long long a, unsigned long long b)
{   
    if(field->characteristic == 2){
        unsigned long long result = 0;
        unsigned long long max_val = 1 << field->power;
        unsigned long long a_copy = a, b_copy = b;
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
    unsigned long long result = 0;
    unsigned long long max_val = pow(field->characteristic, field->power);
    unsigned long long a_copy = a, b_copy = b;
    while(a_copy != 0){
        if(a_copy % field->characteristic != 0){
            for(int i = 0; (unsigned long long)i < a_copy % field->characteristic; ++i){
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

unsigned long long
find_primitive_in_power(struct finite_field* field, int power)
{
    if(power < 0){
        unsigned long long int temp = find_primitive_in_power(field, -power);
        return construct_inverse_element_multiply(field, temp);
    }
    if((unsigned long long)power < field->power){
        if(field->characteristic == 2){
            return (1 << power);
        }
        return pow(field->characteristic, power);
    }else if(power % 2 == 0){
        unsigned long long a = find_primitive_in_power(field, power/ 2);
        return  multiply_in_field(field, a, a);
    }
    return multiply_in_field(field, find_primitive_in_power(field, power - 1), field->characteristic);
}

unsigned long long
construct_inverse_element_add(struct finite_field* field, unsigned long long element)
{
    if(field->characteristic == 2) return element;
    unsigned long long element_copy = element;
    unsigned long long result = 0;
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

unsigned long long
get_primitive_polynomial(struct finite_field* field)
{
    unsigned long long result = pow(field->characteristic, field->power);
    result += construct_inverse_element_add(field, field->primative_in_power_n);
    return result;
}

struct finite_field*
construct_gf_p_m(unsigned long long characteristic, unsigned long long power)
{
    unsigned long long possible_polynomial = (unsigned long long) pow((double)characteristic, (double)power);
    for(unsigned long long i = 1; i < possible_polynomial; ++i){
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
check_if_irreducible(unsigned long long possible_polynomial, unsigned long long p, unsigned long long m)
{
    struct finite_field check_field;
    check_field.characteristic = p;
    check_field.power = m;
    check_field.primative_in_power_n = 100;
    struct extended_polynomial *s0, *t0, *r0, *possible_primitive;
    bool result = 1;
    for(unsigned long long i = 1; i < m; ++i){
        possible_primitive = construct_polynomial_from_field_element(&check_field, possible_polynomial);
        extended_euclidean_algorithm(&check_field, possible_primitive, 
        get_polynomial_for_irruducuble(p, i), &s0, &t0, &r0
        );
        
        free_extended_polynomial(s0);
        free_extended_polynomial(t0);
        if(r0->degree != 1 || r0->coefs[0] != 1){
            result = 0;
            free_extended_polynomial(r0);
            break;
        }
        free_extended_polynomial(r0);
        free_extended_polynomial(possible_primitive);
        
    }
    return result;
}


struct extended_polynomial*
get_polynomial_for_irruducuble(unsigned long long p, unsigned long long d)
{
    struct extended_polynomial* result = make_zero_polynomial(1);
    result->degree = 1;
    result->coefs[0] = 1;
    result = multiply_extended_polynomial_by_x_n(result, pow(p, d), NEED_FREE);
    result->coefs[1] += p - 1;
    return result;
}

int
check_if_primitive(unsigned long long possible_polynomial, unsigned long long p, unsigned long long m)
{
    unsigned long long int n = (unsigned long long)pow((double)p, (double)m);

    int result = 1;

    struct finite_field test_field;
    test_field.characteristic = p;
    test_field.power = m;
    test_field.primative_in_power_n = construct_inverse_element_add(&test_field, possible_polynomial - n);
    int len_primes = 0;
    unsigned long long *primes = find_prime_divisors(n - 1, &len_primes);

    for(int i = 0; i < len_primes; ++i){
        if(find_primitive_in_power(&test_field, n/primes[i]) == 1){
            result = 0;
            break;
        }
    }
    free(primes);
    return result;
}

unsigned long long*
find_prime_divisors(unsigned long long n, int *length)
{
    unsigned long long *primes = malloc(sizeof(unsigned long long) * 64), *new_primes;
    int size = 64, count = 0;

    if (n % 2 == 0) {
        primes[count] = 2;
        ++count;
        if(count >= size){
            new_primes = malloc(sizeof(unsigned long long) * size * 2);
            memcpy(new_primes, primes, size);
            size *= 2;
            free(primes);
            primes = new_primes;
        }
        while (n % 2 == 0) {
            n >>= 1;
        }
    }
    
    unsigned long long limit = (unsigned long long) sqrt(n);
    for (unsigned long long d = 3; d <= limit; d += 2) {
        if (n % d == 0) {
            primes[count] = d;
            if(count >= size){
                new_primes = malloc(sizeof(unsigned long long) * size * 2);
                memcpy(new_primes, primes, size);
                size *= 2;
                free(primes);
                primes = new_primes;
            }
            ++count;
            while (n % d == 0) {
                n /= d;
            }
            limit = (unsigned long long) sqrt(n);
        }
    }

    if (n > 1) {
        primes[count] = n;
        ++count;

    }
    *length = count;
    return primes;
}