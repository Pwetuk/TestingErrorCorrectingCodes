#include <stdbool.h>

#include "polynomials.h"

#ifndef FIN_FIELD  
#define FIN_FIELD

struct finite_field{
    unsigned int characteristic;
    unsigned long long int power;
    unsigned long long int primative_in_power_n;
};

unsigned long long int
add_in_field(struct finite_field* field, unsigned long long int a, unsigned long long int b);

unsigned long long int
multiply_in_field(struct finite_field* field, unsigned long long int a, unsigned long long int b);

unsigned long long int
find_primitive_in_power(struct finite_field* field, int power);

unsigned long long
construct_inverse_element_add(struct finite_field* field, unsigned long long element);

unsigned long long
get_primitive_polynomial(struct finite_field* field);

struct finite_field*
construct_gf_p_m(unsigned long long characteristic, unsigned long long power);

int
check_if_irreducible(unsigned long long possible_polynomial, unsigned long long p, unsigned long long m);

struct extended_polynomial*
get_polynomial_for_irruducuble(unsigned long long p, unsigned long long d);

int
check_if_primitive(unsigned long long possible_polynomial, unsigned long long p, unsigned long long m);

unsigned long long*
find_prime_divisors(unsigned long long n, int *length);

#endif