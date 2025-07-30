#include <stdbool.h>

#include "polynomials.h"

#ifndef FIN_FIELD  
#define FIN_FIELD

struct finite_field{
    unsigned int characteristic;
    uint64_t power;
    uint64_t primative_in_power_n;
};

uint64_t
add_in_field(struct finite_field* field, uint64_t a, uint64_t b);

uint64_t
multiply_in_field(struct finite_field* field, uint64_t a, uint64_t b);

uint64_t
find_primitive_in_power(struct finite_field* field, int power);

uint64_t
construct_inverse_element_add(struct finite_field* field, uint64_t element);

uint64_t
get_primitive_polynomial(struct finite_field* field);

struct finite_field*
construct_gf_p_m(uint64_t characteristic, uint64_t power);

int
check_if_irreducible(uint64_t possible_polynomial, uint64_t p, uint64_t m);

void
get_polynomial_for_irreducible(uint64_t p, uint64_t d, uint64_t *need_check_gcd);

int
check_if_primitive(uint64_t possible_polynomial, uint64_t p, uint64_t m);

uint64_t*
find_prime_divisors(uint64_t n, int *length);

#endif