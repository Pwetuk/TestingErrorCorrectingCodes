#include <stdint.h>

#ifndef POLYNOMIALS  
#define POLYNOMIALS

#include "finite_field.h"

#define MAX_DEGREE 32

typedef uint64_t needed_type;

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

struct cyclotomic_coset
{
    int size;
    int from_el;
};

struct cyclotomic_cosets{
    int number_of_cosets;
    struct cyclotomic_coset** cosets;
};


typedef uint64_t polynomial;

struct extended_polynomial{
    int degree;
    uint64_t *coefs;
};

struct cyclotomic_cosets*
generate_cosets_from_n_needed_elements(struct finite_field* field, int n);


void
free_cosets_struct(struct cyclotomic_cosets*);

void
print_polynomial(struct finite_field* field, polynomial poly);

void
free_extended_polynomial(struct extended_polynomial*);

void
construct_extended_polynomial_from_coefs(uint64_t* coefs, int n, uint64_t res[MAX_DEGREE]);

void
print_extended_polynomial(uint64_t* to_print);

void
add_two_extended_polynomials(struct finite_field*, uint64_t*, uint64_t*);

void
multiply_two_extended_polynomials(struct finite_field*, uint64_t*, uint64_t*);


void
copy_extended_polynomial(uint64_t* poly, uint64_t* res);


void
construct_minimal_polynomial_from_coset(struct finite_field* field, struct cyclotomic_coset* coset, uint64_t res[MAX_DEGREE]);

polynomial
extended_polynomial_to_polynomial(struct finite_field* field, uint64_t* poly);

void
construct_generator_polynomial(struct finite_field* field, int number_of_errors, uint64_t res[MAX_DEGREE]);

void
find_inverse_by_add_for_polynomial(struct finite_field* field, uint64_t a[MAX_DEGREE], uint64_t res[MAX_DEGREE]);

void
difference_of_two_polynomials(struct finite_field* field, uint64_t* a, uint64_t* b, uint64_t res[MAX_DEGREE]);

void
divide_polynomials_with_remainder(struct finite_field* field, uint64_t* dividing, uint64_t* divider,
    uint64_t result[MAX_DEGREE], uint64_t remainder[MAX_DEGREE]);

uint64_t
find_value_from_root(struct finite_field* field, uint64_t poly[MAX_DEGREE], int primitive_power);

int
equal_polynomials(uint64_t* a, uint64_t* b);

struct extended_polynomial*
make_zero_polynomial(int degree);


uint64_t
construct_inverse_element_multiply(struct finite_field* field, uint64_t el);

void
construct_polynomial_from_field_element(struct finite_field* field, uint64_t el, uint64_t[MAX_DEGREE]);

uint64_t
get_field_el_from_polynomial(struct finite_field* field, uint64_t poly[MAX_DEGREE]);

void
multiply_polynomial_by_element(struct finite_field* field, uint64_t poly[MAX_DEGREE], uint64_t el);

void
extended_euclidean_algorithm(struct finite_field* field, uint64_t* el, uint64_t* f, 
    uint64_t s0[MAX_DEGREE], uint64_t t0[MAX_DEGREE], uint64_t r0[MAX_DEGREE]);

int
get_degree(uint64_t poly[MAX_DEGREE]);

#endif