#ifndef POLYNOMIALS  
#define POLYNOMIALS

#include "finite_field.h"

#define NEED_FREE 1

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


typedef unsigned long long int polynomial;

struct extended_polynomial{
    int degree;
    unsigned long long int* coefs;
};

struct cyclotomic_cosets*
generate_cosets_from_n_needed_elements(struct finite_field* field, int n);


void
free_cosets_struct(struct cyclotomic_cosets*);

void
print_polynomial(struct finite_field* field, polynomial poly);

void
free_extended_polynomial(struct extended_polynomial*);

struct extended_polynomial*
construct_extended_polynomial_from_coefs(int* coefs, int n);

void
print_extended_polynomial(struct extended_polynomial* to_print);

struct extended_polynomial*
add_two_extended_polynomials(struct finite_field*, struct extended_polynomial*, struct extended_polynomial*);

struct extended_polynomial*
multiply_two_extended_polynomials(struct finite_field*, struct extended_polynomial*, struct extended_polynomial*);


struct extended_polynomial*
copy_extended_polynomial(struct extended_polynomial* poly);


struct extended_polynomial*
construct_minimal_polynomial_from_coset(struct finite_field* field, struct cyclotomic_coset* coset);

polynomial
extended_polynomial_to_polynomial(struct finite_field* field, struct extended_polynomial* poly);

struct extended_polynomial*
construct_generator_polynomial(struct finite_field* field, int number_of_errors);


struct extended_polynomial*
find_inverse_by_add_for_polynomial(struct finite_field* field, struct extended_polynomial* a);

struct extended_polynomial*
difference_of_two_polynomials(struct finite_field* field, struct extended_polynomial* a, struct extended_polynomial* b);

void
divide_polynomials_with_remainder(struct finite_field* field, struct extended_polynomial* dividing, struct extended_polynomial* divider,
    struct extended_polynomial** result, struct extended_polynomial** remainder);

unsigned long long
find_value_from_root(struct finite_field* field, struct extended_polynomial* poly, int primitive_power);

struct extended_polynomial*
multiply_extended_polynomial_by_x_n(struct extended_polynomial* poly, int n, int need_free);


unsigned long long
construct_inverse_element_multiply(struct finite_field* field, unsigned long long element);

int
equal_polynomials(struct extended_polynomial* a, struct extended_polynomial* b);

struct extended_polynomial*
make_zero_polynomial(int degree);


unsigned long long
construct_inverse_element_multiply(struct finite_field* field, unsigned long long el);

struct extended_polynomial*
construct_polynomial_from_field_element(struct finite_field* field, unsigned long long el);

unsigned long long
get_field_el_from_polynomial(struct finite_field* field, struct extended_polynomial* poly, int need_free);

struct extended_polynomial*
multiply_polynomial_by_element(struct finite_field* field, struct extended_polynomial* poly, unsigned long long int el, int need_free);

#endif