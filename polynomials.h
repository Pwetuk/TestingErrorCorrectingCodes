#include "finite_field.h"

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


polynomial
add_two_polynomials(struct finite_field*, polynomial, polynomial);

polynomial
multiply_two_polynomials(struct finite_field*, polynomial, polynomial);


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
multiplty_two_extended_polynomials(struct finite_field*, struct extended_polynomial*, struct extended_polynomial*);


struct extended_polynomial*
copy_extended_polynomial(struct extended_polynomial* poly);


polynomial
construct_minimal_polynomial_from_coset(struct finite_field* field, struct cyclotomic_coset* coset);

polynomial
extended_polynomial_to_polynomial(struct finite_field* field, struct extended_polynomial* poly);

polynomial
construct_generator_polynomial(struct finite_field* field, int number_of_errors);
