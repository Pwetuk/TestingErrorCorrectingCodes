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


typedef unsigned long long int polynoms;

struct extended_polynom{
    int degree;
    unsigned long long int* coefs;
};

struct cyclotomic_cosets*
generate_cosets_from_n_needed_elements(struct binary_field* field, int n);


void
free_cosets_struct(struct cyclotomic_cosets*);


polynoms
add_two_polynoms(struct binary_field*, polynoms, polynoms);

polynoms
multiply_two_polynoms(struct binary_field*, polynoms, polynoms);


void
print_polynom(struct binary_field* field, polynoms poly);

void
free_extended_polynom(struct extended_polynom*);

struct extended_polynom*
construct_extended_polynom_from_coefs(int* coefs, int n);

void
print_extended_polynom(struct extended_polynom* to_print);

struct extended_polynom*
add_two_extended_polynoms(struct binary_field*, struct extended_polynom*, struct extended_polynom*);

struct extended_polynom*
multiplty_two_extended_polynoms(struct binary_field*, struct extended_polynom*, struct extended_polynom*);


struct extended_polynom*
copy_extended_polynom(struct extended_polynom* poly);


polynoms
construct_minimal_polynom_from_coset(struct binary_field* field, struct cyclotomic_coset* coset);

polynoms
extended_polynom_to_polynoms(struct binary_field* field, struct extended_polynom* poly);

polynoms
construct_generator_polynomial(struct binary_field* field, int number_of_errors);
