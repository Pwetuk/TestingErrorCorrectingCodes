#include "polynomials.h"

struct bch_code{
    struct finite_field* field;
    int number_of_errors;
    int data_length;
    struct extended_polynomial* generator;
};

struct extended_polynomial*
encode_bch(struct bch_code* bch_code_struct, struct extended_polynomial* data);

unsigned long long*
calculate_syndrome(struct finite_field* field, struct extended_polynomial* encoded_message, int number_of_errors);

struct extended_polynomial*
construct_locator_polynomial(struct finite_field* field, unsigned long long* syndrome, int number_of_errors);

int*
chien_search(struct finite_field* field, struct extended_polynomial* locator, int number_of_errors);

void
fix_errors_in_bch_binary(struct extended_polynomial* codeword, int* errors_pos, int number_of_errors);

unsigned long long
decode_bch(struct bch_code* bch_code_struct, struct extended_polynomial* message);

struct bch_code*
init_bch(unsigned long long p, unsigned long long power, int number_of_errors, unsigned long long primitive);

void
free_bch_code(struct bch_code* bch_code_struct);