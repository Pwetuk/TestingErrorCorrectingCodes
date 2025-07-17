#include "polynomials.h"

struct extended_polynomial*
encode_bch(struct finite_field* field, struct extended_polynomial* generator, struct extended_polynomial* data, int n, int degree_data);

unsigned long long*
calculate_syndrome(struct finite_field* field, struct extended_polynomial* encoded_message, int number_of_errors);

struct extended_polynomial*
construct_locator_polynomial(struct finite_field* field, unsigned long long* syndrome, int number_of_errors);

int*
chien_search(struct finite_field* field, struct extended_polynomial* locator, int number_of_errors);

void
fix_errors_in_bch_binary(struct extended_polynomial* codeword, int* errors_pos, int number_of_errors);

unsigned long long
decode_bch(struct finite_field* field, struct extended_polynomial* message, int number_of_errors, int original_data_length);