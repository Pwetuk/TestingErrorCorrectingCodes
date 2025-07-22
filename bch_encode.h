#include "polynomials.h"

struct bch_code{
    struct finite_field* field;
    int number_of_errors;
    int data_length;
    struct extended_polynomial* generator;
};

struct extended_polynomial*
encode_bch(struct bch_code* bch_code_struct, struct extended_polynomial* data);

uint64_t*
calculate_syndrome(struct finite_field* field, struct extended_polynomial* encoded_message, int number_of_errors);

struct extended_polynomial*
construct_locator_polynomial(struct finite_field* field, uint64_t* syndrome, int number_of_errors);

int*
chien_search(struct finite_field* field, struct extended_polynomial* locator, int number_of_errors);

void
fix_errors_in_bch_binary(struct extended_polynomial* codeword, int* errors_pos, int number_of_errors);

struct extended_polynomial*
decode_bch(struct bch_code* bch_code_struct, struct extended_polynomial* message);

struct bch_code*
init_bch(uint64_t p, uint64_t power, int number_of_errors, uint64_t primitive);

void
free_bch_code(struct bch_code* bch_code_struct);

uint64_t*
bch_encode_arr(struct bch_code* bch, uint64_t* to_encode_arr, uint64_t to_encode_length, uint64_t *encoded_length);

uint64_t*
bch_decode_arr(struct bch_code* bch, uint64_t* to_decode_arr, uint64_t to_decode_length, uint64_t *decoded_length);