#include "polynomials.h"

struct bch_code{
    struct finite_field* field;
    int number_of_errors;
    int data_length;
    uint64_t generator[MAX_DEGREE];
    uint64_t n;
};

void
encode_bch(struct bch_code* bch_code_struct, uint64_t data[MAX_DEGREE], uint64_t result[MAX_DEGREE]);

int
calculate_syndrome(struct finite_field* field, uint64_t encoded_message[MAX_DEGREE], int number_of_errors, uint64_t syndrome[MAX_DEGREE]);

void
construct_locator_polynomial(struct finite_field* field, uint64_t* syndrome, int number_of_errors, uint64_t locator[MAX_DEGREE]);

void
chien_search(struct finite_field* field, uint64_t* locator, int number_of_errors, int* errors, uint64_t n);

void
fix_errors_in_bch_binary(uint64_t* codeword, int* errors_pos, int number_of_errors);

void
decode_bch(struct bch_code* bch_code_struct, uint64_t message[MAX_DEGREE], uint64_t result[MAX_DEGREE]);

struct bch_code*
init_bch(uint64_t p, uint64_t power, int number_of_errors, uint64_t primitive);

void
free_bch_code(struct bch_code* bch_code_struct);

uint64_t*
bch_encode_arr(struct bch_code* bch, uint64_t* to_encode_arr, uint64_t to_encode_length, uint64_t *encoded_length);

uint64_t*
bch_decode_arr(struct bch_code* bch, uint64_t* to_decode_arr, uint64_t to_decode_length, uint64_t *decoded_length);