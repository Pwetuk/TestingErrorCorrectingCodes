#include <stdio.h>
#include <string.h>

#include "bch_tests.h"


int main(int argc, char** argv){
    int tests_result = 0;
    int was_testing = 1;
    if(argc > 1 && strcmp(argv[1], "test") == 0){
        tests_result = run_tests();
    }else if(argc > 1 && strcmp(argv[1], "ttest") == 0){
        tests_result = temp_tests();
    }else{
        was_testing = 0;
    }

    if(was_testing == 1){
        if(tests_result == 1){
            printf("All tests succeded\n");
        }else{
            printf("Something went wrong\n");
        }
        return 0;
    }

    unsigned int characteristic;
    unsigned long long power;
    unsigned long long to_encode;
    int number_of_errors;
    struct extended_polynomial* to_encode_poly;

    printf("Input characteristic, power, number of errors polynomial\n");
    scanf("%u %llu %d", &characteristic, &power, &number_of_errors);

    struct bch_code* bch_test = init_bch(characteristic, power, number_of_errors, 0);


    printf("%lluMax value of encoded message is: %llu\n", bch_test->field->primative_in_power_n, 1ULL << (bch_test->data_length - 1));
    scanf("%llu", &to_encode);

    to_encode_poly = construct_polynomial_from_field_element(bch_test->field, to_encode);

    bch_test->data_length = to_encode_poly->degree;

    struct extended_polynomial* result = encode_bch(bch_test, to_encode_poly);
    

    unsigned long long encoded = extended_polynomial_to_polynomial(bch_test->field, result);
    printf("%lluEncoded message: %llu\n", bch_test->field->primative_in_power_n, encoded);

    unsigned long long to_decode;

    printf("Input message to decode\n");
    scanf("%llu", &to_decode);
    
    struct extended_polynomial* to_decode_poly;

    to_decode_poly = construct_polynomial_from_field_element(bch_test->field, to_decode);


    unsigned long long decoded =  decode_bch(bch_test, to_decode_poly);
    
    printf("Decoded message: %lld\n", decoded);

    

    free_extended_polynomial(result);
    free_bch_code(bch_test);

    return 0;
}