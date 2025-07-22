#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "bch_tests.h"

int 
run_tests()
{
    return test_in_characteristics_2_first() &
        test_in_characteristics_2_second() &
        test_in_characteristics_2_third() &
        //test_in_characteristics_n_first() &
        //test_in_characteristics_n_second() &
        //test_in_characteristics_n_third() &
        syndrom_test() &
        //encode_test1() &
        inverse_test() &
        locator_test() &
        chien_test() &
        decode_test() &
        //encode_test2() &
        1;
}

int
temp_tests()
{
    //if(temp_test1() != 1){
    //    printf("Test1\n");
    //    return 0;
    //}
    //if(temp_test2() != 1){
    //    printf("Test2\n");
    //    return 0;
    //}
    //if(temp_test3() != 1){
    //    printf("Test3\n");
    //    return 0;
    //}
    /*
    if(temp_test4() != 1){
        printf("Test4\n");
        return 0;
    }
    if(temp_test5() != 1){
        printf("Test1\n");
        return 0;
    }
    */
    return 1;
}


int
test_in_characteristics_2_first()
{
    struct finite_field test_field;
    test_field.characteristic = 2;
    test_field.power = 5;
    test_field.primative_in_power_n = 5;

    uint64_t true_answer[MAX_DEGREE] = {1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1};
    memset(true_answer + 11, 0, (MAX_DEGREE - 11) * sizeof(uint64_t));
    uint64_t test_result[MAX_DEGREE];
    construct_generator_polynomial(&test_field, 2, test_result);
    int result = equal_polynomials(true_answer, test_result);


    
    if(result != 1){
        printf("Test1 failed\n");
        print_extended_polynomial(test_result);
        print_extended_polynomial(true_answer);
        printf("----------------------------\n");
    }
    

    return result;
}


int
test_in_characteristics_2_second()
{
    struct finite_field test_field;
    test_field.characteristic = 2;
    test_field.power = 6;
    test_field.primative_in_power_n = 3;


    uint64_t true_answer[MAX_DEGREE] = {1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1};;
    
    memset(true_answer + 13, 0, (MAX_DEGREE - 13) * sizeof(uint64_t));


    uint64_t test_result[MAX_DEGREE];
    construct_generator_polynomial(&test_field, 2, test_result);
    int result = equal_polynomials(true_answer, test_result);
    
    if(result != 1){
        printf("Test2 failed\n");
        print_extended_polynomial(test_result);
        print_extended_polynomial(true_answer);
        printf("----------------------------\n");
    }
    

    return result;
}

int
test_in_characteristics_2_third()
{
    struct finite_field test_field;
    test_field.characteristic = 2;
    test_field.power = 4;
    test_field.primative_in_power_n = 3;

    uint64_t true_answer[MAX_DEGREE] = {1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1};

    memset(true_answer + 11, 0, (MAX_DEGREE - 11) * sizeof(uint64_t));

    uint64_t test_result[MAX_DEGREE];
    construct_generator_polynomial(&test_field, 3, test_result);
    int result = equal_polynomials(true_answer, test_result);
    
    if(result != 1){
        printf("Test3 failed\n");
        print_extended_polynomial(test_result);
        print_extended_polynomial(true_answer);
        printf("----------------------------\n");
    }
    

    return result;
}

/*
int
test_in_characteristics_n_first()
{
    struct finite_field test_field;
    test_field.characteristic = 3;
    test_field.power = 5;
    test_field.primative_in_power_n = 25;

    struct extended_polynomial true_answer;
    needed_type coefficients[16] = {2, 2, 1, 2, 2, 1, 0, 2, 1, 1, 1, 1, 1, 0, 1, 1};
    true_answer.degree = 16;
    true_answer.coefs = coefficients;



    struct extended_polynomial* test_result = construct_generator_polynomial(&test_field, 2);
    int result = equal_polynomials(&true_answer, test_result);
    
    if(result != 1){
        printf("Test n 1 failed\n");
        print_extended_polynomial(test_result);
        print_extended_polynomial(&true_answer);
        printf("----------------------------\n");
    }
    
    free_extended_polynomial(test_result);

    return result;
}


int
test_in_characteristics_n_second()
{
    struct finite_field test_field;
    test_field.characteristic = 5;
    test_field.power = 3;
    test_field.primative_in_power_n = 4 + 15;

    struct extended_polynomial true_answer;
    needed_type coefficients[13] = {1, 0, 1, 0, 1, 3, 2, 4, 2, 1, 2, 4, 1};
    true_answer.degree = 13;
    true_answer.coefs = coefficients;



    struct extended_polynomial* test_result = construct_generator_polynomial(&test_field, 2);
    int result = equal_polynomials(&true_answer, test_result);
    
    if(result != 1){
        printf("Test n 2 failed\n");
        print_extended_polynomial(test_result);
        print_extended_polynomial(&true_answer);
        printf("----------------------------\n");
    }
    
    free_extended_polynomial(test_result);

    return result;

}


int
test_in_characteristics_n_third()
{
    struct finite_field test_field;
    test_field.characteristic = 7;
    test_field.power = 4;
    test_field.primative_in_power_n = 228;

    struct extended_polynomial true_answer;
    needed_type coefficients[17] = {4, 0, 0, 6, 4, 3, 5, 2, 6, 3, 4, 6, 0, 2, 4, 2, 1};
    true_answer.degree = 17;
    true_answer.coefs = coefficients;



    struct extended_polynomial* test_result = construct_generator_polynomial(&test_field, 2);
    int result = equal_polynomials(&true_answer, test_result);
    
    if(result != 1){
        printf("Test n 3 failed\n");
        print_extended_polynomial(test_result);
        print_extended_polynomial(&true_answer);
        printf("----------------------------\n");
    }
    
    free_extended_polynomial(test_result);

    return result;
}
*/

int
temp_test1()
{
    struct finite_field test_field;
    test_field.characteristic = 2;
    test_field.power = 5;
    test_field.primative_in_power_n = 5;

    int result = 0;

    uint64_t a1[MAX_DEGREE] = {1, 1, 0, 1}, c[MAX_DEGREE], d[MAX_DEGREE], q[MAX_DEGREE] = {0, 1}, r[MAX_DEGREE] = {1};
    uint64_t a2[MAX_DEGREE] = {1, 1};

    memset(a1 + 4, 0, (MAX_DEGREE - 4) * sizeof(uint64_t));
    memset(a2 + 2, 0, (MAX_DEGREE - 2) * sizeof(uint64_t));

    memset(q + 2, 0, (MAX_DEGREE - 2) * sizeof(uint64_t));
    memset(r + 1, 0, (MAX_DEGREE - 1) * sizeof(uint64_t));



    divide_polynomials_with_remainder(&test_field, a1, a1, c, d);

    if(equal_polynomials(q, c) && equal_polynomials(r, d)) result = 1;

    return result;
}

int
temp_test2()
{
    struct finite_field test_field;
    test_field.characteristic = 2;
    test_field.power = 5;
    test_field.primative_in_power_n = 5;
    uint64_t temp;
    temp = find_primitive_in_power(&test_field, -7);
    return multiply_in_field(&test_field, temp, 16 + 8 + 4 + 2) == 19;
}
    
int
temp_test3()
{
    int result = 0;
    struct finite_field* res = construct_gf_p_m(2, 5);
    if(res->power == 5 && res->primative_in_power_n == 5){
        
        result = 1;
    }
    if(result != 1){
        printf("TEST3\n");
    }
    free(res);
    return result;
}
/*
int
encode_test1()
{
    printf("Started encode test\n");
    struct bch_code* bch = init_bch(2, 5, 3, 0);

    uint64_t message[MAX_DEGREE] = {1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1};
    memset(message + 15, 0, (MAX_DEGREE - 15) * sizeof(uint64_t));
    
    uint64_t result[MAX_DEGREE], res[MAX_DEGREE];
    encode_bch(bch, message, result);
    decode_bch(bch, result, res);
    uint64_t decoded = extended_polynomial_to_polynomial(bch->field, res);
    int test_result = (decoded == 29583);
    
    
    if(test_result != 1){
        printf("Encode test failed\n");
        print_extended_polynomial(result);
    }
    
    
    free_bch_code(bch);
    
    return test_result;
}
*/

int
syndrom_test()
{
    printf("Started syndrome test\n");
    struct finite_field test_field;
    test_field.characteristic = 2;
    test_field.power = 5;
    test_field.primative_in_power_n = 5; 
    
    uint64_t correct_syndrome[6] = {3, 5, 9, 17, 4, 11};

    uint64_t encoded[MAX_DEGREE] = {1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1};
    memset(encoded + 31, 0, (MAX_DEGREE - 31) * sizeof(uint64_t));

    uint64_t syndrome[MAX_DEGREE];
    calculate_syndrome(&test_field, encoded, 3, syndrome);
    
    int result = 1;
    for(int i = 0; i < 6; ++i){
        result &= (correct_syndrome[i] == syndrome[i]);
    }

    if(result != 1){
        printf("Syndrome test failed\n");
        printf("----------------------------\n");
    }
    printf("Ended syndrome test\n");

    return result;
}

int inverse_test()
{
    struct finite_field test_field;
    test_field.characteristic = 2;
    test_field.power = 5;
    test_field.primative_in_power_n = 4 + 1; 
    uint64_t result = construct_inverse_element_multiply(&test_field, 3);
    if( result != 28){
        printf("Inverse test failed\n");
        printf("Result %lu\n", result);
        return 0;
        }
        return 1;
        
}
int
locator_test()
{
    struct finite_field test_field;
    test_field.characteristic = 2;
    test_field.power = 4;
    test_field.primative_in_power_n = 3; 

    uint64_t correct_locator[MAX_DEGREE] = {1, 3, 2};
    memset(correct_locator + 3, 0, (MAX_DEGREE - 3) * sizeof(uint64_t));
    
    uint64_t syndome[4] = {3, 5, 9, 2};

    uint64_t locator[MAX_DEGREE];
    construct_locator_polynomial(&test_field, syndome, 2, locator);

    int result = equal_polynomials(locator, correct_locator);


    if(result != 1){
        printf("Locator test failed\n");
        print_extended_polynomial(locator);
    }

    return result;
}

int
chien_test(){

    struct finite_field test_field;
    test_field.characteristic = 2;
    test_field.power = 5;
    test_field.primative_in_power_n = 5; 
    
    int result = 1;

    int true_errors[3] = {0, 1, -1};

    uint64_t encoded[MAX_DEGREE] = {1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1};
    memset(encoded + 31, 0, (MAX_DEGREE - 31) * sizeof(uint64_t));


    uint64_t syndrome[MAX_DEGREE];
    calculate_syndrome(&test_field, encoded, 3, syndrome);




    uint64_t locator[MAX_DEGREE];
    construct_locator_polynomial(&test_field, syndrome, 3, locator);

    int errors[MAX_DEGREE];
    chien_search(&test_field, locator, 3, errors);

    for(int i = 0; i < 3; ++i){
        result &= (errors[i] == true_errors[i]);   
    }

    return result;
}


int
decode_test(){
    struct bch_code* bch = init_bch(2, 5, 3, 0);

    int result = 1;

    uint64_t encoded[MAX_DEGREE] = {1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1}, res[MAX_DEGREE];
    memset(encoded + 31, 0, (MAX_DEGREE - 31) * sizeof(uint64_t));
    
    decode_bch(bch, encoded, res);
    
    uint64_t data = extended_polynomial_to_polynomial(bch->field, res), true_data = 45967;


    result = (data == 45967);
    if(result != 1){
        printf("Decoding failed\n");
        printf("Data recovered: %lu, data original: %lu\n", data, true_data);
        return 0;
    }
    free_bch_code(bch);
    return 1;
}
/*
int
encode_test2(){
    struct bch_code bch;
    struct finite_field* field = malloc(sizeof(struct finite_field));
    field->characteristic = 2;
    field->power = 4;
    field->primative_in_power_n = 3;

    bch.field = field;
    bch.number_of_errors = 2;
    bch.generator = construct_generator_polynomial(field, 2);
    print_extended_polynomial(bch.generator);

    bch.data_length = 16 - 1 - bch.generator->degree;

    

    return 1;
}
*/