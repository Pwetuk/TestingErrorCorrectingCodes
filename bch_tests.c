#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bch_tests.h"

int 
run_tests()
{
    return test_in_characteristics_2_first() &
        test_in_characteristics_2_second() &
        test_in_characteristics_2_third() &
        test_in_characteristics_n_first() &
        test_in_characteristics_n_second() &
        test_in_characteristics_n_third() &
        syndrom_test() &
        encode_test1() &
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
    if(temp_test1() != 1){
        printf("Test1\n");
        return 0;
    }
    if(temp_test2() != 1){
        printf("Test2\n");
        return 0;
    }
    if(temp_test3() != 1){
        printf("Test3\n");
        return 0;
    }
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

    struct extended_polynomial true_answer;
    needed_type coefficients[11] = {1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1};
    true_answer.degree = 11;
    true_answer.coefs = coefficients;



    struct extended_polynomial* test_result = construct_generator_polynomial(&test_field, 2);
    int result = equal_polynomials(&true_answer, test_result);


    
    if(result != 1){
        printf("Test1 failed\n");
        print_extended_polynomial(test_result);
        print_extended_polynomial(&true_answer);
        printf("----------------------------\n");
    }
    
    free_extended_polynomial(test_result);

    return result;
}


int
test_in_characteristics_2_second()
{
    struct finite_field test_field;
    test_field.characteristic = 2;
    test_field.power = 6;
    test_field.primative_in_power_n = 3;


    struct extended_polynomial true_answer;
    needed_type coefficients[13] = {1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1};
    true_answer.degree = 13;
    true_answer.coefs = coefficients;



    struct extended_polynomial* test_result = construct_generator_polynomial(&test_field, 2);
    int result = equal_polynomials(&true_answer, test_result);
    
    if(result != 1){
        printf("Test2 failed\n");
        print_extended_polynomial(test_result);
        print_extended_polynomial(&true_answer);
        printf("----------------------------\n");
    }
    
    free_extended_polynomial(test_result);

    return result;
}

int
test_in_characteristics_2_third()
{
    struct finite_field test_field;
    test_field.characteristic = 2;
    test_field.power = 4;
    test_field.primative_in_power_n = 3;

    struct extended_polynomial true_answer;
    needed_type coefficients[11] = {1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1};
    true_answer.degree = 11;
    true_answer.coefs = coefficients;



    struct extended_polynomial* test_result = construct_generator_polynomial(&test_field, 3);
    int result = equal_polynomials(&true_answer, test_result);
    
    if(result != 1){
        printf("Test3 failed\n");
        print_extended_polynomial(test_result);
        print_extended_polynomial(&true_answer);
        printf("----------------------------\n");
    }
    
    free_extended_polynomial(test_result);

    return result;
}

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


int
temp_test1()
{
    struct finite_field test_field;
    test_field.characteristic = 2;
    test_field.power = 5;
    test_field.primative_in_power_n = 5;

    int result = 0;


    struct extended_polynomial a, b, *c, *d;

    needed_type a1[16] = {1, 1, 0, 1};
    needed_type a2[16] = {1, 1};

    a.degree = 4;
    a.coefs = a1;
    
    b.degree = 2;
    b.coefs = a2;

    divide_polynomials_with_remainder(&test_field, &a, &b, &c, &d);

    if(c->degree == 3 && c->coefs[0] == 0 && c->coefs[1] == 1 && c->coefs[2] == 1 && d->degree == 1 && d->coefs[0] == 1) result = 1;
    free_extended_polynomial(c);
    free_extended_polynomial(d);

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

int
encode_test1()
{
    struct bch_code* bch = init_bch(2, 5, 3, 0);

    struct extended_polynomial message, *res;

    needed_type coefficients_msg[15] = {1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1};
    message.degree = 15;
    message.coefs = coefficients_msg;
    
    struct extended_polynomial* result = encode_bch(bch, &message);
    res = decode_bch(bch, result);
    uint64_t decoded = extended_polynomial_to_polynomial(bch->field, res);
    free_extended_polynomial(res);
    int test_result = (decoded == 29583);


    if(test_result != 1){
        printf("Encode test failed\n");
        print_extended_polynomial(result);
    }
        
    free_extended_polynomial(result);
    free_bch_code(bch);

    return test_result;
}

int
syndrom_test()
{
    struct finite_field test_field;
    test_field.characteristic = 2;
    test_field.power = 5;
    test_field.primative_in_power_n = 5; 
    
    uint64_t correct_syndrome[6] = {3, 5, 9, 17, 4, 11};

    needed_type coefficients[31] = {1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1};
    struct extended_polynomial encoded;
    encoded.degree = 31;
    encoded.coefs = coefficients;


    uint64_t *syndrome = calculate_syndrome(&test_field, &encoded, 3);
    
    int result = 1;
    for(int i = 0; i < 6; ++i){
        result &= (correct_syndrome[i] == syndrome[i]);
    }

    free(syndrome);

    if(result != 1){
        printf("Syndrome test failed\n");
        printf("----------------------------\n");
    }

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


    needed_type coefs[3] = {1, 3, 2};
    struct extended_polynomial* correct_locator = construct_extended_polynomial_from_coefs(coefs, 3);
    
    uint64_t syndome[4] = {3, 5, 9, 2};

    struct extended_polynomial* locator = construct_locator_polynomial(&test_field, syndome, 2);

    int result = equal_polynomials(locator, correct_locator);


    if(result != 1){
        printf("Locator test failed\n");
        print_extended_polynomial(locator);
    }
    
    free_extended_polynomial(locator);
    free_extended_polynomial(correct_locator);

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

    needed_type coefficients[31] = {1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1};
    struct extended_polynomial encoded;
    encoded.degree = 31;
    encoded.coefs = coefficients;


    uint64_t *syndrome = calculate_syndrome(&test_field, &encoded, 3);




    struct extended_polynomial* locator = construct_locator_polynomial(&test_field, syndrome, 3);

    int *errors = chien_search(&test_field, locator, 3);

    for(int i = 0; i < 3; ++i){
        result &= (errors[i] == true_errors[i]);   
    }

    free(syndrome);
    free(errors);
    free_extended_polynomial(locator);

    return result;
}


int
decode_test(){
    struct bch_code* bch = init_bch(2, 5, 3, 0);

    int result = 1;

    needed_type coefficients[31] = {1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1};
    struct extended_polynomial encoded, *res;
    encoded.degree = 31;
    encoded.coefs = coefficients;
    
    res = decode_bch(bch, &encoded);
    
    uint64_t data = extended_polynomial_to_polynomial(bch->field, res), true_data = 45967;

    free_extended_polynomial(res);

    result = (data == 45967);
    if(result != 1){
        printf("Decoding failed\n");
        printf("Data recovered: %lu, data original: %lu\n", data, true_data);
        return 0;
    }
    free_bch_code(bch);
    return 1;
}

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