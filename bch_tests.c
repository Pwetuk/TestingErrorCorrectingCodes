#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bch_tests.h"
#include "bch_encode.h"

int 
run_tests()
{
    return test_in_characteristics_2_first() &
        test_in_characteristics_2_second() &
        test_in_characteristics_2_third() &
        test_in_characteristics_n_first() &
        test_in_characteristics_n_second() &
        test_in_characteristics_n_third() &
        encode_test1();
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
    if(temp_test4() != 1){
        printf("Test4\n");
        return 0;
    }
    if(temp_test5() != 1){
        printf("Test1\n");
        return 0;
    }
    return 1;
}


int
test_in_characteristics_2_first()
{
    polynomial true_answer = 1897;
    struct finite_field test_field;
    test_field.characteristic = 2;
    test_field.power = 5;
    test_field.primative_in_power_n = 5;
    polynomial test_result = construct_generator_polynomial(&test_field, 2);

    return test_result == true_answer;
}


int
test_in_characteristics_2_second()
{
    polynomial true_answer = 5433;
    struct finite_field test_field;
    test_field.characteristic = 2;
    test_field.power = 6;
    test_field.primative_in_power_n = 3;
    polynomial test_result = construct_generator_polynomial(&test_field, 2);

    return test_result == true_answer;
}

int
test_in_characteristics_2_third()
{
    polynomial true_answer = 1335;
    struct finite_field test_field;
    test_field.characteristic = 2;
    test_field.power = 4;
    test_field.primative_in_power_n = 3;
    polynomial test_result = construct_generator_polynomial(&test_field, 3);

    return test_result == true_answer;
}

int
test_in_characteristics_n_first()
{
    polynomial true_answer = 19930607;
    struct finite_field test_field;
    test_field.characteristic = 3;
    test_field.power = 5;
    test_field.primative_in_power_n = 25;
    polynomial test_result = construct_generator_polynomial(&test_field, 2);
    return test_result == true_answer;

}


int
test_in_characteristics_n_second()
{
    polynomial true_answer = 462072526;
    struct finite_field test_field;
    test_field.characteristic = 5;
    test_field.power = 3;
    test_field.primative_in_power_n = 4 + 15;
    polynomial test_result = construct_generator_polynomial(&test_field, 2);
    return test_result == true_answer;

}


int
test_in_characteristics_n_third()
{
    polynomial true_answer = 45647875710196;
    struct finite_field test_field;
    test_field.characteristic = 7;
    test_field.power = 4;
    test_field.primative_in_power_n = 228;
    polynomial test_result = construct_generator_polynomial(&test_field, 2);
    return test_result == true_answer;

}


int
temp_test1()
{
    polynomial quotient = 2 + 4, remainder = 1;
    struct finite_field test_field;
    test_field.characteristic = 2;
    test_field.power = 5;
    test_field.primative_in_power_n = 5;
    polynomial a = 1 + 2 + 8, b = 1 + 2, c, d;

    divide_polynomials_with_remainder(&test_field, a, b, &c, &d);

    return (c == quotient) & (d == remainder);
}

int
temp_test2()
{
    polynomial quotient = 4, remainder = 1;
    struct finite_field test_field;
    test_field.characteristic = 2;
    test_field.power = 5;
    test_field.primative_in_power_n = 5;
    polynomial a = 1 + 4 + 16, b = 1 + 4, c, d;

    divide_polynomials_with_remainder(&test_field, a, b, &c, &d);

    return (c == quotient) & (d == remainder);
}

int
temp_test3()
{
    polynomial quotient = 2 + 2 * 3, remainder = 0;
    struct finite_field test_field;
    test_field.characteristic = 3;
    test_field.power = 5;
    test_field.primative_in_power_n = 5;
    polynomial a = 2 + 9, b = 1 + 6, c, d;

    divide_polynomials_with_remainder(&test_field, a, b, &c, &d);

    return (c == quotient) & (d == remainder);
}

int
temp_test4()
{
    polynomial quotient = 3 * 5, remainder = 2;
    struct finite_field test_field;
    test_field.characteristic = 5;
    test_field.power = 5;
    test_field.primative_in_power_n = 5;
    polynomial a = 125 * 3 + 4 * 5 + 2, b = 3 + 25, c, d;

    divide_polynomials_with_remainder(&test_field, a, b, &c, &d);

    return (c == quotient) & (d == remainder);
}

int
temp_test5()
{
    polynomial quotient = 5 + 4 * 7 + 4 * 49, remainder = 14;
    struct finite_field test_field;
    test_field.characteristic = 7;
    test_field.power = 5;
    test_field.primative_in_power_n = 5;
    polynomial a = 5 + 3 * 7 + 6 * 49 + 49 * 49, b = 1 + 5 * 7 + 2 * 49, c, d;

    divide_polynomials_with_remainder(&test_field, a, b, &c, &d);

    return (c == quotient) & (d == remainder);
}

int
encode_test1()
{
    struct finite_field test_field;
    test_field.characteristic = 2;
    test_field.power = 5;
    test_field.primative_in_power_n = 5;  
    printf("---------------\n");
    polynomial generator = construct_generator_polynomial(&test_field, 3);
    polynomial message = 45967;
    polynomial encoded = 1506255096;
    polynomial result = encode_bch(&test_field, generator, message, pow(test_field.characteristic, test_field.power) - 1, polynomial_degree(&test_field, message));
    printf("Encoded message: %lld\n", result);
    return (result == encoded);
}