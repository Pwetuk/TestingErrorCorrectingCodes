#include <stdio.h>
#include <stdlib.h>

#include "bch_tests.h"
#include "finite_field.h"
#include "polynomials.h"

int 
run_tests()
{
    return test_in_characteristics_2_first() &
        test_in_characteristics_2_second() &
        test_in_characteristics_2_third() &
        test_in_characteristics_n_first() &
        test_in_characteristics_n_second() &
        test_in_characteristics_n_third();
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
    printf("----------------\n");
    printf("%lld\n", test_result);
    print_polynomial(&test_field, test_result);
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
    printf("----------------\n");
    printf("%lld\n", test_result);
    print_polynomial(&test_field, test_result);
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
    printf("----------------\n");
    printf("%lld\n", test_result);
    print_polynomial(&test_field, test_result);
    return test_result == true_answer;

}
