#include <stdio.h>
#include <stdlib.h>

#include "bch_tests.h"
#include "finite_field.h"
#include "polynoms.h"

int 
run_tests()
{

    //return test1() | test2() | test4() | test5() | test6() | test7() | test8();
    return test8();
}


int
test1(){
    struct binary_field test_field;
    test_field.power = 3;
    test_field.primative_in_power_n = 3;
    printf("Result: %lld\n", multiply_in_binary_field(&test_field, 7, 2));
    printf("Result: %lld\n", add_in_binary_field(&test_field, 7, 2));
    return 0;
}

int test2(){
    struct binary_field test_field;
    test_field.power = 4;
    test_field.primative_in_power_n = 3;
    struct cyclotomic_cosets* to_free = generate_cosets_from_n_needed_elements(&test_field, 6);
    free_cosets_struct(to_free);
    return 0;
}

int test3(){
    struct binary_field test_field;
    test_field.power = 16;
    test_field.primative_in_power_n = 45;
    for(unsigned long long int i = 0; i < (1UL << 16); ++i){
        printf("x^%lld = %lld\n",i, find_primitive_in_power(&test_field, i));
    }
    return 0;
}

int test4(){
    struct binary_field test_field;
    test_field.power = 4;
    test_field.primative_in_power_n = 3;
    polynoms a;
    polynoms b;
    a = 2 + 4 + 8 + (1 << 16);
    b = 2 + 4;
    /*
    print_polynom(&test_field, a);
    print_polynom(&test_field, b);
    printf("---------------------------------\n");
    */
    print_polynom(&test_field, add_two_polynoms(&test_field, a, b));
    print_polynom(&test_field, multiply_two_polynoms(&test_field, a, b));
    return 0;
}

int
test5()
{
    struct binary_field test_field;
    test_field.power = 4;
    test_field.primative_in_power_n = 3;
    int coefs1[4] = {3, 0, 15, 1};
    int coefs2[5] = {2, 0, 2, 0, 1};
    struct extended_polynom *a, *b, *c;
    a = construct_extended_polynom_from_coefs(coefs1, 4);
    b = construct_extended_polynom_from_coefs(coefs2, 5);

    c = add_two_extended_polynoms(&test_field, a, b);

    print_extended_polynom(c);



    free_extended_polynom(a);
    free_extended_polynom(b);
    free_extended_polynom(c);
    return 0;
}

int
test6()
{
    printf("Multiplication test\n");
    struct binary_field test_field;
    test_field.power = 4;
    test_field.primative_in_power_n = 3;
    int coefs1[4] = {3, 0, 15, 1};
    int coefs2[5] = {2, 0, 2, 0, 1};
    struct extended_polynom *a, *b, *c;
    a = construct_extended_polynom_from_coefs(coefs1, 4);
    b = construct_extended_polynom_from_coefs(coefs2, 5);

    c = multiplty_two_extended_polynoms(&test_field, a, b);

    print_extended_polynom(c);



    free_extended_polynom(a);
    free_extended_polynom(b);
    free_extended_polynom(c);
    return 0;
}

int
test7()
{
    printf("Minimal polynomial test");
    struct binary_field test_field;
    test_field.power = 4;
    test_field.primative_in_power_n = 3;
    struct cyclotomic_cosets* to_free = generate_cosets_from_n_needed_elements(&test_field, 8);


    for(int i = 0; i < to_free->number_of_cosets; ++i){
        polynoms minimal_polinomial = construct_minimal_polynom_from_coset(&test_field, to_free->cosets[i]);
        
        print_polynom(&test_field, minimal_polinomial);
    }

    free_cosets_struct(to_free);
    return 0;
}

int
test8()
{
    printf("Generator polynomial test\n");
    struct binary_field test_field;
    test_field.power = 6;
    test_field.primative_in_power_n = 3;

    printf("Generator polynomial: %lld\n", construct_generator_polynomial(&test_field, 2));

    return 0;
}