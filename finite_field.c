#include <stdio.h>
#include <math.h>

#include "finite_field.h"

unsigned long long
pow_placeholder(unsigned long long int base, int power)
{
    if(power == 0){
        return 1;
    }
    if(power == 1){
        return base;
    }
    if(power % 2 == 0){
        unsigned long long int a = pow_placeholder(base, power / 2);
        return a * a;
    }
    return base * pow_placeholder(base, power - 1);
}

unsigned long long
add_in_field(struct finite_field* field, unsigned long long a, unsigned long long b)
{
    if(field->characteristic == 2) return a ^ b;
    unsigned long long int a_copy = a, b_copy = b, result = 0;
    int i = 0;
    while (a_copy != 0)
    {
        result += pow((double)field->characteristic, (double)i) * ((a_copy + b_copy) % field->characteristic);
        a_copy /= (field->characteristic);
        b_copy /= (field->characteristic);
        ++i;
    }
    result += pow(field->characteristic, i) * ((a_copy + b_copy));
    return result;
}

unsigned long long
multiply_in_field(struct finite_field* field, unsigned long long a, unsigned long long b)
{   
    if(field->characteristic == 2){
        unsigned long long result = 0;
        unsigned long long max_val = 1 << field->power;
        unsigned long long a_copy = a, b_copy = b;
        while (a_copy) {
            if(a_copy % 2 == 1){
                result ^= b_copy;
            }
            a_copy >>= 1;
            b_copy <<= 1;
            if(b_copy >= max_val){
                b_copy ^= max_val;
                b_copy ^= field->primative_in_power_n;
            }
        }
        return result;
    }
    unsigned long long result = 0;
    unsigned long long max_val = pow(field->characteristic, field->power);
    unsigned long long a_copy = a, b_copy = b;
    while(a_copy != 0){
        if(a_copy % field->characteristic != 0){
            for(int i = 0; (unsigned long long)i < a_copy % field->characteristic; ++i){
                result = add_in_field(field, result, b_copy);
            }
        }
        a_copy /= field->characteristic;
        b_copy *= field->characteristic;
        while(b_copy >= max_val){
            b_copy -= max_val;
            b_copy = add_in_field(field, b_copy, field->primative_in_power_n);
        }
    }
    return result;
}

unsigned long long
find_primitive_in_power(struct finite_field* field, unsigned long long power)
{
    if(power < field->power){
        if(field->characteristic == 2){
            return (1 << power);
        }
        return pow(field->characteristic, power);
    }else if(power % 2 == 0){
        unsigned long long a = find_primitive_in_power(field, power/ 2);
        return  multiply_in_field(field, a, a);
    }
    return multiply_in_field(field, find_primitive_in_power(field, power - 1), field->characteristic);
}

unsigned long long
construct_inreverse_element(struct finite_field* field, unsigned long long element)
{
    if(field->characteristic == 2) return element;
    unsigned long long element_copy = element;
    unsigned long long result = 0;
    int i = 0;
    while(element_copy != 0){
        if(element_copy % field->characteristic != 0){
            result += (field->characteristic - (element_copy % field->characteristic)) * pow(field->characteristic, i);
        }
        element_copy /= field->characteristic;
        ++i;
    }
    return result;
}