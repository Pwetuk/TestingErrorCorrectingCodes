#include <stdio.h>
#include "finite_field.h"

unsigned long long
add_in_binary_field(struct binary_field*, unsigned long long a, unsigned long long b)
{
    unsigned long long c = a ^ b;
    return c;
}

unsigned long long
multiply_in_binary_field(struct binary_field* field, unsigned long long a, unsigned long long b)
{
    unsigned long long result = 0;
    unsigned long long max_val = 1 << field->power;
    while (a) {
        if(a % 2 == 1){
            result ^= b;
        }
        a >>= 1;
        b <<= 1;
        if(b >= max_val){
            b ^= max_val;
            b ^= field->primative_in_power_n;
        }
    }
    return result;
}

unsigned long long
find_primitive_in_power(struct binary_field* field, unsigned long long power)
{
    if(power < field->power){
        return (1 << power);
    }else if(power % 2 == 0){
        unsigned long long a = find_primitive_in_power(field, power/ 2);
        return  multiply_in_binary_field(field, a, a);
    }
    return multiply_in_binary_field(field, find_primitive_in_power(field, power - 1), 2);
}