#include <math.h>

#include "bch_encode.h"


polynomial
encode_bch(struct finite_field* field, polynomial generator, polynomial data, int n, int degree_data)
{
    polynomial result;
    data *= pow(field->characteristic, n - degree_data);
    polynomial quotient, remainder;
    divide_polynomials_with_remainder(field, data, generator, &quotient, &remainder);
    result = difference_of_two_polynomials(field, data, remainder);
    return result;
}