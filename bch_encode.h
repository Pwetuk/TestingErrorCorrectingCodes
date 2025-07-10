#include "polynomials.h"

polynomial
encode_bch(struct finite_field* field, polynomial generator, polynomial data, int n, int degree_data);