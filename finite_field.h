struct finite_field{
    unsigned int characteristic;
    unsigned long long int power;
    unsigned long long int primative_in_power_n;
};

struct finite_field*
generate_field(unsigned long long int power);

unsigned long long int
add_in_field(struct finite_field* field, unsigned long long int a, unsigned long long int b);

unsigned long long int
multiply_in_field(struct finite_field* field, unsigned long long int a, unsigned long long int b);

unsigned long long int
find_primitive_in_power(struct finite_field* field, unsigned long long int power);


unsigned long long
pow_placeholder(unsigned long long int base, int power);

unsigned long long
construct_inverse_element_add(struct finite_field* field, unsigned long long element);

unsigned long long
construct_inverse_element_multiply(struct finite_field* field, unsigned long long element);