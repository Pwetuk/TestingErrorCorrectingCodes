struct field{
    unsigned int characteristic;
    unsigned long long int power;
    unsigned long long int primative_in_power_n;
};

struct field*
generate_field(unsigned long long int power);

unsigned long long int
add_in_field(struct field* field, unsigned long long int a, unsigned long long int b);

unsigned long long int
multiply_in_field(struct field* field, unsigned long long int a, unsigned long long int b);

unsigned long long int
find_primitive_in_power(struct field* field, unsigned long long int power);