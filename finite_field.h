struct binary_field{
    unsigned long long int power;
    unsigned long long int primative_in_power_n;
};


struct binary_field*
generate_binary_field(unsigned long long int power);

unsigned long long int
add_in_binary_field(struct binary_field* field, unsigned long long int a, unsigned long long int b);

unsigned long long int
multiply_in_binary_field(struct binary_field* field, unsigned long long int a, unsigned long long int b);

unsigned long long int
find_primitive_in_power(struct binary_field* field, unsigned long long int power);