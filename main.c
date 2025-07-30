#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <math.h>

#include "bch_encode.h"

#define PAGE_SIZE 512
#define NUM_PAGES 1

void
put_error(uint64_t *arr, uint64_t len_arr){
    for(uint64_t i = 0; i < len_arr; ++i){
        arr[i] ^= (1UL << (rand() % 16));
    }
}

static uint64_t read_total_jiffies(void) {
    FILE *f = fopen("/proc/stat","r");
    uint64_t user, nice, sys, idle;
    fscanf(f, "cpu  %lu %lu %lu %lu", &user, &nice, &sys, &idle);
    fclose(f);
    return user + nice + sys + idle;
}

static uint64_t read_proc_jiffies(void) {
    FILE *f = fopen("/proc/self/stat","r");
    for(int i=0; i<13; i++) fscanf(f, "%*s");
    uint64_t utime, stime;
    fscanf(f, "%lu %lu", &utime, &stime);
    fclose(f);
    return utime + stime;
}

uint64_t*
prepare_arr(uint64_t input_arr[NUM_PAGES * PAGE_SIZE], int data_len, int *len_res){
    printf("Dt_len: %d\n", data_len);
    *len_res = NUM_PAGES * PAGE_SIZE * sizeof(uint64_t) / data_len + 1;
    int pos = 0;
    uint64_t *result = malloc(sizeof(uint64_t) * (*len_res));
    uint64_t buff = 0, el;
    int deg = 0;
    for(uint64_t i = 0; i < NUM_PAGES * PAGE_SIZE; ++i){
        el = input_arr[i];
        for(size_t j = 0; j < sizeof(uint64_t); ++j){
            if((el & 1) == 1){
                buff += (1UL << deg);
            }
            ++deg;
            if(deg >= data_len){
                deg = 0;
                result[pos] = buff;
                ++pos;
                buff = 0;
                
            }
            el >>= 1;
        }
    }
    result[pos] = buff;
    return result;
}

void
measure_bch(int p, int power, int t){
    clock_t start = clock();
    struct bch_code* bch_test = init_bch(p, power, t, 0);
    clock_t end = clock();
    double cpu_time_sec = 1000 * (double)(end - start) / CLOCKS_PER_SEC;
    printf("Init bch: %f seconds, %d\n", cpu_time_sec, 0);
    printf("Gen: %d\n", get_degree(bch_test->generator));

    int len_res;
    uint64_t to_encode, encoding_poly[MAX_DEGREE], encoded[MAX_DEGREE], decoded[MAX_DEGREE], to_decode, decoding_poly[MAX_DEGREE], temp;
    memset(encoding_poly, 0, MAX_DEGREE * sizeof(uint64_t));
    memset(decoding_poly, 0, MAX_DEGREE * sizeof(uint64_t));

    uint64_t input_arr[PAGE_SIZE * NUM_PAGES], *encoded_arr, *decoded_arr;
    for(int i = 0; i < PAGE_SIZE * NUM_PAGES; ++i){
        input_arr[i] = (((uint64_t) rand() << 32) | (uint64_t) rand()) % (1UL << (bch_test->data_length - 1));
    }
    
    start = clock();
    uint64_t *prepared = prepare_arr(input_arr, bch_test->data_length, &len_res);
    end = clock();
    cpu_time_sec = 1000 * (double)(end - start) / CLOCKS_PER_SEC;
    printf("Normalize data: %f seconds, %d\n", cpu_time_sec, 0);

    encoded_arr = malloc(sizeof(uint64_t) * len_res);
    decoded_arr = malloc(sizeof(uint64_t) * len_res);

    uint64_t sys1 = read_total_jiffies(), proc1 = read_proc_jiffies();
    start = clock();

    for(int i = 0; i < len_res; ++i){
        to_encode = prepared[i];
        for(int i = 0; i < bch_test->data_length; ++i){
            encoding_poly[i] = to_encode % bch_test->field->characteristic;
            to_encode >>= 1;
        }
        encode_bch(bch_test, encoding_poly, encoded);
        encoded_arr[i] = extended_polynomial_to_polynomial(bch_test->field, encoded);
    }
    end = clock();
    cpu_time_sec = 1000 * (double)(end - start) / CLOCKS_PER_SEC;
    printf("Encoded: %f seconds, %d\n", cpu_time_sec, 0);

    for(int n = 0; n <= bch_test->number_of_errors; ++n){
        start = clock();
        for(int i = 0; i < len_res; ++i){
            to_decode = encoded_arr[i];
            for(int i = 0; i < bch_test->n; ++i){
                decoding_poly[i] = to_decode % bch_test->field->characteristic;
                to_decode >>= 1;
            }
            decode_bch(bch_test, decoding_poly, decoded);
            decoded_arr[i] = extended_polynomial_to_polynomial(bch_test->field, decoded);
            if(decoded_arr[i] != prepared[i]) printf("PANIC");
        }
        end = clock();
        double cpu_time_sec = 1000 * (double)(end - start) / CLOCKS_PER_SEC;
        printf("Decoded with %d errors: %f seconds\n", n, cpu_time_sec);
        put_error(encoded_arr, len_res);
    }
    free_bch_code(bch_test);
    free(prepared);
    free(encoded_arr);
    free(decoded_arr);
}

int main(int argc, char** argv){
    
    
    measure_bch(2, 5, 2);
    measure_bch(2, 5, 3);
    measure_bch(2, 6, 3);
    

    return 0;
}