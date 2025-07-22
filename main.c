#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <sys/resource.h>

#include "bch_tests.h"

#define PAGE_SIZE 512
#define NUM_PAGES 1

void
put_error(uint64_t *arr, uint64_t len_arr){
    for(uint64_t i = 0; i < len_arr; ++i){
        arr[i] ^= (1UL << rand());
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

int main(int argc, char** argv){
    int tests_result = 0;
    int was_testing = 1;
    if(argc > 1 && strcmp(argv[1], "test") == 0){
        tests_result = run_tests();
    }else if(argc > 1 && strcmp(argv[1], "ttest") == 0){
        tests_result = temp_tests();
    }else{
        was_testing = 0;
    }

    if(was_testing == 1){
        if(tests_result == 1){
            printf("All tests succeded\n");
        }else{
            printf("Something went wrong\n");
        }
        return 0;
    }
    /*
    unsigned int characteristic;
    uint64_t power;
    uint64_t res_len, decoded_len;
    int number_of_errors;
    
    
    printf("Input characteristic, power, number of errors polynomial\n");
    scanf("%u %lu %d", &characteristic, &power, &number_of_errors);
    
    clock_t start = clock();
    struct bch_code* bch_test = init_bch(characteristic, power, number_of_errors, 0);
    clock_t end = clock();
    double cpu_time_sec = 1000 * (double)(end - start) / CLOCKS_PER_SEC;
    
    printf("Inited bch:%f seconds, %lu\n", cpu_time_sec, find_primitive_in_power(bch_test->field, 1023));
    print_extended_polynomial(bch_test->generator);
    
    uint64_t input_arr[PAGE_SIZE * NUM_PAGES];
    
    for(int i = 0; i < PAGE_SIZE * NUM_PAGES; ++i){
        input_arr[i] = ((uint64_t) rand() << 32) | rand();
    }
    printf("\n");
    uint64_t sys1 = read_total_jiffies(), proc1 = read_proc_jiffies();
    start = clock();
    uint64_t *encoded = bch_encode_arr(bch_test, input_arr, PAGE_SIZE * NUM_PAGES, &res_len);
    end = clock();
    uint64_t sys2 = read_total_jiffies(), proc2 = read_proc_jiffies();
    double cpu_pct = 100.0 * (proc2 - proc1) / (double)(sys2 - sys1);
    printf("Process CPU load: %.2f%%\n", cpu_pct);
    
    cpu_time_sec = 1000 * (double)(end - start) / CLOCKS_PER_SEC;
    printf("Encoded: %f seconds, %d\n", cpu_time_sec, res_len);
    //printf("\n\n----------------------\n\n");
    //for(int i = 0; i < NUM_PAGES * PAGE_SIZE; ++i){
        //    uint64_t el = input_arr[i];
        //    for(int j = 0; j < 64; ++j){
            //        printf("%lu", el % 2);
            //        el >>= 1;
            //    }
            //    
            //}
            //printf("\n\n----------------------\n\n");
            
            for(int i = 0; i < 2; ++i){
                start = clock();
                uint64_t *decoded = bch_decode_arr(bch_test, encoded, res_len, &decoded_len);
                end = clock();
                cpu_time_sec = 1000 * (double)(end - start) / CLOCKS_PER_SEC;
                printf("Decoded %f seconds with %d errors\n", cpu_time_sec, i);
                put_error(encoded, res_len);
                int result = 1;
                for(uint64_t i = 0; i < NUM_PAGES * PAGE_SIZE; ++i){
                    result &= (decoded[i] == input_arr[i]);
                    if(result == 0){
                        printf("Expected: %lu, got %lu in ind: %lu\n", input_arr[i], decoded[i], i);
                        break;
                    }
                }
                
                printf("Result is: %d\n", result);
                
                free(decoded);
            }
            
    free(encoded);
    */
    /*
    start = clock();
    for(int i = 0; i < PAGE_SIZE * NUM_PAGES; ++i){
        to_encode = input_arr[i];
        to_encode_poly = construct_polynomial_from_field_element(bch_test->field, to_encode);
    
        struct extended_polynomial* result = encode_bch(bch_test, to_encode_poly);
        input_arr[i] = extended_polynomial_to_polynomial(bch_test->field, result);
        free_extended_polynomial(result);
    }
    end = clock();
    cpu_time_sec = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Encoded: %f seconds\n", cpu_time_sec);

    uint64_t to_decode;

    for(int num_err = 0; num_err < number_of_errors; ++num_err){
        
        start = clock();
        for(int i = 0; i < PAGE_SIZE * NUM_PAGES; ++i){
            to_encode = input_arr[i];
            to_encode_poly = construct_polynomial_from_field_element(bch_test->field, to_encode + 1);
        
            decode_bch(bch_test, to_encode_poly);
            free_extended_polynomial(to_encode_poly);
        }
        end = clock();
        cpu_time_sec = (double)(end - start) / CLOCKS_PER_SEC;
        printf("Decoded %f seconds with %d errors\n", cpu_time_sec, num_err);
        
        for(int i = 0; i < PAGE_SIZE * NUM_PAGES; ++i){
            int err_pos = rand() % (bch_test->data_length + bch_test->generator->degree - 2);
            input_arr[i] ^= (1ULL << err_pos);
        }

    }

    
    free_bch_code(bch_test);
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    printf("Peak memory usage: %ld\n", ru.ru_maxrss);
    */

    return 0;
}