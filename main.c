#include <stdio.h>
#include <string.h>

#include "bch_tests.h"


int main(int argc, char** argv){
    int tests_result = 0;
    if(argc > 1 && strcmp(argv[1], "test") == 0){
        tests_result = run_tests();
    }else if(argc > 1 && strcmp(argv[1], "ttest") == 0){
        tests_result = temp_tests();
    }
    if(tests_result == 1){
        printf("All tests succeded\n");
    }else{
        printf("Something went wrong\n");
    }
    
    return 0;
}