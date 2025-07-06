#include <stdio.h>
#include <string.h>

#include "bch_tests.h"


int main(int argc, char** argv){
    if(argc > 1 && strcmp(argv[1], "test") == 0){
        int tests_result = run_tests();
        if(tests_result == 0){
            printf("All tests succeded\n");
        }else{
            printf("Something went wrong\n");
        }
    }
    
    return 0;
}