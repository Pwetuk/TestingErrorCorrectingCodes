CC = gcc
CFLAGS = -Wall -Wextra -Ipaths -fstack-usage
OBj_DIR = obj

OBJS = $(OBj_DIR)/finite_field.o $(OBj_DIR)/polynomials.o $(OBj_DIR)/main.o $(OBj_DIR)/bch_tests.o $(OBj_DIR)/bch_encode.o

create_dir:
		@mkdir -p $(OBj_DIR)

all: create_dir main.out

$(OBj_DIR)/finite_field.o: finite_field.c finite_field.h
		$(CC) $(CFLAGS) -c -o $@ $<  -lm

$(OBj_DIR)/polynomials.o: polynomials.c $(OBj_DIR)/finite_field.o polynomials.h
		$(CC) $(CFLAGS) -c -o $@ $< -lm

$(OBj_DIR)/bch_encode.o: bch_encode.c bch_encode.h $(OBj_DIR)/polynomials.o $(OBj_DIR)/finite_field.o
		$(CC) $(CFLAGS) -c -o $@ $< -lm

$(OBj_DIR)/bch_tests.o: bch_tests.c bch_tests.h $(OBj_DIR)/finite_field.o $(OBj_DIR)/polynomials.o $(OBj_DIR)/bch_encode.o
		$(CC) $(CFLAGS) -c -o $@ $< -lm

$(OBj_DIR)/main.o: main.c $(OBj_DIR)/polynomials.o $(OBj_DIR)/finite_field.o $(OBj_DIR)/bch_encode.o $(OBj_DIR)/bch_tests.o
		$(CC) $(CFLAGS) -c -o $@ $< -lm

main.out: $(OBJS)
		$(CC) $(CFLAGS) -o $@ $^ -lm

clean:
		rm -f $(OBj_DIR)/*.o main.out

make_clean: all clean
