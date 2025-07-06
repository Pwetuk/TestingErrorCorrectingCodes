CC = gcc
CFLAGS = -Wall -Wextra -Ipaths 
OBj_DIR = obj

OBJS = $(OBj_DIR)/finite_field.o $(OBj_DIR)/polynoms.o $(OBj_DIR)/main.o $(OBj_DIR)/bch_tests.o

create_dir:
		@mkdir -p $(OBj_DIR)

all: create_dir main.out

$(OBj_DIR)/finite_field.o: finite_field.c finite_field.h
		$(CC) $(CFLAGS) -c -o $@ $<

$(OBj_DIR)/polynoms.o: polynoms.c $(OBj_DIR)/finite_field.o polynoms.h
		$(CC) $(CFLAGS) -c -o $@ $<

$(OBj_DIR)/bch_tests.o: bch_tests.c bch_tests.h $(OBj_DIR)/finite_field.o $(OBj_DIR)/polynoms.o
		$(CC) $(CFLAGS) -c -o $@ $<

$(OBj_DIR)/main.o: main.c $(OBj_DIR)/polynoms.o $(OBj_DIR)/finite_field.o $(OBj_DIR)/bch_tests.o
		$(CC) $(CFLAGS) -c -o $@ $<

main.out: $(OBJS)
		$(CC) $(CFLAGS) -o $@ $^

clean:
		rm -f $(OBj_DIR)/*.o main.out

make_clean: all clean
