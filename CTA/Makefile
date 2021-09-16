CC=gcc
CFLAGS=-O3 -march=native -Wall

OBJS=fastrandombytes.o gaussian_avx.o poly_mult.o poly_q.o randombytes.o sammat.o sample.o merkle_tree_commit.o merkle_tree_response.o merkle_tree_verify.o test.o

PROGRAM=cta

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CC) $(OBJS) -o $(PROGRAM)

clean:;
	rm -f *.o $(PROGRAM)
