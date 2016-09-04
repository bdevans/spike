CC     = gcc-4.2 -arch i386 -O3 -lm

CFLAGS  = -I/usr/lib/gcc/i686-apple-darwin9/4.2.1/include/ -I/usr/local/include/gsl/
LDFLAGS = -L/usr/lib/gcc/i686-apple-darwin9/4.2.1/ -L/usr/local/lib/

LIBS    = -lgomp -lgsl

Spike: main.o array_utils.o utils.o read_parameters.o rng.o spike.o
	$(CC) -o Spike $(CFLAGS) $(LDFLAGS) main.o array_utils.o utils.o read_parameters.o rng.o spike.o $(LIBS)

main.o: main.c
	$(CC) -fopenmp $(CFLAGS) -c $(^) $(LDFLAGS)

array_utils.o: array_utils.c
	$(CC) -fopenmp $(CFLAGS) -c $(^) $(LDFLAGS)

utils.o: utils.c
	$(CC) -fopenmp $(CFLAGS) -c $(^) $(LDFLAGS)

read_parameters.o: read_parameters.c
	$(CC) -fopenmp $(CFLAGS) -c $(^) $(LDFLAGS)

rng.o: rng.c
	$(CC) -fopenmp $(CFLAGS) -c $(^) $(LDFLAGS)

spike.o: spike.c
	$(CC) -fopenmp $(CFLAGS) -c $(^) $(LDFLAGS)

clean:
	rm -f *.o