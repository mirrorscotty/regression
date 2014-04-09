CC=gcc
CFLAGS=-Imatrix -lm -ggdb
VPATH=matrix
.PHONY: matrix

all: kF gab

matrix:
	$(MAKE) -C matrix

kf.o:

gab.o:

fitnlm.o: fitnlm.h

regress.o: regress.h

gab: fitnlm.o gab.o matrix
	$(CC) $(CFLAGS) -o gab gab.o fitnlm.o matrix/*.o

kF: fitnlm.o kf.o matrix
	$(CC) $(CFLAGS) -o kF kf.o fitnlm.o matrix/*.o 

clean:
	rm *.o matrix/*.o kF gab

