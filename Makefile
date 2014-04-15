CC=gcc
CFLAGS=-Imatrix -Imaterial-data/pasta -lm -ggdb
VPATH=matrix material-data/pasta
.PHONY: matrix pasta

all: kF gab fitdiff

matrix:
	$(MAKE) -C matrix

pasta:
	$(MAKE) -C material-data/pasta

kf.o:

gab.o:

fitdiff.o: diffusivity.h isotherms.h constants.h

fitnlm.o: fitnlm.h

regress.o: regress.h

gab: fitnlm.o gab.o matrix
	$(CC) $(CFLAGS) -o gab gab.o fitnlm.o matrix/*.o

kF: fitnlm.o kf.o matrix
	$(CC) $(CFLAGS) -o kF kf.o fitnlm.o matrix/*.o 

fitdiff: fitdiff.o regress.o matrix pasta
	$(CC) $(CFLAGS) -o fitdiff fitdiff.o regress.o matrix/*.o material-data/pasta/binding.o material-data/pasta/diffusivity.o material-data/pasta/isotherms.o material-data/pasta/choi-okos.o material-data/pasta/composition.o material-data/pasta/gas.o material-data/pasta/fluid.o

doc: Doxyfile
	doxygen Doxyfile

clean:
	rm -rf *.o matrix/*.o doc kF gab fitdiff

