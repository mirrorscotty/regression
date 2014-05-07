CC=gcc
CFLAGS=-Imatrix -Imaterial-data/choi-okos -Imaterial-data/pasta -I. -ggdb -Wall
LDFLAGS=-lm
VPATH=matrix material-data material-data/pasta programs programs/kF

all: kF gab fitdiff

matrix.a:
	$(MAKE) -C matrix matrix.a
	cp matrix/matrix.a .

material-data.a: matrix.a
	cp matrix/matrix.a material-data
	$(MAKE) -C material-data material-data.a
	cp material-data/material-data.a .

calc.o: kf.h

crank.o: kf.h

io.o: kf.h

Xe.o: kf.h

De.o: kf.h

flux.o: kf.h

kFmain.o: kf.h

gab.o:

fitdiff.o: diffusivity.h isotherms.h constants.h

fitnlm.o: regress.h

regress.o: regress.h

gab: fitnlm.o gab.o matrix.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

kF: calc.o crank.o io.o Xe.o kFmain.o fitnlm.o regress.o De.o flux.o matrix.a material-data.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

fitdiff: fitdiff.o regress.o matrix.a material-data.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

doc: Doxyfile
	doxygen Doxyfile

clean:
	rm -rf *.o *.a doc kF gab fitdiff
	$(MAKE) -C material-data clean
	$(MAKE) -C matrix clean
