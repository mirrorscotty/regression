CC=gcc
CFLAGS=-Imatrix -Imaterial-data/pasta -I. -ggdb
LDFLAGS=-lm
VPATH=matrix material-data/pasta programs
.PHONY:

all: kF gab fitdiff

matrix.a:
	$(MAKE) -C matrix matrix.a
	cp matrix/matrix.a .

material-data.a: matrix.a
	cp matrix.a material-data
	$(MAKE) -C material-data material-data.a
	cp material-data/material-data.a .

kf.o:

gab.o:

fitdiff.o: diffusivity.h isotherms.h constants.h

fitnlm.o: regress.h

regress.o: regress.h

gab: fitnlm.o gab.o matrix.a
	$(CC) $(CFLAGS) -o $@ $? $(LDFLAGS)

kF: fitnlm.o kf.o matrix.a
	$(CC) $(CFLAGS) -o $@ $? $(LDFLAGS)

fitdiff: fitdiff.o regress.o matrix.a material-data.a
	$(CC) $(CFLAGS) -o $@ $? $(LDFLAGS)

doc: Doxyfile
	doxygen Doxyfile

clean:
	rm -rf *.o *.a doc kF gab fitdiff
	$(MAKE) -C material-data clean
	$(MAKE) -C matrix clean
