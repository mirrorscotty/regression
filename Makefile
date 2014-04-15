CC=gcc
CFLAGS=-Imatrix -Imaterial-data/pasta -I. -ggdb
LDFLAGS=-lm
VPATH=matrix material-data/pasta programs
.PHONY:

all: kF gab fitdiff

matrix.a:
	$(MAKE) -C matrix
	cp matrix/matrix.a .

material-data.a:
	$(MAKE) -C material-data
	cp material-data/material-data.a .

kf.o:

gab.o:

fitdiff.o: diffusivity.h isotherms.h constants.h

fitnlm.o: fitnlm.h

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
	rm -rf *.o doc kF gab fitdiff
	$(MAKE) -C material-data clean
	$(MAKE) -C matrix clean
