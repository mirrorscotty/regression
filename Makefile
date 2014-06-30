CC=gcc
CFLAGS=-Imatrix -Imaterial-data/choi-okos -Imaterial-data/pasta -I. -ggdb -Wall
LDFLAGS=-lm
VPATH=matrix material-data material-data/pasta programs programs/kF programs/modulus

all: kF gab fitdiff modulus

# Make stuff from other projects using their makefile
matrix.a:
	$(MAKE) -C matrix matrix.a
	cp matrix/matrix.a .

material-data.a: matrix.a
	cp matrix/matrix.a material-data
	$(MAKE) -C material-data material-data.a
	cp material-data/material-data.a .

# FOr the kF program
calc.o: kf.h
crank.o: kf.h
io.o: kf.h
Xe.o: kf.h
De.o: kf.h
L.o: kf.h
flux.o: kf.h
kFmain.o: kf.h
kF: calc.o crank.o io.o Xe.o L.o kFmain.o fitnlm.o regress.o De.o flux.o matrix.a material-data.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# GAB program
gab.o:
gab: fitnlm.o gab.o matrix.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# fitdiff program
fitdiff.o: diffusivity.h isotherms.h constants.h
fitdiff: fitdiff.o regress.o matrix.a material-data.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# modulus program
modulus.o: stress-strain.h pasta.h matrix.h
stress-strain.o: stress-strain.h pasta.h matrix.h
modulus: fitnlm.o modulus.o stress-strain.o matrix.a material-data.a 
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# regression files in the main directory
fitnlm.o: regress.h
regress.o: regress.h

doc: Doxyfile
	doxygen Doxyfile

clean:
	rm -rf *.o *.a doc kF gab fitdiff modulus
	$(MAKE) -C material-data clean
	$(MAKE) -C matrix clean
