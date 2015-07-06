CC=gcc
CFLAGS=-Imatrix -Imaterial-data -I. -ggdb -Wall
LDFLAGS=-lm
VPATH=matrix material-data material-data/pasta programs programs/kF programs/modulus
SRC=$(wildcard *.c) \
	$(wildcard programs/*.c) \
	$(wildcard programs/kF/*.c) \
	$(wildcard programs/modulus/*.c)

all: kF gab fitdiff modulus

# Make stuff from other projects using their makefile
matrix/matrix.a:
	$(MAKE) -C matrix matrix.a

material-data/material-data.a: 
	$(MAKE) -C material-data material-data.a

kF: programs/kF/calc.o programs/kF/crank.o programs/kF/io.o programs/kF/Xe.o programs/kF/L.o programs/kF/kFmain.o fitnlm.o regress.o programs/kF/De.o programs/kF/flux.o matrix/matrix.a material-data/material-data.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# GAB program
gab: fitnlm.o programs/gab.o matrix.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# fitdiff program
fitdiff: programs/fitdiff.o regress.o matrix/matrix.a material-data/material-data.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# modulus program
modulus: fitnlm.o programs/modulus/modulus.o programs/modulus/stress-strain.o matrix/matrix.a material-data/material-data.a 
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# fitburgers program
fitburgers: programs/fitburgers.o fitnlmM.o matrix.a

fitcreep: programs/fitcreep.o regress.o matrix/matrix.a
nlin-fitcreep: programs/nlin-fitcreep.o fitnlm.o material-data/material-data.a matrix/matrix.a
nlin-fitcreepv2: programs/nlin-fitcreepv2.o fitnlmP.o material-data/material-data.a matrix/matrix.a
creep-table: programs/creep-table.o fitnlmP.o material-data/material-data.a matrix/matrix.a

doc: Doxyfile
	doxygen Doxyfile

clean:
	rm -rf doc kF gab fitdiff modulus
	rm -rf $(SRC:.c=.o)
	rm -rf $(SRC:.c=.d)
	rm -rf *.a
	$(MAKE) -C material-data clean
	$(MAKE) -C matrix clean

%.o: %.c
	$(CC) -c $(CFLAGS) $*.c -o $*.o
	$(CC) -MM $(CFLAGS) $*.c > $*.d
	@mv -f $*.d $*.d.tmp
	@sed -e 's|.*:|$*.o:|' < $*.d.tmp > $*.d
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
          sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.d.tmp

-include $(SRC:.c=.d)

