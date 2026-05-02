# makefile for genhas

CFLAGS  = -std=c99 -Wall -O2 -pedantic
LDLIBS  = -lm

all     : genhas
genhas  : genhas.o

clean   :
	rm -f genhas genhas.exe *.o

# regenerate the constant tables embedded in src/has.c
test    :
	./genhas > has_rs_tables.c

# self-check only: rebuild G and verify against the ICD Annex C D-matrix
check   :
	./genhas -c
