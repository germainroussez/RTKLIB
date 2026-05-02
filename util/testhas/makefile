# makefile for testhas
#
# Each test_has_*.c is a self-contained standalone validation driver
# linking against src/has.c. Built without -DTRACE so the trace() macros
# expand to no-ops (see rtklib.h), keeping the test free of any link
# dependency on src/trace.c.

SRC    = ../../src
CFLAGS = -std=c99 -Wall -O2 -pedantic -I$(SRC) \
         -DENAGLO -DENAQZS -DENACMP -DENAGAL -DENAIRN \
         -DNFREQ=3 -DNEXOBS=3 -DSVR_REUSEADDR
LDLIBS = -lm

all         : test_has_rs

test_has_rs : test_has_rs.o has.o

has.o       : $(SRC)/rtklib.h $(SRC)/has.c
	$(CC) -c $(CFLAGS) $(SRC)/has.c

test        :
	./test_has_rs

clean       :
	rm -f test_has_rs test_has_rs.exe *.o
