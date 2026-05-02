# makefile for testhas
#
# Each test_has_*.c is a self-contained standalone validation driver
# linking against src/has.c. Built without -DTRACE so the trace() macros
# expand to no-ops (see rtklib.h), keeping the test free of any link
# dependency on src/trace.c. The page-parser and MT1 sections of has.c
# call getbitu() / timediff(), so rtkcmn.c is linked in alongside.

SRC    = ../../src
CFLAGS = -std=c99 -Wall -O2 -pedantic -I$(SRC) \
         -DENAGLO -DENAQZS -DENACMP -DENAGAL -DENAIRN \
         -DNFREQ=3 -DNEXOBS=3 -DSVR_REUSEADDR
LDLIBS = -lm

all             : test_has_rs test_has_decode test_has_mt1

test_has_rs     : test_has_rs.o     has.o rtkcmn.o
test_has_decode : test_has_decode.o has.o rtkcmn.o
test_has_mt1    : test_has_mt1.o    has.o rtkcmn.o

has.o           : $(SRC)/rtklib.h $(SRC)/has.c
	$(CC) -c $(CFLAGS) $(SRC)/has.c
rtkcmn.o        : $(SRC)/rtklib.h $(SRC)/rtkcmn.c
	$(CC) -c $(CFLAGS) $(SRC)/rtkcmn.c

test            :
	./test_has_rs && ./test_has_decode && ./test_has_mt1

clean           :
	rm -f test_has_rs     test_has_rs.exe \
	      test_has_decode test_has_decode.exe \
	      test_has_mt1    test_has_mt1.exe *.o
