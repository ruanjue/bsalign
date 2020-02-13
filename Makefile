VERSION=1.0.0
RELEASE=20200214

CC  := gcc
BIN := /usr/local/bin

ifeq (0, ${MAKELEVEL})
TIMESTAMP=$(shell date)
endif

ifeq (1, ${DEBUG})
CFLAGS=-g3 -W -Wall -Wno-unused-but-set-variable -O0 -DDEBUG=1 -DVERSION="$(VERSION)" -DRELEASE="$(RELEASE)" -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE -mpopcnt -msse4.2
else
CFLAGS=-g3 -W -Wall -Wno-unused-but-set-variable -O4 -DVERSION="$(VERSION)" -DRELEASE="$(RELEASE)" -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE -mpopcnt -msse4.2
endif

GLIBS=-lm -lrt -lpthread -lz
GENERIC_SRC=mem_share.h chararray.h sort.h list.h pgzf.h sort.h dna.h thread.h filereader.h bitvec.h hashset.h

PROGS=bsalign

all: $(PROGS)

bsalign: $(GENERIC_SRC) bsalign.h main.c
	$(CC) $(CFLAGS) -o $@ main.c $(GLIBS)

clean:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out $(PROGS)

clear:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out

install: $(PROGS)
	mkdir -p $(BIN) && cp -fvu $(PROGS) $(BIN)
