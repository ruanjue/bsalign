VERSION=1.2.1
RELEASE=20210520

CC  := gcc
BIN := /usr/local/bin
HEADCORE := .

ifeq (0, ${MAKELEVEL})
TIMESTAMP=$(shell date)
endif

ifeq (1, ${AVX2})
SIMD_FLAGS=-mavx2 -mavx
else
SIMD_FLAGS=-msse4.2
endif

ifeq (1, ${DEBUG})
CFLAGS=-g3 -W -Wall -Wno-unused-but-set-variable -O0 -DDEBUG=1 -DVERSION="$(VERSION)" -DRELEASE="$(RELEASE)" -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE -mpopcnt $(SIMD_FLAGS)
else
CFLAGS=-g3 -W -Wall -Wno-unused-but-set-variable -O4 -DVERSION="$(VERSION)" -DRELEASE="$(RELEASE)" -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE -mpopcnt $(SIMD_FLAGS)
endif

GLIBS=-lm -lrt -lpthread -lz
GENERIC_SRC=$(HEADCORE)/mem_share.h $(HEADCORE)/chararray.h $(HEADCORE)/sort.h $(HEADCORE)/list.h $(HEADCORE)/pgzf.h $(HEADCORE)/sort.h $(HEADCORE)/dna.h \
	$(HEADCORE)/thread.h $(HEADCORE)/filereader.h $(HEADCORE)/bitvec.h $(HEADCORE)/hashset.h

PROGS=bsalign

all: $(PROGS)

bsalign: $(GENERIC_SRC) bsalign.h bspoa.h main.c
	$(CC) $(CFLAGS) -o $@ main.c $(GLIBS)

clean:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out $(PROGS)

clear:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out

install: $(PROGS)
	mkdir -p $(BIN) && cp -fvu $(PROGS) $(BIN)
