export CFLAGS = -Wall -O2 -std=c89

.PHONY: all
all:
	make -C src scamp

.PHONY: clean
clean:
	make -C src clean

.PHONY: check
check: all
	./src/scamp tests/data1.fits.cat
