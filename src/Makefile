SOURCES = main.c kdtree.c catalog.c
HEADERS = scamp.h kdtree.h catalog.h
OBJECTS = $(SOURCES:.c=.o)

CFLAGS ?= -Wall
CFLAGS += $(shell pkg-config --cflags wcslib)
CFLAGS += $(shell pkg-config --cflags cfitsio)
LDFLAGS += $(shell pkg-config --libs cfitsio)
LDFLAGS += $(shell pkg-config --libs wcslib)

scamp: $(OBJECTS) $(HEADERS)
	$(CC) -o scamp $(OBJECTS) $(CFLAGS) $(LDFLAGS)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	$(RM) scamp $(OBJECTS)
