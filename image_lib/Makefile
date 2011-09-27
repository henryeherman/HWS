CC=gcc
CFLAGS=-c
LDFLAGS=
SOURCES=image_lib.c
LIBS=-lm
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLES=cal hw1


all: $(SOURCE) $(EXECUTABLES)

cal: intro.c $(OBJECTS)
		$(CC) $(LDFLAGS) $(OBJECTS) $(LIBS) $< -o $@

hw1: hw1.c $(OBJECTS)
		$(CC) $(LDFLAGS) $(OBJECTS) $(LIBS) $< -o $@
		
.c.o:
		$(CC) $(CFLAGS) $< -o $@

clean:
		-rm *.o *.exe
		