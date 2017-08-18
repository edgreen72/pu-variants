CC=gcc
#CFLAGS=-O2
CFLAGS=-g

pileup.o: pileup.c pileup.h
	echo "Making pileup.o..."
	$(CC) $(CFLAGS) -c pileup.c

pu-variants: pu-variants.c pileup.o
	echo "Making pu-variants..."
	$(CC) $(CFLAGS) -o pu-variants pileup.o pu-variants.c

clean:
	echo "Removing targets..."
	rm pileup.o pu-variants

