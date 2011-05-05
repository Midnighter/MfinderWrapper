# 
CC=gcc
# 
CFLAGS=-c -Wall -DUNIX=1 -O3 -DNDEBUG

LDFLAGS= -lm -O3 -DNDEBUG

FILES=main clustering  grassberger  hash  list mat  metropolis  motif_ids  output  permutation  prob  random  results  role  stubs  switches
SOURCES=$(FILES:=.c)
OBJECTS=$(FILES:=.o)
EXECUTABLES=mfinder

all: $(SOURCES) $(EXECUTABLES)
	
$(EXECUTABLES): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $(EXECUTABLES)

$(FILES):
	$(CC) $(CFLAGS) $@.c -o $@.o

clean:
	rm -f *o *~ mfinder
