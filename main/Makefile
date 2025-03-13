CC=gcc
CFLAGS=-Wall -ggdb -ffast-math 
LDLIBS=-lgsl -lgslcblas -lm
#VPATH=../Library
COBJS=exec.o functions.o odes.o vacua.o

all: exec

functions.o: ../Library/functions.c ../Library/functions.h 
	$(CC) $(CFLAGSSAFE) -c $<

odes.o: ../Library/odes.c ../Library/odes.h 
	$(CC) $(CFLAGSSAFE) -c $<

vacua.o: ../Library/vacua.c ../Library/vacua.h 
	$(CC) $(CFLAGSSAFE) -c $<

exec: $(COBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

exec.o: exec.c
	$(CC) $(CFLAGS) -c $<

clean: 
	rm -f  exec $(COBJS) *.rpo

final:
	rm -f $(COBJS) *.rpo
