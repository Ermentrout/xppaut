CFLAGS=-g -O -I/usr/local/include -I/usr/local/include/sbml
LDFLAGS=-L/usr/local/lib
SRC=sbml2xpp.c
OBJ=sbml2xpp.o
LIBS= -lsbml -lexpat  
CC= gcc
s2x:	$(OBJ)
	$(CC) -o s2x $(OBJ) $(CFLAGS) $(LDFLAGS) $(LIBS)


