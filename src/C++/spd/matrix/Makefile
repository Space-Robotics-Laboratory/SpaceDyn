CC = gcc
CFLAGS = -O3 -Wall

LDLIBS = -lm -L./ -lmatrix -lgsl -lgslcblas -lstdc++ 
INLIBS = -I/usr/local/include

SRC = test.cpp
#OBJ = matrix.o vector.o
TARGET = test

$(TARGET) : $(SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(INLIBS) $(LDLIBS)

#.c.o :
#	$(CC) $(CFLAGS) -o $@ -c $<
clean :
	-@rm *.o $(TARGET)
