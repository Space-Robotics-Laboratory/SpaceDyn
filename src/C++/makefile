CC = g++

CFLAGS = -static -Wall -O2 -fno-exceptions -fno-rtti 

#for CYGWIN in Satoko's laptop or in Tohoku Univ.
LIB = /usr/local/lib
LDLIBS = -lm -L ./spd/lib -lspd -L ${LIB} -lgsl -lgslcblas -lstdc++
INLIBS = -I/usr/local/include -I./spd/include

SRC = main_test.cpp

OBJ2 = main_test.o \

TARGET = run

$(TARGET) : $(OBJ2)
	$(CC) $(CFLAGS) -o $@ $^ $(INLIBS) $(LDLIBS)

.cpp.o :
	$(CC) $(CFLAGS) $(INLIBS) -o $@ -c $<

clean :
	-@rm -f $(OBJ2) $(TARGET)

