# Makefile
#
# s.abiko [2007.5]
CC = g++
#AR = ar

CFLAG = -O3 -Wall

#INLIBS = -I/usr/local/include
LDLIBS = -lm -L../matrix -lmatrix -lgsl -lgslcblas -lstdc++ 
#LDLIBS = 

SRC = spn.cpp rot.cpp
OBJ = spn.o rot.o

TARGET = test
#TARGET = libspn.a

$(TARGET) : $(OBJ)
	$(CC) $(CFLAG) -o $@ $^ $(INLIBS) $(LDLIBS)

#$(TARGET) : $(OBJ)
#	$(AR) -r $(TARGET) $(OBJ)

spn.o:spn.cpp
	$(CC) $(CFLAG) -o spn.o -c spn.cpp $(INLIBS) $(LDLIBS)

rot.o:rot.cpp
	$(CC) $(CFLAG) -o rot.o -c rot.cpp $(INLIBS) $(LDLIBS)

# # .cpp.o : 
# # 	$(CC) $(CFLAG) -o -$@ -c $< $(INLIBS) $(LDLIBS)

# # .c.o : 
# # 	$(CC) $(CFLAG) -o -$@ -c $<

# # spn.o: spn.h

clean :
	-@rm *.o $(TARGET)

# --- EOF ---
