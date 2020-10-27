CC = g++
CFLAGS = -Wall -O2
AR = ar

## IKUTA's Laptop set
INCLUDE = -I/usr/local/include -I../include

## DLR PC for Satoko set
#INCLUDE = -I/home/abiko/project/DYNLIB/gsl/include -I../include

SRC_DIR = src/
MAT_DIR = matrix/
SPN_DIR = spn/
LIB_DIR = lib/

SRC =$(SRC_DIR)model.cpp\
	$(SRC_DIR)model_init.cpp\
	$(SRC_DIR)calc_SP.cpp\
	$(SRC_DIR)calc_SPN.cpp\
	$(SRC_DIR)calc_hh.cpp\
	$(SRC_DIR)calc_Rg.cpp\
	$(SRC_DIR)calc_Je.cpp\
	$(SRC_DIR)calc_Jb.cpp\
	$(SRC_DIR)calc_GJe.cpp\
	$(SRC_DIR)calc_Lg.cpp\
	$(SRC_DIR)calc_JJ.cpp\
	$(SRC_DIR)calc_SS.cpp\
	$(SRC_DIR)calc_Xup_I.cpp\
	$(SRC_DIR)j_num.cpp\
	$(SRC_DIR)f_kin_e.cpp\
	$(SRC_DIR)f_kin_j.cpp\
	$(SRC_DIR)i_dyn.cpp\
	$(SRC_DIR)i_dyn_fix.cpp\
	$(SRC_DIR)f_dyn.cpp\
	$(SRC_DIR)f_dyn_fix.cpp\
	$(SRC_DIR)inner_force.cpp\
	$(MAT_DIR)matrix.cpp\
	$(MAT_DIR)vector.cpp\
	$(SPN_DIR)spn.cpp\
	$(SPN_DIR)rot.cpp\
	$(SRC_DIR)rref.cpp\
	$(SRC_DIR)aw\
	$(SRC_DIR)w2dQtn.cpp\
	$(SRC_DIR)f_dyn_rk.cpp\
	$(SRC_DIR)model_test.cpp\
	$(SRC_DIR)calc_GJb.cpp\
	$(SRC_DIR)calc_GJe_r.cpp

OBJ =$(SRC_DIR)model.o\
	$(SRC_DIR)model_init.o\
	$(SRC_DIR)calc_SP.o\
	$(SRC_DIR)calc_SPN.o\
	$(SRC_DIR)calc_hh.o\
	$(SRC_DIR)calc_Rg.o\
	$(SRC_DIR)calc_Je.o\
	$(SRC_DIR)calc_Jb.o\
	$(SRC_DIR)calc_GJe.o\
	$(SRC_DIR)calc_Lg.o\
	$(SRC_DIR)calc_JJ.o\
	$(SRC_DIR)calc_SS.o\
	$(SRC_DIR)calc_Xup_I.o\
	$(SRC_DIR)j_num.o\
	$(SRC_DIR)f_kin_e.o\
	$(SRC_DIR)f_kin_j.o\
	$(SRC_DIR)i_dyn.o\
	$(SRC_DIR)i_dyn_fix.o\
	$(SRC_DIR)f_dyn.o\
	$(SRC_DIR)f_dyn_fix.o\
	$(SRC_DIR)inner_force.o\
	$(MAT_DIR)matrix.o\
	$(MAT_DIR)vector.o\
	$(SPN_DIR)spn.o\
	$(SPN_DIR)rot.o\
	$(SRC_DIR)rref.o\
	$(SRC_DIR)aw.o\
	$(SRC_DIR)w2dQtn.o\
	$(SRC_DIR)f_dyn_rk.o\
	$(SRC_DIR)model_test.o\
	$(SRC_DIR)calc_GJb.o\
	$(SRC_DIR)calc_GJe_r.o

TARGET = libspd.a

$(LIB_DIR)$(TARGET) : $(OBJ)
	$(AR) -r $(LIB_DIR)$(TARGET) $(OBJ)

.cpp.o :
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $<

.c.o :
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $<

#clean :
#	-@rm -f $(OBJ) $(LIB_DIR)$(TARGET)
clean :
	-@rm -f $(OBJ)
