####################################################################
#
#   Makefile
#
####################################################################
#
#
ifndef OSTYPE
OSTYPE=ubuntu
endif

ROOT = $(PWD)/..
include $(ROOT)/$(OSTYPE).mk


OUT = ./libcommon.a

SRC =  file_utils.c  hdf_helper.c fft_wrap.c laplace_solver.c par_helper.c util_alloc.c util_inifiles.c utilities.c

OBJ = $(SRC:.c=.o)

HEADERFILES = ../include/*.h

.SUFFIXES=.c

.c.o:
	$(CC) $(CFLAGS) -c $< -o$@	

$(OUT):$(OBJ)
	ar rcs $(OUT) $(OBJ)

depend: dep

dep: 
	makedepend --$(CFLAGS) --$(HEADERFILES) $(SRC)

clean:
	rm -f $(OBJ) $(OUT) *.bak 

install:
	make; cp -v $(OUT) $(HOME)/lib 
