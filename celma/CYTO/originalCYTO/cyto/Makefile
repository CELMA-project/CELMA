####################################################################
#
#   Makefile
#
#################################################################### 
#
# 
ROOT = ..
HERE = $(PWD)
.SUFIXES = .o .c 

ifndef OSTYPE
OSTYPE=ubuntu
endif 

include $(ROOT)/$(OSTYPE).mk
LINCLUDE = $(HERE)/../include

# For real fast optimized programms switch back to original definitions
# of cflag makro

LINT = splint

#SPLINT_CMD = splint -f splint.opt $(ALL_DEFS) $(ALL_SPLINT_INC) -I$(<D)


F90FLAGS = $(F90OPTIM) -c

FFLAGS  = $(FOPTIM) -c

CCFLAGS = $(CFLAGS) -DCPP

LIBS   	= $(LOPTS)   $(FFT) $(HDF)  $(CLIBS) $(FLIBS)

LIBSOFFT= $(LOPTS)   $(HDF)  $(CLIBS) 

COWNF  	=   -o $(HERE)/$@.$(OSTYPE) 

COWNFOFFT  =  -o $(HERE)/$@.$(OSTYPE) $(LIBSOFFT) 

OBJDIR = .

BINDIR = $(HOME)/bin
#BINDIR = ../bin

# Keep up with changes of hidden dependencies (like eg. *.h files)



HEADERFILES = $(LINCLUDE)/utilities.h $(LINCLUDE)/hdf_helper.h $(LINCLUDE)/file_utils.h 

.KEEP_STATE:
.FAILED:


all:	
	make clean; make cyto; make install
clean:
	-rm  -f *.o *.$(OSTYPE) *.lo 

install :
	cp -v  cyto.$(OSTYPE) $(BINDIR)


CYTO =  cyto.o 
CYTO_TE  =  cyto_te.o 
CYTO_NEUTRAL = cyto_neutral.o Neutral_Plasma_Interaction.o cyto_initial_condition.o cyto_geometry.o


######################################################################

cyto: $(CYTO) $(HEADERFILES)  
	cd ../common  && $(MAKE) install
	$(PCC)   $(CYTO) $(COMMON) $(FLIBS)  $(COWNF) $(LIBS) $(MPI) $(FFTW) 

cyto_te: $(CYTO_TE) $(HEADERFILES)  
	cd ../common  && $(MAKE) install
	$(PCC)   $(CYTO_TE) $(COMMON) $(FLIBS)  $(COWNF) $(LIBS) $(MPI) $(FFTW) 

cyto_neutral: $(CYTO_NEUTRAL) $(HEADERFILES)  
	cd ../common  && $(MAKE) install
	$(PCC)   $(CYTO_NEUTRAL) $(COMMON) $(FLIBS)  $(COWNF) $(LIBS) $(MPI) $(FFTW) 

cyto_te.o: cyto.c
	$(PCC)  $(CFLAGS) -DTE -o cyto_te.o $<
######## General rules 


%.o : %.c
	$(PCC)  $(CFLAGS) $< -o $(@)

%.o:  %.f
	$(FC) -$(FFLAGS) $< -o $(@) 
.F.f:

	$(RM) $@
	$(CC) -P -D$(OSTYPE) $*.F
	$(MV) $*.i $@
