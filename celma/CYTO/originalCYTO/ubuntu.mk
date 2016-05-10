DBF = -DDATABASEPATH='"'$(HOME)'/database"'

#ICC = 1
GCC = 1
#DEBUG =1 
GNU = 1

#Includes to be found here

ifdef GNU
MPI_L_PATH =/usr/lib/openmpi/lib
MPI_I_PATH =/usr/lib/openmpi/include
endif

ifdef ICC
MPI_L_PATH =  /opt/intel//impi/4.1.1.036/intel64/lib
MPI_I_PATH =  /opt/intel//impi/4.1.1.036/intel64/include
endif


MKL_PATH =


HDF_I_PATH =/usr/include/hdf
HDF_L_PATH =/home/vona/lib

FFTW_I_PATH =/usr/include
FFTW_L_PATH =/usr/lib64	
#FFTW_I_PATH =$(HOME)/include
#FFTW_L_PATH =$(HOME)/lib

COMMON_PATH =$(HOME)



FFTLIB=FFTW 
# or make it:  Intel_MPL


INCLUDE = -I$(ROOT)/include -I$(HERE)  -I$(ROOT)/lib/include  -I/usr/include -I/usr/local/include -I/usr/local/include/plplot  -I$(MPI_I_PATH)  -I$(MKL_PATH)/include  -I$(HDF_I_PATH) -I$(FFTW_I_PATH) -I$(COMMON_PATH)/include


#the gcc compiler

ifdef GCC
CC = gcc
CCFLAG =  -std=gnu99 $(DBF) -c  -Wstrict-aliasing=2  -Wall -pedantic -W  -Wstrict-prototypes -Wshadow -Wpointer-arith  -Wnested-externs -Wcast-align -Wwrite-strings   -D$(OSTYPE) -D$(FFTLIB) $(INCLUDE)
CDEBUG = -g
COPTIM = -O4
CPROF = -O2 -p


PCC    = mpicc 
PCCFLAG =   $(CCFLAG) 
PCDEBUG = -g
PCOPTIM =  -O4
PCPROF  = -G
endif

ifdef ICC
#the intel compiler 
CC = icc
CCFLAG = $(DBF) -D$(OSTYPE) -Dlinux -D$(FFTLIB) -Wstrict-prototypes $(INCLUDE)  -c 
CCOPTIM = -O2   
CCPROF =  -qp
CCDEBUG = -O0 -w2 -g


PCC    =  mpicc  
PCCFLAG =  $(CCFLAG) 
PCOPTIM =  -O2  
PCDEBUG = -O0 -w2 -g
PCPROF  = -qp 
endif

LLOPTS =  -L/soft/intel/composerxe-2011.3.174/compiler/lib/intel64/ -lifcoremt -Bstatic -lirc -Bdynamic
LLOPTS=  -O4  -L$(HDF_L_PATH) -L$(FFTW_L_PATH) -L$(COMMON_PATH)/lib  -L/usr/lib64 -L/usr/X11R6/lib -L/usr/local/lib -L$(MPI_L_PATH)  -Bstatic  -Bdynamic  

#LLOPTS=  -std=gnu99  -O2  -L$(HDF_L_PATH) -L$(FFTW_L_PATH) -L$(COMMON_PATH)/lib  -L/usr/lib64 -L/usr/X11R6/lib -L/usr/local/lib -L$(MPI_L_PATH)
LPROF = -p $(LOPTS)


ifdef DEBUG
CFLAGS =   $(CDEBUG) $(CCFLAG)
PCFLAGS =  $(PCDEBUG) $(PCCFLAG) 
LOPTS  = $(CDEBUG) $(LLOPTS)
else
CFLAGS =   $(CCOPTIM) $(CCFLAG)
PCFLAGS =  $(PCOPTIM) $(PCCFLAG) 
LOPTS  = $(CCOPTIM) $(LLOPTS)
endif 


FC     = ifort
FCFLAG = -c
FOPTIM = -O2 -axP  -tpp7 -w  -extend_source 
FDEBUG = -g
FPROF  = -fast -O4 -p
FDOUBLE= 
LF    = ifort


# local libraries / obj's
XLIBS     =  -lXaw -lXext -lXt -lXmu -lX11
CLIBS = -lm
FLIBS =
FFT =
MKL = 
FFTW=  -lfftw3
HDF = -lmfhdf -ldf -ljpeg -lz
PLPLOT = `pkg-config --cflags --libs plplotd`
COMMON = $(COMMON_PATH)/lib/libcommon.a
TCLTK =  -ltk8.5  -ltcl8.5  

MV = mv
RM = rm -f
CP = cp
