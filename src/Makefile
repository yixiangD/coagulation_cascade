# cheetah = ORNL IBM machine, mpCC, native MPI, FFTW

SHELL = /bin/sh
.SUFFIXES: .cpp .u
.IGNORE:


# Definitions
                                                                                                             
ROOT =  lmp
EXE =   $(ROOT)_stam
SRC =   $(wildcard *.cpp)
INC =   $(wildcard *.h)
OBJ =   $(SRC:.cpp=.o)



# System-specific settings

CC =            mpicxx
CCFLAGS =   -O3 -g -DFFT_NONE -DGZIP -std=c++11
DEPFLAGS =	-M
LINK =	        mpicxx
LINKFLAGS =	-O3 -g
LIBDIR = 
USRLIB = 
SYSLIB = #-L$(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -liomp5 -lmkl_core -lmkl_intel_thread -lpthread -lm
SIZE =		size

# Link rule

$(EXE):	$(OBJ)
	$(LINK) $(LINKFLAGS) -o $(EXE) $(OBJ) $(USRLIB) $(SYSLIB) 
	$(SIZE) $(EXE)

# --------- old section -------------

# Compilation rules

#.cpp.o:
#	$(CC) $(CCFLAGS) -c $<

# Individual dependencies

#$(OBJ):	     $(INC)

# --------- new section -------------

# Compilation rules

%.o:%.cpp
	$(CC) $(CCFLAGS) -c $<

#%.u:%.cpp
#	$(CC) $(CCFLAGS) $(DEPFLAGS) -c $<

# Individual dependencies

#DEPENDS = $(OBJ:.o=.u)
#include $(DEPENDS)

clean:
	rm -f *.o $(EXE)
