#=======================================================================
# User definitions
#

MAIN=ddevp

#======================================================================
# Libraries
#
ifndef NETCDF
nonetcdf:
	@echo ""
	@echo "!!!NETCDF variable does not set!!!"
	@echo ""
	@echo "  You need to specify the NETCDF root directory as below:"
	@echo "    for bash: export NETCDF=<netcdf root directory>"
	@echo "    for  csh: setenv NETCDF <netcdf root directory>"
	@echo ""
endif

MPI = /opt/intel/impi/4.0.3.008
INCS = -I$(NETCDF)/include #-I$(MPI)/include64
LIBS = -L$(NETCDF)/lib   -lm -lz -lcurl #-L$(MPI)/lib64 -lmpi 

include make.defs
F90FILES = $(AFILES) $(OFILES) $(RFILES) $(SFILES) $(TFILES) $(MFILES) $(GFILES) $(PFILES)



#=======================================================================
# Standard definitions
#

OFILESF90 = $(F90FILES:.f90=.o)
OFILESC = $(CFILES:.c=.o)


VPATH = obj
.SUFFIXES:
.SUFFIXES: .out .o .F90 .f90 .f .c .inc .h


#=======================================================================
# Flags 
#
FLAGS = $(INCS) 
CFLAGS = -O2

#=======================================================================
# Compiler Option
#
include make.compiler


#=======================================================================
# Targets and dependencies
#
default: opt
all: $(MAIN) 
debug:
	@make all "FLAGS = -C -g -traceback $(GENERAL_FLAGS) -D_DEBUG_"
#prof:
#	@make all "LIBS = -p $(LIBS)"
opt:
	@make all "FLAGS = $(GENERAL_FLAGS)"

$(MAIN): mod obj $(OFILESF90) $(OFILESC)
	@echo "$(LD) $(OFILESF90) $(OFILESC) $(LIBS)"
	@cd obj; $(LD) $(LDOPTION) $(OFILESF90) $(OFILESC) $(LIBS) -o $(@); mv $(@) ../;

mod:
	mkdir mod
obj:
	mkdir obj
clean:
	rm -rf obj mod *.mod *.o $(MAIN)
run: 
	./$(MAIN) 
distclean:
#	@make clean;rm -rf make.defs *.f90 *.c *.h

#=======================================================================
# Compilation rules
#

.f90.o:
	$(F90) $(FLAGS) -c $*.f90 -o obj/$*.o
.c.o:
	$(CC) $(CFLAGS) -c $*.c -o obj/$*.o

