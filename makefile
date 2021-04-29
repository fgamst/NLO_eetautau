# Makefile to eventually compile all of NLO eetautau.
# was previously compiled using the following commands:
#
# icc -c main.c
# icc main.o ../vegas/nvegas.o ../Madgraph/emep_taptam.o ../Madgraph/emep_taptama.o ../Madgraph/MGrWrapper.a -lifcore -limf -lgsl -lgslcblas
#
# using -limf being the math library called via -lm in gcc
# using -lgsl -lgslcblas to call gsl needed for the dilogarithm/spences' function. (in gcc?)
# (-l)ifcore is an intel fortran library needed to link fortran code containing 
# STOP statements. (Quote Stackoverflow: "C does not like the STOP statement")
# otherwise there will be errors of the type "for_...._stop not found".

#
# 'make depend' uses makedepend to automatically generate dependencies 
#               (dependencies are added to end of Makefile)
# 'make'        build executable file
# 'make clean'  removes all .o and executable files
#

# define the C compiler to use
CC = gcc

# define the C compiler to use
CPC = g++

# define the FORTRAN compiler to use
FC = gfortran

# define any C compile-time flags
# -g does what?
CFLAGS = 
# -Wall -g

# define any FORTRAN compile-time flags
# -Wall does what?
# -g does what?
FFLAGS = 

# define any directories containing header files other than /usr/include
INCLUDES = -I./lib/qcdloop/src -I./lib/HELAS

# define library paths in addition to /usr/lib
LFLAGS = -L./lib

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -lx 
#   option:
LIBS = -lm -ldl -lgsl -lgslcblas -lqcdloop -lquadmath -lhelas -lgfortran
#  spot to salvage temporary unused libraries


# define the FORTRAN90 source files
 SRCSF90 = Madgraph/MGWrapper.f90

# define the FORTRAN source files
 SRCSF = Madgraph/emep_taptama.f Madgraph/emep_taptam_gZ.f
# emep_taptam.f no need for this guy atm I checked my matrix element to be correct
 


# define the C++ source files
SRCSCC = main.cc 
# wrapper.cc 

# define the C source files
# SRCSC = vegas/nvegas.c

# define the object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
OBJS = $(SRCSF90:.f90=.o) $(SRCSF:.f=.o)
# $(SRCSCC:.cc=.o)
#Madgraph/MGWrapper.o emep_taptama.o emep_taptam_gZ.o 



# define the executable file 
MAIN = bin/eetautau

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

THIS = .
#VPATH = $(THIS):$(THIS)/vegas

.PHONY: depend clean

all:    $(MAIN)
	@echo Successfully compiled eetautau

$(MAIN): $(OBJS) 
	$(CPC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(SRCSCC) $(OBJS) $(LFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)

%.o: %.f90
	$(CPC) $(FFLAGS) -c $< -o $@
.f.o:
	$(CPC) $(FFLAGS) $(INCLUDES) -c $< -o $@

%.o: %.cc
	$(CPC) $(CFLAGS) $(INCLUDES) -c $< -o $@
.c.o:
	$(CPC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJS) $(MAIN)
	#*.o *~ 

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
