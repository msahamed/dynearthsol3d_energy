# -*- Makefile -*-
#
# Makefile for DynEarthSol3D

# Author: Eh Tan <tan2@earth.sinica.edu.tw>
# Energy module Author: sabber ahamed<sabbers@gmail.com>

## Execute "make" if making production run. Or "make opt=0 openmp=0" for debugging run.
##
## ndims = 3: 3D code; 2: 2D code
## opt = 1 ~ 3: optimized build; others: debugging build
## openmp = 1: enable OpenMP
## useadapt = 1: use libadaptivity for mesh optimization during remeshing
## energy = 1: use full energy balance equation and its associated functions like :
## 			update_thermal_energy

ndims=2
opt=2
openmp=0
useadapt=0
energy=1

## Select C++ compiler
ifeq ($(useadapt), 1)
	CXX = mpic++-mpich-gcc49 # g++-mp-4.7
else
	CXX=g++
	CXX_BACKEND = ${CXX}
endif

## Boost location and library name
BOOST_ROOT_DIR = ${HOME}/Library/Boost


########################################################################
## Select compiler and linker flags
## (Usually you won't need to modify anything below)
########################################################################

OSNAME := $(shell uname -s)

BOOST_LDFLAGS = -lboost_program_options
ifdef BOOST_ROOT_DIR
	# check existence of stage/ directory
	has_stage_dir = $(wildcard $(BOOST_ROOT_DIR)/stage)
	ifeq (, $(has_stage_dir))
		# no stage dir, BOOST_ROOT_DIR is the installation directory
		BOOST_CXXFLAGS = -I$(BOOST_ROOT_DIR)/include
		BOOST_LIB_DIR = $(BOOST_ROOT_DIR)/lib
	else
		# with stage dir, BOOST_ROOT_DIR is the build directory
		BOOST_CXXFLAGS = -I$(BOOST_ROOT_DIR)
		BOOST_LIB_DIR = $(BOOST_ROOT_DIR)/stage/lib
	endif
	BOOST_LDFLAGS += -L$(BOOST_LIB_DIR)
	ifneq ($(OSNAME), Darwin)  # Apple's ld doesn't support -rpath
		BOOST_LDFLAGS += -Wl,-rpath=$(BOOST_LIB_DIR)
	endif
endif

ifneq (, $(findstring g++, $(CXX_BACKEND))) # if using any version of g++
	CXXFLAGS = -g -std=c++0x
	LDFLAGS = -lm

	ifeq ($(opt), 1)
		CXXFLAGS += -O1
	else ifeq ($(opt), 2)
		CXXFLAGS += -O2
	else ifeq ($(opt), 3) # experimental, use at your own risk :)
		CXXFLAGS += -march=native -O3 -ffast-math -funroll-loops
	else # debugging flags
		CXXFLAGS += -O0 -Wall -Wno-unused-variable -Wno-unused-function -Wno-unknown-pragmas -fbounds-check -ftrapv
	endif

	ifeq ($(openmp), 1)
		CXXFLAGS += -fopenmp -DUSE_OMP
		LDFLAGS += -fopenmp
	endif

	ifeq ($(useadapt), 1)
		CXXFLAGS += -I$(VTK_INCLUDE)
	endif

else ifneq (, $(findstring icpc, $(CXX_BACKEND))) # if using intel compiler, tested with v14
	CXXFLAGS = -g -std=c++0x
	LDFLAGS = -lm

	ifeq ($(opt), 1)
		CXXFLAGS += -O1
	else ifeq ($(opt), 2)
		CXXFLAGS += -O2
	else ifeq ($(opt), 3) # experimental, use at your own risk :)
		CXXFLAGS += -fast -fast-transcendentals -fp-model fast=2
	else # debugging flags
		CXXFLAGS += -O0 -check=uninit -check-pointers=rw -check-pointers-dangling=all -fp-trap-all=all
	endif

	ifeq ($(openmp), 1)
		CXXFLAGS += -fopenmp -DUSE_OMP
		LDFLAGS += -fopenmp
	endif

	ifeq ($(useadapt), 1)
		CXXFLAGS += -I$(VTK_INCLUDE)
	endif

else
# the only way to display the error message in Makefile ...
all:
	@echo "Unknown compiler, check the definition of 'CXX' in the Makefile."
	@false
endif

## Is this a mercurial repository?
HAS_HG := $(shell hg log -r tip --template '{node}' 2>/dev/null)

##
ifeq ($(useadapt), 1)
	REMESHING_FILE = remeshing_adapt.cxx
else
	REMESHING_FILE = remeshing.cxx
endif

## Energy balance files
ifeq ($(energy), 1)
	ENERGY_FILE = energy/energy.cxx
endif

SRCS =	\
	barycentric-fn.cxx \
	brc-interpolation.cxx \
	bc.cxx \
	binaryio.cxx \
	dynearthsol.cxx \
	$(ENERGY_FILE)\
	fields.cxx \
	geometry.cxx \
	ic.cxx \
    ic-read-temp.cxx \
	input.cxx \
	matprops.cxx \
	mesh.cxx \
	nn-interpolation.cxx \
	output.cxx \
	phasechanges.cxx \
	$(REMESHING_FILE) \
	rheology.cxx \
	markerset.cxx

INCS =	\
	array2d.hpp \
	barycentric-fn.hpp \
	binaryio.hpp \
	constants.hpp \
	parameters.hpp \
	matprops.hpp \
	sortindex.hpp \
	utils.hpp \
	mesh.hpp \
	markerset.hpp \
	output.hpp

OBJS = $(SRCS:.cxx=.$(ndims)d.o)

EXE = dynearthsol$(ndims)d


## Libraries

TET_SRCS = tetgen/predicates.cxx tetgen/tetgen.cxx
TET_INCS = tetgen/tetgen.h
TET_OBJS = $(TET_SRCS:.cxx=.o)

TRI_SRCS = triangle/triangle.c
TRI_INCS = triangle/triangle.h
TRI_OBJS = $(TRI_SRCS:.c=.o)

M_SRCS = $(TRI_SRCS)
M_INCS = $(TRI_INCS)
M_OBJS = $(TRI_OBJS)

ifeq ($(ndims), 3)
	M_SRCS += $(TET_SRCS)
	M_INCS += $(TET_INCS)
	M_OBJS += $(TET_OBJS)
	CXXFLAGS += -DTHREED
endif

ifeq ($(energy), 1)
	CXXFLAGS += -DENERGY
endif

C3X3_DIR = 3x3-C
C3X3_LIBNAME = 3x3

ANN_DIR = ann
ANN_LIBNAME = ANN
CXXFLAGS += -I$(ANN_DIR)/include

ifeq ($(useadapt), 1)
	LIBADAPTIVITY_DIR = ./libadaptivity
	LIBADAPTIVITY_INC = $(LIBADAPTIVITY_DIR)/include
	LIBADAPTIVITY_LIB = $(LIBADAPTIVITY_DIR)/lib
	LIBADAPTIVITY_LIBNAME = adaptivity
	VTK_INC = /home/staff/echoi/opt/vtk-5.6.1/include/vtk-5.6
	CXXFLAGS += -I$(LIBADAPTIVITY_INC) -I$(VTK_INC) -DADAPT -DHAVE_VTK=1 \
    	    -I$(LIBADAPTIVITY_DIR)/adapt3d/include -I$(LIBADAPTIVITY_DIR)/metric_field/include \
        	-I$(LIBADAPTIVITY_DIR)/load_balance/include

	LIBADAPTIVITY_LIBS = $(LIBADAPTIVITY_LIB)/libadaptivity.a  -llapack -lblas -lvtkIO -lvtkGraphics -lvtkFiltering -lvtkexpat -lvtkzlib -lvtkCommon -ldl -lpthread -lm -lstdc++   -L/home/staff/echoi/opt/vtk-5.6.1/lib/vtk-5.6 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../.. -lgfortranbegin -lgfortran -lm -L/home/staff/echoi/opt/vtk-5.6.1/lib/vtk-5.6 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../.. -lgfortranbegin -lgfortran -lm -L./lib -L/home/staff/echoi/opt/vtk-5.6.1/lib/vtk-5.6 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../.. -lgfortranbegin -lgfortran -lm -lmpi_f77
endif

## Action

.PHONY: all clean take-snapshot

all: $(EXE) take-snapshot

ifeq ($(useadapt), 1)
$(EXE): $(M_OBJS) $(OBJS) $(C3X3_DIR)/lib$(C3X3_LIBNAME).a $(ANN_DIR)/lib/lib$(ANN_LIBNAME).a $(LIBADAPTIVITY_LIB)/libadaptivity.a
		$(CXX) $(M_OBJS) $(OBJS) $(LDFLAGS) $(BOOST_LDFLAGS) \
			-L$(C3X3_DIR) -l$(C3X3_LIBNAME) -L$(ANN_DIR)/lib -l$(ANN_LIBNAME) \
			$(LIBADAPTIVITY_LIBS) \
			-o $@
else
$(EXE): $(M_OBJS) $(OBJS) $(C3X3_DIR)/lib$(C3X3_LIBNAME).a $(ANN_DIR)/lib/lib$(ANN_LIBNAME).a
		$(CXX) $(M_OBJS) $(OBJS) $(LDFLAGS) $(BOOST_LDFLAGS) \
			-L$(C3X3_DIR) -l$(C3X3_LIBNAME) -L$(ANN_DIR)/lib -l$(ANN_LIBNAME) \
			-o $@
endif

take-snapshot:
	@# snapshot of the code for building the executable
	@echo Flags used to compile the code: > snapshot.diff
	@echo '  '  CXX=$(CXX) opt=$(opt) openmp=$(openmp) >> snapshot.diff
	@echo '  '  PATH=$(PATH) >> snapshot.diff
	@echo '  '  LD_LIBRARY_PATH=$(LD_LIBRARY_PATH) >> snapshot.diff
ifneq ($(HAS_HG),)
	@echo >> snapshot.diff
	@echo >> snapshot.diff
	@echo '==== Summary of the code ====' >> snapshot.diff
	@hg summary >> snapshot.diff
	@echo >> snapshot.diff
	@echo >> snapshot.diff
	@echo '== Code modification (not checked-in) ==' >> snapshot.diff
	@hg diff >> snapshot.diff
	@echo >> snapshot.diff
	@echo >> snapshot.diff
	@echo '== Code modification (checked-in but not public) ==' >> snapshot.diff
	@hg log --patch -r "draft()" >> snapshot.diff
else
	@echo \'hg\' is not in path, cannot take code snapshot. >> snapshot.diff
endif

$(OBJS): %.$(ndims)d.o : %.cxx $(INCS)
	$(CXX) $(CXXFLAGS) $(BOOST_CXXFLAGS) -c $< -o $@

$(TRI_OBJS): %.o : %.c $(TRI_INCS)
	@# Triangle cannot be compiled with -O2
	$(CXX) $(CXXFLAGS) -O1 -DTRILIBRARY -DREDUCED -DANSI_DECLARATORS -c $< -o $@

tetgen/predicates.o: tetgen/predicates.cxx $(TET_INCS)
	@# Compiling J. Shewchuk predicates, should always be
	@# equal to -O0 (no optimization). Otherwise, TetGen may not
	@# work properly.
	$(CXX) $(CXXFLAGS) -DTETLIBRARY -O0 -c $< -o $@

tetgen/tetgen.o: tetgen/tetgen.cxx $(TET_INCS)
	$(CXX) $(CXXFLAGS) -DNDEBUG -DTETLIBRARY -Wno-unused-but-set-variable -Wno-int-to-pointer-cast -c $< -o $@

$(C3X3_DIR)/lib$(C3X3_LIBNAME).a:
	@+$(MAKE) -C $(C3X3_DIR)

$(ANN_DIR)/lib/lib$(ANN_LIBNAME).a:
	@+$(MAKE) -C $(ANN_DIR) linux-g++

deepclean:
	@rm -f $(TET_OBJS) $(TRI_OBJS) $(OBJS) $(EXE)
	@+$(MAKE) -C $(C3X3_DIR) clean

clean:
	@rm -f $(OBJS) $(EXE)
