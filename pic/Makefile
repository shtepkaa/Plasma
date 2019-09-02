# Makefile

.SUFFIXES:
.SUFFIXES: .cc .cpp .cu .o

#==============================================================================#
#  Load options
#==============================================================================#
include Makefile.in

# set shell
SHELL := /bin/bash

### # compiler standard
### STD = --std=c++11

# compiler settings
CXX = nvcc
CFLAGS = -c --compiler-options -Wall
CXXFLAGS = -c $(STD) --compiler-options -Wall
COPT = -O3

# find CUDA compute capability
SMS ?= $(shell nvcc -h | grep "compute" | sed -e "s/['=, \d40\d41]/\n/g" \
	   | grep compute | sort | uniq | cut -c 9-)

ifeq ($(GENCODE_FLAGS),)
HIGHEST_SM = $(strip $(lastword $(SMS)))
ifneq ($(HIGHEST_SM),)
GENCODE_FLAGS += -gencode arch=compute_$(CUDAARCH),code=sm_$(CUDAARCH) \
	-gencode arch=compute_$(HIGHEST_SM),code=compute_$(HIGHEST_SM) \
	-gencode arch=compute_30,code=compute_30  # default arch for compatibility
endif
endif

### NOTE: code=sm_xx means cuda binary generation for xx arch,
### code=compute_xx - PTX generation, PTX for compute_30 is compatible
### with higher arches, but binary is not

#==============================================================================#
#  Definitions
#==============================================================================#
BIN_DIR = bin
BUILD_DIR = build
CONF_DIR = conf
SRC_DIR = src
TEST_DIR = test

UTILITY_DIR = $(SRC_DIR)/utility
BLAS_DIR = $(SRC_DIR)/blas
TENSOR_DIR = $(SRC_DIR)/tensor
OPERATOR_DIR = $(SRC_DIR)/operator
SCHEME_DIR = $(SRC_DIR)/scheme
SOLVER_DIR = $(SRC_DIR)/solver

SOLVER = $(wildcard $(SOLVER_DIR)/*.cpp)
SRCS = $(filter-out $(SOLVER), $(shell find $(SRC_DIR)/ -name '*.cpp'))
OBJS = $(patsubst %.cpp, build/%.o, $(notdir $(SRCS)))

LIBTT = $(BUILD_DIR)/libtt.a

# lib
AR = ar
LIBPATH = ./lib/lib.a

SRCDIR = ./src

# define sources
CUSOURCES = $(filter-out $(SRCDIR)/test.cu $(SRCDIR)/autolykos.cu, \
			$(wildcard $(SRCDIR)/*.cu))
CPPSOURCES = $(wildcard $(SRCDIR)/*.cc) $(wildcard $(SRCDIR)/bip39/*.cc)
CSOURCES = $(wildcard $(SRCDIR)/*.c)

# define objects
OBJECTS = $(CUSOURCES:.cu=.o) $(CPPSOURCES:.cc=.o) $(CSOURCES:.c=.o)

# define executables
AUTOEXEC = auto.out
TESTEXEC = test.out

# compile objects
%.o: %.cu
	$(CXX) $(COPT) $(CXXFLAGS) $(GENCODE_FLAGS) --maxrregcount $(MAXREG) \
		-DBLOCK_DIM=$(BLOCKDIM) -DNONCES_PER_ITER=$(WORKSPACE) $< -o $@
%.o: %.cc
	$(CXX) $(COPT) $(CXXFLAGS) $(EMBED) $< -o $@
%.o: %.c
	$(CXX) $(COPT) $(CFLAGS) $< -o $@

# default (miner executable)
all: clean lib autoexec

# lib
lib: $(OBJECTS)
	mkdir -p ./lib;
	$(AR) rc $(LIBPATH) $(OBJECTS)
	ranlib $(LIBPATH)

# test executable if lib made
testexec:
	$(CXX) $(SRCDIR)/test.cu $(LIBPATH) $(LIBS) $(COPT) $(STD) \
		$(GENCODE_FLAGS) -DBLOCK_DIM=$(BLOCKDIM) \
		-DNONCES_PER_ITER=$(WORKSPACE) -o $(TESTEXEC)

#==============================================================================#
#   Default
#==============================================================================#
all: libtt

#==============================================================================#
#   Blas
#==============================================================================#
$(BLAS_DIR)/blas.h: $(UTILITY_DIR)/definitions.h
	touch $@

$(BUILD_DIR)/%.o: $(BLAS_DIR)/%.cpp $(BLAS_DIR)/blas.h
	$(CXX) $(CFLAGS) $(COPT) $< -o $@

#==============================================================================#
#	Libtt
#==============================================================================#
build:
	mkdir -p $(BUILD_DIR)

libtt: build $(OBJS)
	$(AR) rc $(LIBTT) $(OBJS)
	ranlib $(LIBTT)

#==============================================================================#
#	Binary fragmentation test
#==============================================================================#
BINARY_FRAG_TEST = $(TEST_DIR)/binary_frag_test/binary_frag_test

$(BINARY_FRAG_TEST).h:
	touch $@

binary_frag_test: libtt $(BINARY_FRAG_TEST).h
	mkdir -p $(BIN_DIR)
	$(CXX) $(BINARY_FRAG_TEST).cpp $(LIBTT) $(LIBS) $(COPT) -o \
		$(BIN_DIR)/binary_frag_test.out

#==============================================================================#
#	Integration test
#==============================================================================#
INTEGRATION_TEST = $(TEST_DIR)/integration_test/integration_test

$(INTEGRATION_TEST).h:
	touch $@

integration_test: libtt $(INTEGRATION_TEST).h
	mkdir -p $(BIN_DIR)
	$(CXX) $(INTEGRATION_TEST).cpp $(LIBTT) $(LIBS) $(COPT) -o \
		$(BIN_DIR)/integration_test.out

#==============================================================================#
#	Solver
#==============================================================================#
$(CONF_DIR)/user_definitions.h:
	touch $@

solver: libtt $(CONF_DIR)/user_definitions.h
	mkdir -p $(BIN_DIR)
	$(CXX) $(SOLVER_DIR)/solver.cpp $(LIBTT) $(LIBS) $(COPT) -o \
		$(BIN_DIR)/solver.out

#==============================================================================#
#	Brill Solver
#==============================================================================#
$(CONF_DIR)/brill_definitions.h:
	touch $@

brill_solver: libtt $(CONF_DIR)/brill_definitions.h
	mkdir -p $(BIN_DIR)
	$(CXX) $(SOLVER_DIR)/brill_solver.cpp $(LIBTT) $(LIBS) $(COPT) -o \
		$(BIN_DIR)/brill_solver.out

#==============================================================================#
#	Clean
#==============================================================================#
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

purge: clean
	find ./ -name data -exec rm -rf {} +;
	find ./ -name checkpoints -exec rm -rf {} +

.PHONY: binary_frag_test brill_solver build clean integration_test libtt purge \
	solver

# Makefile
