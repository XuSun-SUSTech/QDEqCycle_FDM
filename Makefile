GPU := on

CCHOME := /usr
CUDAHOME := /public/software/cuda-11.5

CC := $(CCHOME)/bin/gcc -pipe
CXX := $(CCHOME)/bin/g++ -pipe
GC := $(CUDAHOME)/bin/nvcc -rdc=true -maxrregcount=127 -arch=sm_80 #-Xptxas=-v

GCC := $(GC)

INCS := -I ./src
LIBS := -L$(CUDAHOME)/lib64 -lcudart -lcublas
INCS += -I$(CUDAHOME)/include

LIBS += -lm



CFLAGS := -x cu -c -O2 -std=c++11
LFLAGS := -O2

GCFLAGS := 

SRCDIR := ./src
OBJDIR := ./obj
BINDIR := ./bin

vpath

vpath % $(SRCDIR)
vpath % $(OBJDIR)
vpath % $(BINDIR)

DFLAGS_LIST := GPU

DFLAGS := $(foreach flag,$(DFLAGS_LIST),$(if $($(flag)),-D$(flag)))

OBJS := SetParams.o odeSolver.o Alloc.o coord.o \
        fault.o eq_cycle.o

OBJS := $(addprefix $(OBJDIR)/,$(OBJS))

$(BINDIR)/eq_cycle : $(OBJS)
	$(GCC) $(LFLAGS) $(LIBS) $^ -o $@

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(GCC) $(CFLAGS) $(DFLAGS) $(INCS) $^ -o $@

clean:
	-rm $(OBJDIR)/* -rf
	-rm $(BINDIR)/* -rf
	-rm output -rf
