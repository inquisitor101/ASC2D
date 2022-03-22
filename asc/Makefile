MAKE = make -j 8



# Define compiler with OpenMP.
CXX = g++ -std=c++17 -fopenmp

# Define compiler serial version.
CXXSERIAL = g++ -std=c++17

# Program executable.
EXE := ACS3

# Source files directory.
SRC_DIR := src
# Object files directory.
OBJ_DIR := obj
# Binary directory.
BIN_DIR := bin
# Dependency directory.
DEP_DIR := Dependencies

# Extract source files.
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
# Specify object files.
OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))   
# Dependency files.
DEPS := $(patsubst $(SRC_DIR)/%.cpp,$(DEP_DIR)/%.d,$(SRCS))

# Compiler include files.
INCLUDE = -Iinclude

# Compiler debug flags.
CXX_DBGFLAGS = -g -Wall -DENABLE_NAN_CHECK -Werror 

# Optimization flags
CXX_OPTFLAGS = -O2 -mavx2 -mfma -ffast-math -funroll-loops -ftree-vectorize -ftree-vectorizer-verbose=2 -march=native 


# OpenMP flags.
CXX_OMP = -DHAVE_OPENMP

# Compiler flags.
CXXFLAGS = $(CXX_OPTFLAGS) $(CXX_STDFLAGS) $(CXX_OMP) $(MPI_FLAGS)  
#CXXFLAGS =  #-fconcepts 

# Linker flags.
LDFLAGS = -Llib -flto

# Blas/lapack flags.
LD_BLAS = -llapack -lblas

# Linker libraries to link.
LDLIBS = -lm $(LD_BLAS) 

# Linked.
LD := $(CXX)

# Dependency flags.
#DEPFLAGS = -MT $@ -MMD -MP -MF $(DEP_DIR)/$*.d
DEPFLAGS = -MMD -MT $@ -MP -MF $(DEP_DIR)/$*.d

# Target executable.
TARGET := $(BIN_DIR)/$(EXE)

# Preprocessor flags.
CPPFLAGS = $(INCLUDE)

# Link object files to binary.
LINK.o = $(LD) $(LDFLAGS) $(LDLIBS) 

# Compile commands.
COMPILE.cc = $(CXX) $(CXXFLAGS) $(CPPFLAGS) $(DEPFLAGS) 


#$(BIN_DIR)/$(EXE): $(OBJS) 
#	$(LINK.o) $^ -o $@

.PHONY: all clean debug
default: all


debug: CXXFLAGS = -DDEBUG $(CXX_DBGFLAGS)
debug: $(TARGET)

debugSerial: COMPILE.cc = $(CXXSERIAL) -DDEBUG $(CXX_DBGFLAGS) $(CPPFLAGS) $(DEPFLAGS) -w
debugSerial: $(TARGET)


standard: COMPILE.cc = $(CXXSERIAL) $(CPPFLAGS) $(DEPFLAGS) 
standard: $(TARGET)	


all:
	$(MAKE) $(TARGET)

$(TARGET): $(OBJS) 
	$(LINK.o) $^ -o $@


$(OBJ_DIR)/%.o:	$(SRC_DIR)/%.cpp $(DEP_DIR)/%.d | $(DEP_DIR)
	$(COMPILE.cc) -o $@ -c $<



$(DEPS):

$(DEP_DIR): ; @mkdir -p $@

.PRECIOUS: $(DEP_DIR)/%.d
	


clean:
	$(RM) $(OBJS) $(DEPS) $(TARGET)

-include $(wildcard $(DEPS))











