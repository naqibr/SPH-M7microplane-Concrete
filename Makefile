# Compiler options
CC = g++
FC = gcc


DEBUG=no
LGFORTRAN= -lgfortran -lquadmath
# LGFORTRAN= 
DEBUGCF = -Wall -std=c++11 -g $(LGFORTRAN) -fopenmp
DEBUGFF = -Wall -g -fPIC $(LGFORTRAN)

RUNCF = -Wall -std=c++11 -fopenmp -Ofast $(LGFORTRAN)
RUNFF = -Wall -Ofast -fPIC $(LGFORTRAN)

ifeq ($(DEBUG),yes)
	CFLAGS = $(DEBUGCF)
	FFLAGS = $(DEBUGFF)
else
	CFLAGS = $(RUNCF)
	FFLAGS = $(RUNFF)
endif

LDFLAGS = -fopenmp $(LGFORTRAN)

# C++ and Fortran source files
CPP_SRCS = SPHSystem.cpp Main.cpp Vec3DMatrix3D.cpp
FORTRAN_SRCS = m7fmaterial.f

# Object files directory
OBJ_DIR = obj

# Create a list of object files from the C++ source files
CPP_OBJS = $(addprefix $(OBJ_DIR)/,$(notdir $(CPP_SRCS:.cpp=.o)))

# Create a list of object files from the Fortran source files
FORTRAN_OBJS = $(addprefix $(OBJ_DIR)/,$(notdir $(FORTRAN_SRCS:.f=.o)))

# Main target
all: Concrete

# Linking the C++ code with the Fortran object files
Concrete: $(CPP_OBJS) $(FORTRAN_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "Linking complete. Executable 'Concrete' created successfully."

# Compiling C++ source files
$(OBJ_DIR)/%.o: %.cpp | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Compiling Fortran source files with -fPIC flag
$(OBJ_DIR)/%.o: %.f | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

# Create the object directory if it doesn't exist
$(OBJ_DIR):
	mkdir $(OBJ_DIR)

# Clean up the object files and the executable
clean:
	rm -r $(OBJ_DIR)
	rm Concrete
	@echo "Cleanup complete. Removed object files and 'Concrete' executable."
