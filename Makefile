# Compiler and flags
FC = gfortran
FLAGS = -llapack
# Output directory
OUT_DIR = bin

# Source directory and files
SRC_DIR = src/f90_routines
# SRC_FILES = $(wildcard $(SRC_DIR)/*.f90)
SRC_FILES = $(SRC_DIR)/diagonalize.f90

# Executable name
EXECUTABLE = $(OUT_DIR)/diagonalize

# Creats (root)/bin folder and runs compiler
all: $(EXECUTABLE)

$(EXECUTABLE): $(SRC_FILES) | $(OUT_DIR)
	$(FC) -g -o $@ $(SRC_FILES) $(FLAGS)

$(OUT_DIR):
	mkdir -p $(OUT_DIR)

# Cleans the output directory
clean:
	rm -rf $(OUT_DIR)