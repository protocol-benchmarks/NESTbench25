# Makefile for engine_ising.c
#
# USAGE:
#   make standalone    # Cleans and compiles as stand-alone executable 'sim'
#   make library       # Cleans and compiles as static library 'libengine_ising.a'
#
# NOTES:
#   - Automatically removes old object, library, and executable files before building
#   - 'standalone' builds with main() uncommented
#   - 'library' builds object file and static library for linking

CXX = g++
CXXFLAGS = -Wall -O
AR = ar
ARFLAGS = rcs
LIBNAME = libengine_ising.a
SRC = engine_ising.c
OBJ = engine_ising.o
EXEC = sim

standalone:
	rm -f $(OBJ) $(LIBNAME) $(EXEC)
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(SRC) -lm

library:
	rm -f $(OBJ) $(LIBNAME) $(EXEC)
	$(CXX) -c $(SRC) -o $(OBJ)
	$(AR) $(ARFLAGS) $(LIBNAME) $(OBJ)
