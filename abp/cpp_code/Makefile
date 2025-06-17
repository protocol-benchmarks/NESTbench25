# Makefile for engine_abp.c
#
# USAGE:
#   make standalone    # Cleans and compiles as stand-alone executable 'sim'
#   make library       # Cleans and compiles as static library 'libengine_abp.a'
#
# NOTES:
#   - Automatically removes old object, library, and executable files before building
#   - 'standalone' builds with main() uncommented
#   - 'library' builds object file and static library for linking

CXX = g++
CXXFLAGS = -Wall -O
AR = ar
ARFLAGS = rcs
LIBNAME = libengine_abp.a
SRC = engine_abp.cpp
OBJ = engine_abp.o
EXEC = sim

standalone:
	rm -f $(OBJ) $(LIBNAME) $(EXEC)
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(SRC) -lm

library:
	rm -f $(OBJ) $(LIBNAME) $(EXEC)
	$(CXX) -c $(SRC) -o $(OBJ)
	$(AR) $(ARFLAGS) $(LIBNAME) $(OBJ)
