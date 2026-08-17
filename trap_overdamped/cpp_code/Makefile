# Makefile for engine_trap_overdamped.cpp
#
# USAGE:
#   make standalone    # Cleans and compiles as stand-alone executable 'sim'
#   make library       # Cleans and compiles as static library 'libengine_trap_overdamped.a'
#   make ga            # Builds the genetic-algorithm driver 'ga' and worker 'ga_process'
#   make all           # standalone + ga
#
# NOTES:
#   - Automatically removes old object, library, and executable files before building
#   - 'sim' accepts command-line options; run './sim -h' for usage
#   - 'library' builds object file and static library for linking

CXX = g++
CXXFLAGS = -Wall -O
AR = ar
ARFLAGS = rcs
LIBNAME = libengine_trap_overdamped.a
SRC = engine_trap_overdamped.cpp
OBJ = engine_trap_overdamped.o
EXEC = sim

all: standalone ga

standalone:
	rm -f $(OBJ) $(LIBNAME) $(EXEC)
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(SRC) -lm

library:
	rm -f $(OBJ) $(LIBNAME) $(EXEC)
	$(CXX) $(CXXFLAGS) -DENGINE_LIBRARY -c $(SRC) -o $(OBJ)
	$(AR) $(ARFLAGS) $(LIBNAME) $(OBJ)

# -mcmodel=medium is required for the large static arrays on x86_64 Linux,
# but is unsupported on other targets (e.g. ARM macOS)
UNAME_M := $(shell uname -m)
ifeq ($(UNAME_M),x86_64)
MCMODEL = -mcmodel=medium
else
MCMODEL =
endif

ga: library
	$(CXX) $(CXXFLAGS) -o ga ga.cpp -lm
	$(CXX) $(CXXFLAGS) $(MCMODEL) -o ga_process ga_process.cpp -L. -lengine_trap_overdamped -lm

clean:
	rm -f $(OBJ) $(LIBNAME) $(EXEC) ga ga_process report_*.dat

.PHONY: all standalone library ga clean
