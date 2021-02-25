## Define compiler
CC=ifort

## Module file
MODULE= global.f90

## Define which code to compile
SOURCE= $(filter-out global.f90, $(wildcard *.f90))

## Define array checking flags
CCFLAGS=-w -fbounds-check -o

## Define array checking executable
CCEXEC=../bin/MC_check.x

## Define optimized compiling flags
OPFLAGS=-w -O2 -o

## Define optimized compiling executable
OPEXEC=MC_v2.2.x

all:
	$(CC) $(MODULE) $(SOURCE) $(OPFLAGS) $(OPEXEC)
check:
	$(CC) $(MODULE) $(SOURCE) $(CCFLAGS) $(CCEXEC)
clean:
	rm *.mod
	rm *.o
