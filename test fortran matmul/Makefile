#
#     Copyright (c) 2017, NVIDIA CORPORATION.  All rights reserved.
#
# NVIDIA CORPORATION and its licensors retain all intellectual property
# and proprietary rights in and to this software, related documentation
# and any modifications thereto.  Any use, reproduction, disclosure or
# distribution of this software and related documentation without an express
# license agreement from NVIDIA CORPORATION is strictly prohibited.
# 

FNAME = test
FEXT = f95

NTHREADS=4
OBJ=o
EXE=out

#UNAME := $(shell uname -a)
#ifeq ($(findstring CYGWIN_NT, $(UNAME)), CYGWIN_NT)
#	OBJ=obj
#	EXE=exe
#endif

all: build run verify

build: $(FNAME).$(FEXT)
	@echo ------------------------------------ building $(FNAME).$(EXE)
#	pgf90 -fast -mp -Minfo -DPGF90 -DSYSCLK -o $(FNAME).$(EXE) $(FNAME).$(FEXT)
	pgfortran -fast -m64 -mp -Minfo -Mallocatable=03 -DSYSCLK -o $(FNAME).$(EXE) $(FNAME).$(FEXT) -lblas -llapack
# -L/home/brian/openblas/pgiopenblas/OpenBLAS-0.3.0 -lopenblas -o 

run: $(FNAME).$(EXE)
	@echo '------------- executing $(FNAME).$(EXE) with $(NTHREADS) threads'
	export OMP_NUM_THREADS=$(NTHREADS); ./$(FNAME).$(EXE)

verify: run

clean:
	@echo 'Cleaning up...'
	@rm -rf $(FNAME).$(EXE) *.dwf *.pdb prof
