# h/t to @jimhester and @yihui for this parse block:
# https://github.com/yihui/knitr/blob/dc5ead7bcfc0ebd2789fe99c527c7d91afb3de4a/Makefile#L1-L4
# Note the portability change as suggested in the manual:
# https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Writing-portable-packages
PKGNAME = `sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION`
PKGVERS = `sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION`

ifeq ($(OS), Windows_NT)
  R_EXEC=/c/Program\ Files/R/R-4.0.2/bin/R.exe
  MINGW_DIR=/c/msys64/mingw64
  CC = $(MINGW_DIR)/bin/gcc.exe
  CXX = $(MINGW_DIR)/bin/g++.exe
  FC = $(MINGW_DIR)/bin/gfortran.exe
  AR = $(MINGW_DIR)/bin/ar.exe rv
  LD = $(MINGW_DIR)/bin/ld.exe
else
  R_EXEC = R
  CC = gcc
  CXX = g++
  FC = gfortran
  AR = ar rv
  LD = ld
endif

POLYCHORD_DIR = inst/PolyChordLite

export CC CXX FC AR LD

default: build_PolyChord build_BayesMR

build_PolyChord: 
	cd $(POLYCHORD_DIR) && $(MAKE) MPI=

build_BayesMR:
	make -C posterior

clean: clean_BayesMR clean_PolyChord

# This only cleans the BayesMR application
clean_BayesMR:
	cd posterior && rm -rf *.o && rm libBayesMR_likelihood.a

# This cleans PolyChord
clean_PolyChord:
	cd $(POLYCHORD_DIR) && $(MAKE) veryclean

