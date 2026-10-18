# Master Makefile for d*TREK 
#
# Copyright (C) 2017 Rigaku Americas Corporation
#                    9009 New Trails Drive
#                    The Woodlands, TX, USA  77381
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#    * Redistributions of source code must retain the above copyright
#      notice(s), this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright
#      notice(s), this list of conditions and the following disclaimer in the
#      documentation and/or other materials provided with the distribution.
#    * Neither the name of the Rigaku Americas Corporation nor the
#      names of its contributors may be used to endorse or promote products
#      derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL RIGAKU AMERICAS CORPORATION BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA OR PROFITS; OR BUSINESS INTERUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.
#

# Edit the next lines for your site
#
# Directory where OpenMotif is installed.
MOTIFHOME  := /opt/homebrew/opt/openmotif
# Directory where CBF library is located.
CBFHOME    := /Users/aaronfinke/CBFlib-0.9.6
# Directory where PTypes library is located.
PTYPESHOME := $(CURDIR)/ptypes/ptypes-2.1.1
# Directory where HDF5 is installed (NeXus TOF-Laue files; brew install hdf5).
HDF5HOME   ?= /opt/homebrew/opt/hdf5
# Directory where the gemmi C++ library is installed (dtexport writes MTZ
# files; brew install gemmi).
GEMMIHOME  ?= /opt/homebrew/opt/gemmi

# Lines below should not need to be modified

$(info )
$(info d*TREK master Makefile)

ifeq (,$(PTYPESHOME))
  $(error PTYPESHOME not found - edit Makefile with correct directory path)
endif
ifeq (,$(CBFHOME))
  $(error CBFHOME not found - edit Makefile with correct directory path)
endif

SHELL := /bin/sh
#OS    := $(echo osx)
OS    := $(shell [ "Darwin" == "`uname`" ] && echo osx || echo linux)
ARCH  := $(shell uname -m)

$(info )
$(info Operating system is $(OS))
$(info Architecture is $(ARCH))

OSTYPE_EXT.osx   = darwin
OSTYPE_EXT.linux = linux
GMAKE.osx        = /usr/bin/make
GMAKE.linux      = /usr/bin/gmake

OSTYPE     = _$(OSTYPE_EXT.$(OS))
GMAKE      = $(GMAKE.$(OS))
MFILE      = Makefile.$(OS)
INSTALLDIR = $(CURDIR)/bin$(OSTYPE)
HOMEBREW   = /opt/homebrew

CC  = /usr/bin/clang++
CXX = /usr/bin/clang++
# CC  = /usr/local/bin/g++-14
# CXX = /usr/local/bin/g++-14

$(info )
$(info Motif directory is $(MOTIFHOME))
$(info CBF library directory is $(CBFHOME))
$(info PTypes library directory is $(PTYPESHOME))
$(info HDF5 directory is $(HDF5HOME))
$(info gemmi library directory is $(GEMMIHOME))
$(info Installation directory is $(INSTALLDIR))

export OSTYPE
export CC
export CXX
export MOTIFHOME
export CBFHOME
export PTYPESHOME
export INSTALLDIR
export HOMEBREW
export HDF5HOME
export GEMMIHOME

$(info )

# Libraries that must be present, as "description|file to test"
DEPENDENCIES = \
	"OpenMotif (MOTIFHOME, brew install openmotif)|$(MOTIFHOME)/include/Xm/Xm.h" \
	"CBFlib (CBFHOME)|$(CBFHOME)/include/cbf.h" \
	"PTypes (PTYPESHOME)|$(PTYPESHOME)/include/pport.h" \
	"HDF5 (HDF5HOME, brew install hdf5)|$(HDF5HOME)/include/hdf5.h" \
	"gemmi (GEMMIHOME, brew install gemmi)|$(GEMMIHOME)/include/gemmi/mtz.hpp" \
	"gemmi library (GEMMIHOME)|$(GEMMIHOME)/lib/libgemmi_cpp.*"

.PHONY: all check clean cleanbin install

all:	check
	cd DTTREK;    $(GMAKE) -f $(MFILE) lib
	cd DTDISPLAY; $(GMAKE) -f $(MFILE) lib
	cd DTDEV;     $(GMAKE) -f $(MFILE) lib
	cd DTDISPLAY; $(GMAKE) -f $(MFILE)
	cd DTTREK;    $(GMAKE) -f $(MFILE)
	cd DTPROCESS; $(GMAKE) -f $(MFILE)

# Check that the libraries are where the lines at the top say
check:
	@status=0; \
	for dep in $(DEPENDENCIES); do \
		name=$${dep%%|*}; file=$${dep#*|}; \
		if ls $$file > /dev/null 2>&1; then echo "  found   $$name: $$file"; \
		else echo "  MISSING $$name: $$file"; status=1; fi; \
	done; \
	if [ $$status -ne 0 ]; then \
		echo "ERROR: missing libraries; install them or edit the lines at the top of the Makefile"; \
		exit 1; \
	fi

clean:
	cd DTTREK;    $(GMAKE) -f $(MFILE) clean
	cd DTDISPLAY; $(GMAKE) -f $(MFILE) clean
	cd DTDEV;     $(GMAKE) -f $(MFILE) clean
	cd DTPROCESS; $(GMAKE) -f $(MFILE) clean

cleanbin:
	cd DTTREK;    $(GMAKE) -f $(MFILE) cleanbin
	cd DTDISPLAY; $(GMAKE) -f $(MFILE) cleanbin
	cd DTPROCESS; $(GMAKE) -f $(MFILE) cleanbin

install:
	cd DTTREK;    $(GMAKE) -f $(MFILE) install
	cd DTDISPLAY; $(GMAKE) -f $(MFILE) install
	cd DTPROCESS; $(GMAKE) -f $(MFILE) install

