#!/bin/bash
#
# makeall.sh: clean, build and install all d*TREK libraries and programs.
#
# Needs the libraries named at the top of Makefile: OpenMotif, CBFlib,
# PTypes, HDF5 and gemmi (on macOS: brew install openmotif hdf5 gemmi).
# They are checked before anything is cleaned.  Set HDF5HOME or GEMMIHOME
# in the environment, or edit Makefile, if they are not in the default
# (Homebrew) places.

[ "Darwin" = "`uname`" ] && MAKE=make || MAKE=gmake
[ "Darwin" = "`uname`" ] && BINDIR=bin_darwin || BINDIR=bin_linux
[ "Darwin" = "`uname`" ] && OS=osx || OS=linux

homeDir=`dirname $0`
[ ! "\." = "${homeDir}" ] && cd $homeDir

echo
echo Running ${MAKE} in $PWD
echo using Makefile $OS variants
echo placing results in $BINDIR
echo

${MAKE} -f Makefile check
level=$? && [ $level -ne 0 ] && exit $level
${MAKE} -f Makefile clean
level=$? && [ $level -ne 0 ] && exit $level
${MAKE} -f Makefile
level=$? && [ $level -ne 0 ] && exit $level
${MAKE} -f Makefile install
level=$? && [ $level -ne 0 ] && exit $level


exit 0
