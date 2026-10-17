//
// dtisawub.cpp
//    Put the orientation (UB) matrix of an ISAW/Mantid UB file (SaveIsawUB)
//    into a d*TREK header, as the unit cell and orientation angles.
//
//    The ISAW file holds the transpose of UB in the IPNS frame (x along the
//    beam, z up) with |Q| = 1/d, for the goniometer at zero.  d*TREK's
//    laboratory frame is the NeXus frame (z along the beam, y up) rotated
//    by 180 degrees about y, so a reciprocal vector (x, y, z) in IPNS
//    coordinates is (-y, z, -x) in d*TREK coordinates.
//
//    Usage: dtisawub UB.mat in.head out.head
//

//+Include files

#include <stdio.h>
#include <fstream>
#include "Dtrek.h"
#include "dtreksys.h"
#include "Cimage_header.h"
#include "Ccrystal.h"

#if !defined(VC6)
using std::cout;
using std::flush;
#endif

int
main(int argc, char *argv[])
{
  vDtrekSetModuleName("dtisawub");
  if (4 != argc)
    {
      cout << "\nUsage: dtisawub UB.mat in.head out.head\n\n"
           << "  Sets CRYSTAL_UNIT_CELL and CRYSTAL_ORIENT_ANGLES of in.head from an\n"
           << "  ISAW/Mantid UB matrix file (the goniometer setting is that of the\n"
           << "  header) and writes out.head.\n" << flush;
      return (1);
    }

  // First three lines: rows of UB^T, i.e. the columns of UB (IPNS frame)
  std::ifstream oIn(argv[1]);
  if (!oIn)
    {
      cout << "ERROR: cannot open " << argv[1] << '\n' << flush;
      return (2);
    }
  double a3x3dCols[3][3];
  for (int j = 0; j < 3; j++)
    if (!(oIn >> a3x3dCols[j][0] >> a3x3dCols[j][1] >> a3x3dCols[j][2]))
      {
        cout << "ERROR: " << argv[1] << " is not an ISAW UB file.\n" << flush;
        return (2);
      }

  // Column j of UB in d*TREK coordinates: IPNS (x, y, z) -> (-y, z, -x)
  double a9dMat[9];
  for (int j = 0; j < 3; j++)
    {
      a9dMat[3 * j + 0] = -a3x3dCols[j][1];
      a9dMat[3 * j + 1] =  a3x3dCols[j][2];
      a9dMat[3 * j + 2] = -a3x3dCols[j][0];
    }

  Cstring sIn = argv[2];
  Cimage_header oHeader(sIn);
  if (!oHeader.bIsAvailable())
    {
      cout << "ERROR: cannot read " << argv[2] << '\n' << flush;
      return (3);
    }
  Ccrystal oCrystal;
  int nSpacegroup = 1;
  if (0 == oHeader.nGetValue(Cstring("CRYSTAL_SPACEGROUP"), &nSpacegroup))
    oCrystal.m_poSpacegroup->vSet(nSpacegroup);
  if (0 != oCrystal.nSetOrientMatrix(a9dMat, 1.0))
    {
      cout << "ERROR: the UB matrix is singular.\n" << flush;
      return (4);
    }
  oCrystal.nUpdateHeader(&oHeader);
  if (0 != oHeader.nWrite(Cstring(argv[3])))
    {
      cout << "ERROR: cannot write " << argv[3] << '\n' << flush;
      return (5);
    }
  oCrystal.nList(1);
  cout << "Orientation from " << argv[1] << " written to " << argv[3] << '\n' << flush;
  return (0);
}
