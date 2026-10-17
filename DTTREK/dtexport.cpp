//
// dtexport.cpp
//    Convert LAUENORM input files ("afout" files, as written by dtref2afout
//    or by Mantid's SaveLauenorm) into one unmerged multi-record MTZ file for
//    AIMLESS, including its Laue version:
//
//       afout:  h k l lambda theta I sigma     (%5d%5d%5d%10.5f%10.5f%10d%10d)
//       MTZ:    H K L M/ISYM BATCH I SIGI LAMBDA [ROT]
//
//    Each afout file becomes one batch.  The indices are reduced to the
//    asymmetric unit with M/ISYM set, and the file is sorted, as POINTLESS
//    writes unmerged data.  ROT is written as 0 (a stationary crystal), so
//    that AIMLESS does not take the batch number as the rotation angle.
//    The afout format holds no cell or space group: they come from a d*TREK
//    header (-head) or from -cell and -spacegroup.
//
//    With the refined d*TREK header of each run (given among the input
//    files, in the same order as the afout files), each batch header also
//    gets the crystal orientation, so that AIMLESS can apply its secondary
//    beam and absorption corrections.  The header frame is the d*TREK
//    laboratory frame:
//       UMAT    the rotation U with U B = G (U B)_dtrek, B being the AIMLESS
//               (Busing & Levy) B matrix of the space-group-constrained cell
//               and G the crystal goniometer matrix of the run; U is made
//               orthogonal, which absorbs a refined cell that does not obey
//               the space group
//       SOURCE  the d*TREK source vector S0 (anti-parallel to the beam)
//       E1      the outermost goniometer axis; one axis at datum 0
//    with a nominal rotation range of -0.5 to 0.5 degrees, since AIMLESS
//    ignores the orientation of a batch without one.  The wavelengths the
//    header predicts for the reflections of its afout file are checked, so
//    a header paired with the wrong file is caught.
//
//    Uses the gemmi library to write the MTZ file.
//
//    Usage: dtexport [options] [run1.head ...] run1.afout [run2.afout ...]
//

//+Include files

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <gemmi/mtz.hpp>
#include <gemmi/symmetry.hpp>
#include <gemmi/unitcell.hpp>
#include "Dtrek.h"
#include "dtreksys.h"
#include "dtrekvec.h"
#include "Cimage_header.h"
#include "Ccrystal.h"
#include "Cgoniometer.h"
#include "Csource.h"

using std::cout;
using std::flush;
using std::string;
using std::vector;

static void
vUsage(void)
{
  cout << "\nUsage: dtexport [options] [run1.head ...] run1.afout [run2.afout ...]\n\n"
       << "  Writes the afout files (h k l lambda theta I sigma, the LAUENORM input\n"
       << "  written by dtref2afout) as one unmerged MTZ file for AIMLESS, with one\n"
       << "  batch per file and the columns H K L M/ISYM BATCH I SIGI LAMBDA ROT.\n\n"
       << "  d*TREK headers among the input files (the refined header of each run,\n"
       << "  one per afout file and in the same order) put the crystal orientation\n"
       << "  of the run into its batch header, for the AIMLESS secondary beam and\n"
       << "  absorption corrections.  The wavelengths each header predicts for its\n"
       << "  afout file are checked.\n\n"
       << "  -out file          Output MTZ file (default dtexport.mtz).\n"
       << "  -head file         d*TREK header with the cell (CRYSTAL_UNIT_CELL) and\n"
       << "                     space group (CRYSTAL_SPACEGROUP).  Default: the\n"
       << "                     first orientation header.\n"
       << "  -cell a b c al be ga\n"
       << "                     Unit cell (overrides -head).  It is made to obey\n"
       << "                     the space group (e.g. a = b, gamma = 120 for\n"
       << "                     hexagonal), with a warning if it changes.\n"
       << "  -spacegroup sg     Space group number or name, e.g. 179 or P6522\n"
       << "                     (overrides -head).\n"
       << "  -batch n           Batch number of the first file (default 1); the\n"
       << "                     files get consecutive numbers in the given order.\n"
       << "  -wave lmin lmax    Keep only these wavelengths (A).\n"
       << "  -scale f           Multiply I and sigma by f (default 1).\n"
       << "  -norot             Do not write the ROT column.\n"
       << "  -project name, -crystal name, -dataset name\n"
       << "                     MTZ dataset names (default dtrek, xtal, laue).\n"
       << "  -title text        MTZ title.\n\n"
       << "  Reflections with sigma <= 0 are skipped.  dtref2afout writes negative\n"
       << "  intensities as 0 unless given -keepneg; use -keepneg for AIMLESS.\n"
       << flush;
}

// One afout record.  Whitespace-separated fields first (some legacy files
// were edited by hand), then the fixed columns (wide integers can run into
// each other).
static bool
bParseAfout(const string& rsLine, int a3nHKL[3], double *pdLambda,
            double *pdTheta, double *pdI, double *pdSigma)
{
  std::istringstream oIn(rsLine);
  string sExtra;
  if (   (oIn >> a3nHKL[0] >> a3nHKL[1] >> a3nHKL[2] >> *pdLambda >> *pdTheta >> *pdI >> *pdSigma)
      && !(oIn >> sExtra) )
    return (true);

  // %5d%5d%5d%10.5f%10.5f%10d%10d
  static const int a7nStart[7] = {0, 5, 10, 15, 25, 35, 45};
  static const int a7nWidth[7] = {5, 5, 5, 10, 10, 10, 10};
  if (55 > (int) rsLine.size())
    return (false);
  double a7dValue[7];
  for (int i = 0; i < 7; i++)
    {
      string sField = rsLine.substr(a7nStart[i], a7nWidth[i]);
      char  *pcEnd  = NULL;
      a7dValue[i] = strtod(sField.c_str(), &pcEnd);
      if (pcEnd == sField.c_str())
        return (false);
    }
  for (int i = 0; i < 3; i++)
    a3nHKL[i] = (int) a7dValue[i];
  *pdLambda = a7dValue[3];
  *pdTheta  = a7dValue[4];
  *pdI      = a7dValue[5];
  *pdSigma  = a7dValue[6];
  return (true);
}

// Is the file a d*TREK header ("{" then HEADER_BYTES)?
static bool
bIsDtrekHeader(const string& rsFile)
{
  std::ifstream oIn(rsFile.c_str(), std::ios::binary);
  char a64cStart[64] = {0};
  oIn.read(a64cStart, sizeof(a64cStart) - 1);
  string sStart(a64cStart, (size_t) oIn.gcount());
  size_t nPos = sStart.find_first_not_of(" \t\r\n");
  return (   (string::npos != nPos) && ('{' == sStart[nPos])
          && (string::npos != sStart.find("HEADER_BYTES")) );
}

// Value of a keyword of a d*TREK header ("KEY=value;"), or "".
static string
sHeaderValue(const string& rsHeader, const string& rsKey)
{
  size_t nPos = 0;
  while (string::npos != (nPos = rsHeader.find(rsKey + "=", nPos)))
    {
      if ( (0 == nPos) || ('\n' == rsHeader[nPos - 1]) || ('\r' == rsHeader[nPos - 1]) )
        {
          size_t nStart = nPos + rsKey.size() + 1;
          size_t nEnd   = rsHeader.find(';', nStart);
          if (string::npos == nEnd)
            return ("");
          return (rsHeader.substr(nStart, nEnd - nStart));
        }
      nPos += rsKey.size();
    }
  return ("");
}

static int
nReadHeader(const string& rsFile, double a6dCell[6], bool *pbCell, string *psSG)
{
  std::ifstream oIn(rsFile.c_str(), std::ios::binary);
  if (!oIn)
    {
      cout << "ERROR: cannot open " << rsFile << '\n' << flush;
      return (1);
    }
  string sHeader(65536, '\0');
  oIn.read(&sHeader[0], (std::streamsize) sHeader.size());
  sHeader.resize((size_t) oIn.gcount());
  size_t nEnd = sHeader.find('}');
  if ( (string::npos == sHeader.find("HEADER_BYTES")) || (string::npos == nEnd) )
    {
      cout << "ERROR: " << rsFile << " is not a d*TREK header\n" << flush;
      return (1);
    }
  sHeader.resize(nEnd);

  string sCell = sHeaderValue(sHeader, "CRYSTAL_UNIT_CELL");
  if ("" != sCell)
    {
      std::istringstream oCell(sCell);
      int i = 0;
      while ( (i < 6) && (oCell >> a6dCell[i]) )
        i++;
      *pbCell = (6 == i);
    }
  string sSG = sHeaderValue(sHeader, "CRYSTAL_SPACEGROUP");
  if ("" != sSG)
    {
      std::istringstream oSG(sSG);
      oSG >> *psSG;
    }
  return (0);
}

// Make the cell obey the crystal system of the space group.
static gemmi::UnitCell
oConstrainCell(const double a6dCell[6], const gemmi::SpaceGroup& roSG)
{
  double a = a6dCell[0], b = a6dCell[1], c = a6dCell[2];
  double al = a6dCell[3], be = a6dCell[4], ga = a6dCell[5];
  switch (roSG.crystal_system())
    {
    case gemmi::CrystalSystem::Triclinic:
      break;
    case gemmi::CrystalSystem::Monoclinic:
      {
        char cAxis = roSG.monoclinic_unique_axis();
        if ('a' != cAxis) al = 90.0;
        if ('b' != cAxis) be = 90.0;
        if ('c' != cAxis) ga = 90.0;
      }
      break;
    case gemmi::CrystalSystem::Orthorhombic:
      al = be = ga = 90.0;
      break;
    case gemmi::CrystalSystem::Tetragonal:
      a = b = 0.5 * (a + b);
      al = be = ga = 90.0;
      break;
    case gemmi::CrystalSystem::Trigonal:
      if ('R' == roSG.ext)        // Rhombohedral axes
        {
          a = b = c = (a + b + c) / 3.0;
          al = be = ga = (al + be + ga) / 3.0;
          break;
        }
      [[fallthrough]];          // Hexagonal axes
    case gemmi::CrystalSystem::Hexagonal:
      a = b = 0.5 * (a + b);
      al = be = 90.0;
      ga = 120.0;
      break;
    case gemmi::CrystalSystem::Cubic:
      a = b = c = (a + b + c) / 3.0;
      al = be = ga = 90.0;
      break;
    }
  return gemmi::UnitCell(a, b, c, al, be, ga);
}

// 3x3 matrices below are row-major, M[row][column], M v = sum_j M[i][j] v[j].

static void
vMatMul(const double A[3][3], const double B[3][3], double C[3][3])
{
  double T[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      T[i][j] = A[i][0] * B[0][j] + A[i][1] * B[1][j] + A[i][2] * B[2][j];
  memcpy(C, T, sizeof(T));
}

static void
vMatVec(const double A[3][3], const double v[3], double w[3])
{
  for (int i = 0; i < 3; i++)
    w[i] = A[i][0] * v[0] + A[i][1] * v[1] + A[i][2] * v[2];
}

static double
dMatInv(const double A[3][3], double B[3][3])
{
  double dDet =   A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
                - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
                + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
  if (0.0 == dDet)
    return (0.0);
  B[0][0] =  (A[1][1] * A[2][2] - A[1][2] * A[2][1]) / dDet;
  B[0][1] = -(A[0][1] * A[2][2] - A[0][2] * A[2][1]) / dDet;
  B[0][2] =  (A[0][1] * A[1][2] - A[0][2] * A[1][1]) / dDet;
  B[1][0] = -(A[1][0] * A[2][2] - A[1][2] * A[2][0]) / dDet;
  B[1][1] =  (A[0][0] * A[2][2] - A[0][2] * A[2][0]) / dDet;
  B[1][2] = -(A[0][0] * A[1][2] - A[0][2] * A[1][0]) / dDet;
  B[2][0] =  (A[1][0] * A[2][1] - A[1][1] * A[2][0]) / dDet;
  B[2][1] = -(A[0][0] * A[2][1] - A[0][1] * A[2][0]) / dDet;
  B[2][2] =  (A[0][0] * A[1][1] - A[0][1] * A[1][0]) / dDet;
  return (dDet);
}

// The AIMLESS (Scell::Bmat, Busing & Levy) B matrix: a* along x, in 1/A
static void
vBmatrix(const gemmi::UnitCell& roCell, double B[3][3])
{
  double dSinBetaR  = sqrt(1.0 - roCell.cos_betar * roCell.cos_betar);
  double dSinGammaR = sqrt(1.0 - roCell.cos_gammar * roCell.cos_gammar);
  double dCosAlpha  = cos(roCell.alpha * 3.14159265358979323846 / 180.0);
  B[0][0] = roCell.ar;  B[0][1] = roCell.br * roCell.cos_gammar;  B[0][2] = roCell.cr * roCell.cos_betar;
  B[1][0] = 0.0;        B[1][1] = roCell.br * dSinGammaR;         B[1][2] = -roCell.cr * dSinBetaR * dCosAlpha;
  B[2][0] = 0.0;        B[2][1] = 0.0;                            B[2][2] = 1.0 / roCell.c;
}

// Nearest rotation to a nearly orthogonal matrix (polar decomposition by
// the iteration R <- (R + R^-T) / 2).  Returns max |M^T M - I| of the input.
static double
dOrthonormalize(const double M[3][3], double R[3][3])
{
  double dDev = 0.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      {
        double dDot = M[0][i] * M[0][j] + M[1][i] * M[1][j] + M[2][i] * M[2][j];
        dDev = std::max(dDev, fabs(dDot - ((i == j) ? 1.0 : 0.0)));
      }
  memcpy(R, M, sizeof(double) * 9);
  for (int nIter = 0; nIter < 50; nIter++)
    {
      double Inv[3][3];
      if (0.0 == dMatInv(R, Inv))
        break;
      double dChange = 0.0;
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          {
            double dNew = 0.5 * (R[i][j] + Inv[j][i]);
            dChange = std::max(dChange, fabs(dNew - R[i][j]));
            R[i][j] = dNew;
          }
      if (1.0e-12 > dChange)
        break;
    }
  return (dDev);
}

// Laue wavelength of reciprocal vector X for source vector S0 (toward the
// source): lambda = 2 X.S0 / |X|^2
static double
dLaueLambda(const double X[3], const double S0[3])
{
  double dXX = X[0] * X[0] + X[1] * X[1] + X[2] * X[2];
  if (0.0 >= dXX)
    return (0.0);
  return (2.0 * (X[0] * S0[0] + X[1] * S0[1] + X[2] * S0[2]) / dXX);
}

// Crystal orientation of one run, from its refined d*TREK header
struct tOrient
{
  string sFile;
  double a3x3dA[3][3];   // Reciprocal vector of h in the lab frame: A h (1/A)
  double a3dS0[3];       // Source vector (unit, anti-parallel to the beam)
  double a3dE1[3];       // Outermost goniometer axis (unit)
  string sAxis;          // Its name
  double a6dCell[6];
  int    nSG;
};

static int
nReadOrientation(const string& rsFile, tOrient *poOrient)
{
  poOrient->sFile = rsFile;
  Cimage_header oHeader(Cstring(rsFile.c_str()));
  if (!oHeader.bIsAvailable())
    {
      cout << "ERROR: cannot read the d*TREK header " << rsFile << '\n' << flush;
      return (1);
    }
  Ccrystal oCrystal(oHeader);
  if (!oCrystal.bIsAvailable())
    {
      cout << "ERROR: no crystal orientation in " << rsFile << '\n' << flush;
      return (1);
    }
  Cgoniometer oGonio(oHeader, Ccrystal::ms_sCrystalPrefix);
  if (!oGonio.bIsAvailable())
    {
      cout << "ERROR: no crystal goniometer in " << rsFile << '\n' << flush;
      return (1);
    }
  Csource oSource(oHeader);
  if (!oSource.bIsAvailable())
    {
      cout << "ERROR: no source in " << rsFile << '\n' << flush;
      return (1);
    }

  // d*TREK: X = G (U B) h, with lambda = 1 so that X is in 1/A.  Columns of
  // A from the d*TREK routines, so that their storage order does not matter.
  double a3x3dUB[3][3], a3x3dG[3][3];
  (void) oCrystal.nCalcOrientMatrix(1.0);
  oCrystal.vGetOrientMatrix(&a3x3dUB[0][0]);
  oGonio.vCalcGetRotMatrix(&a3x3dG[0][0], -1);
  for (int j = 0; j < 3; j++)
    {
      double a3dE[3] = {0.0, 0.0, 0.0}, a3dT[3], a3dX[3];
      a3dE[j] = 1.0;
      vMulMat3DVec3D(a3x3dUB, a3dE, a3dT);
      vMulMat3DVec3D(a3x3dG, a3dT, a3dX);
      for (int i = 0; i < 3; i++)
        poOrient->a3x3dA[i][j] = a3dX[i];
    }

  oSource.vCalcGetS0(poOrient->a3dS0);
  float a3fE1[3] = {0.0f, 0.0f, 0.0f};
  if (0 != oGonio.nGetRotVector(0, a3fE1))
    {
      cout << "ERROR: no goniometer axis in " << rsFile << '\n' << flush;
      return (1);
    }
  double dS0 = 0.0, dE1 = 0.0;
  for (int i = 0; i < 3; i++)
    {
      poOrient->a3dE1[i] = a3fE1[i];
      dS0 += poOrient->a3dS0[i] * poOrient->a3dS0[i];
      dE1 += poOrient->a3dE1[i] * poOrient->a3dE1[i];
    }
  if ( (0.0 == dS0) || (0.0 == dE1) )
    {
      cout << "ERROR: zero source vector or goniometer axis in " << rsFile << '\n' << flush;
      return (1);
    }
  for (int i = 0; i < 3; i++)
    {
      poOrient->a3dS0[i] /= sqrt(dS0);
      poOrient->a3dE1[i] /= sqrt(dE1);
    }
  poOrient->sAxis = oGonio.sGetName(0).string();
  for (size_t k = 0; k < poOrient->sAxis.size(); k++)
    poOrient->sAxis[k] = (char) toupper((unsigned char) poOrient->sAxis[k]);
  if ("" == poOrient->sAxis)
    poOrient->sAxis = "OMEGA";

  float a6fCell[6];
  oCrystal.vGetCell(a6fCell);
  for (int i = 0; i < 6; i++)
    poOrient->a6dCell[i] = a6fCell[i];
  poOrient->nSG = oCrystal.m_poSpacegroup->nGet();
  return (0);
}

struct tBatchFile
{
  string sFile;
  int    nBatch;
  long   lNumRead;
  long   lNumUsed;
  long   lNumZero;
  long   lNumNeg;
  double dLambdaSum;
  double a2dLambda[2];
  vector<int>    vnHKL;      // Original indices and wavelengths, for the
  vector<double> vdLambda;   // check against the orientation header

  // Orientation (if a header was given)
  bool   bOrient;
  double a3x3dU[3][3];
  double dOrthDev;           // max |M^T M - I| before orthogonalization
  double dLamRms;            // rms relative wavelength error: header UB
  double dLamRmsU;           // ... with U and the constrained cell
};

int
main(int argc, char *argv[])
{
  vDtrekSetModuleName("dtexport");
  string sOut     = "dtexport.mtz";
  string sHead    = "";
  string sSG      = "";
  string sProject = "dtrek";
  string sCrystal = "xtal";
  string sDataset = "laue";
  string sTitle   = "";
  double a6dCell[6];
  bool   bCell    = false;
  int    nBatch0  = 1;
  double a2dWave[2] = {0.0, 0.0};
  double dScale   = 1.0;
  bool   bRot     = true;
  vector<string> vsIn, vsHeads;

  for (int i = 1; i < argc; i++)
    {
      string sArg  = argv[i];
      int    nMore = argc - 1 - i;
      if ( ("-out" == sArg) && (1 <= nMore) )
        sOut = argv[++i];
      else if ( ("-head" == sArg) && (1 <= nMore) )
        sHead = argv[++i];
      else if ( ("-cell" == sArg) && (6 <= nMore) )
        {
          for (int j = 0; j < 6; j++)
            a6dCell[j] = atof(argv[++i]);
          bCell = true;
        }
      else if ( ("-spacegroup" == sArg) && (1 <= nMore) )
        sSG = argv[++i];
      else if ( ("-batch" == sArg) && (1 <= nMore) )
        nBatch0 = atoi(argv[++i]);
      else if ( ("-wave" == sArg) && (2 <= nMore) )
        {
          a2dWave[0] = atof(argv[++i]);
          a2dWave[1] = atof(argv[++i]);
          if (a2dWave[0] > a2dWave[1])
            std::swap(a2dWave[0], a2dWave[1]);
        }
      else if ( ("-scale" == sArg) && (1 <= nMore) )
        dScale = atof(argv[++i]);
      else if ("-norot" == sArg)
        bRot = false;
      else if ( ("-project" == sArg) && (1 <= nMore) )
        sProject = argv[++i];
      else if ( ("-crystal" == sArg) && (1 <= nMore) )
        sCrystal = argv[++i];
      else if ( ("-dataset" == sArg) && (1 <= nMore) )
        sDataset = argv[++i];
      else if ( ("-title" == sArg) && (1 <= nMore) )
        sTitle = argv[++i];
      else if ( ("-help" == sArg) || ("-h" == sArg) )
        {
          vUsage();
          return (0);
        }
      else if ('-' == sArg[0])
        {
          cout << "ERROR: unknown or incomplete option " << sArg << '\n' << flush;
          vUsage();
          return (1);
        }
      else if (bIsDtrekHeader(sArg))
        vsHeads.push_back(sArg);
      else
        vsIn.push_back(sArg);
    }
  if (vsIn.empty())
    {
      vUsage();
      return (1);
    }
  if ( !vsHeads.empty() && (vsHeads.size() != vsIn.size()) )
    {
      cout << "ERROR: " << vsHeads.size() << " orientation headers for " << vsIn.size()
           << " afout files; give one header per afout file, in the same order.\n" << flush;
      return (1);
    }

  vector<tOrient> voOrient(vsHeads.size());
  for (size_t n = 0; n < vsHeads.size(); n++)
    if (0 != nReadOrientation(vsHeads[n], &voOrient[n]))
      return (2);

  // Cell and space group: -cell and -spacegroup, then -head, then the
  // first orientation header
  if ("" != sHead)
    {
      double a6dHeadCell[6];
      bool   bHeadCell = false;
      string sHeadSG   = "";
      if (0 != nReadHeader(sHead, a6dHeadCell, &bHeadCell, &sHeadSG))
        return (2);
      if (!bCell && bHeadCell)
        {
          std::copy(a6dHeadCell, a6dHeadCell + 6, a6dCell);
          bCell = true;
        }
      if ("" == sSG)
        sSG = sHeadSG;
    }
  if (!voOrient.empty())
    {
      if (!bCell)
        {
          std::copy(voOrient[0].a6dCell, voOrient[0].a6dCell + 6, a6dCell);
          bCell = true;
        }
      if ( ("" == sSG) && (0 < voOrient[0].nSG) )
        sSG = std::to_string(voOrient[0].nSG);
    }
  if (!bCell || ("" == sSG))
    {
      cout << "ERROR: the cell and space group are needed (-head, or -cell and -spacegroup)\n"
           << flush;
      return (1);
    }
  const gemmi::SpaceGroup *poSG = NULL;
  char *pcEnd = NULL;
  long lSG = strtol(sSG.c_str(), &pcEnd, 10);
  if ( ('\0' == *pcEnd) && (0 < lSG) )
    poSG = gemmi::find_spacegroup_by_number((int) lSG);
  else
    poSG = gemmi::find_spacegroup_by_name(sSG);
  if (NULL == poSG)
    {
      cout << "ERROR: unknown space group " << sSG << '\n' << flush;
      return (1);
    }
  gemmi::UnitCell oCell = oConstrainCell(a6dCell, *poSG);
  double a6dUsed[6] = {oCell.a, oCell.b, oCell.c, oCell.alpha, oCell.beta, oCell.gamma};
  double dMaxRel = 0.0, dMaxAng = 0.0;
  for (int j = 0; j < 3; j++)
    {
      dMaxRel = std::max(dMaxRel, fabs(a6dUsed[j] - a6dCell[j]) / a6dCell[j]);
      dMaxAng = std::max(dMaxAng, fabs(a6dUsed[j + 3] - a6dCell[j + 3]));
    }
  if ( (1.0e-5 < dMaxRel) || (1.0e-4 < dMaxAng) )
    {
      char a200cLine[200];
      snprintf(a200cLine, sizeof(a200cLine),
               "WARNING: cell changed to obey %s (%s): %.4f %.4f %.4f %.3f %.3f %.3f\n",
               poSG->xhm().c_str(), poSG->crystal_system_str(),
               oCell.a, oCell.b, oCell.c, oCell.alpha, oCell.beta, oCell.gamma);
      cout << a200cLine << flush;
      if ( (0.01 < dMaxRel) || (1.0 < dMaxAng) )
        cout << "WARNING: that is a change of more than 1 % or 1 degree; check the cell.\n"
             << flush;
    }

  // Read the afout files
  const int nCols = bRot ? 9 : 8;
  vector<float>      vfData;
  vector<tBatchFile> voFiles;
  for (size_t n = 0; n < vsIn.size(); n++)
    {
      std::ifstream oIn(vsIn[n].c_str());
      if (!oIn)
        {
          cout << "ERROR: cannot open " << vsIn[n] << '\n' << flush;
          return (3);
        }
      tBatchFile oF;
      oF.sFile        = vsIn[n];
      oF.nBatch       = nBatch0 + (int) n;
      oF.lNumRead     = oF.lNumUsed = oF.lNumZero = oF.lNumNeg = 0;
      oF.dLambdaSum   = 0.0;
      oF.a2dLambda[0] = 1.0e10;
      oF.a2dLambda[1] = 0.0;
      oF.bOrient      = !voOrient.empty();
      oF.dOrthDev = oF.dLamRms = oF.dLamRmsU = 0.0;
      string sLine;
      int    nLine = 0;
      while (std::getline(oIn, sLine))
        {
          nLine++;
          if (sLine.find_first_not_of(" \t\r") == string::npos)
            continue;
          int    a3nHKL[3];
          double dLambda, dTheta, dI, dSigma;
          if (!bParseAfout(sLine, a3nHKL, &dLambda, &dTheta, &dI, &dSigma))
            {
              cout << "ERROR: " << vsIn[n] << " line " << nLine
                   << " is not an afout record:\n  " << sLine << '\n' << flush;
              return (3);
            }
          oF.lNumRead++;
          if (0.0 >= dSigma)
            continue;
          if ( (0.0 < a2dWave[1]) && ( (dLambda < a2dWave[0]) || (dLambda > a2dWave[1]) ) )
            continue;
          if ( (0 == a3nHKL[0]) && (0 == a3nHKL[1]) && (0 == a3nHKL[2]) )
            continue;
          if (0.0 == dI)
            oF.lNumZero++;
          else if (0.0 > dI)
            oF.lNumNeg++;
          oF.lNumUsed++;
          oF.dLambdaSum  += dLambda;
          oF.a2dLambda[0] = std::min(oF.a2dLambda[0], dLambda);
          oF.a2dLambda[1] = std::max(oF.a2dLambda[1], dLambda);
          if (oF.bOrient)
            {
              oF.vnHKL.insert(oF.vnHKL.end(), a3nHKL, a3nHKL + 3);
              oF.vdLambda.push_back(dLambda);
            }
          vfData.push_back((float) a3nHKL[0]);
          vfData.push_back((float) a3nHKL[1]);
          vfData.push_back((float) a3nHKL[2]);
          vfData.push_back(1.0f);                       // M/ISYM, set below
          vfData.push_back((float) oF.nBatch);
          vfData.push_back((float) (dI * dScale));
          vfData.push_back((float) (dSigma * dScale));
          vfData.push_back((float) dLambda);
          if (bRot)
            vfData.push_back(0.0f);
        }
      if (0 == oF.lNumUsed)
        {
          cout << "ERROR: no usable reflections in " << vsIn[n] << '\n' << flush;
          return (3);
        }
      voFiles.push_back(oF);
    }

  // Orientation of each batch: U = orthogonalized A B^-1, and the
  // wavelengths it predicts for the afout reflections
  if (!voOrient.empty())
    {
      double a3x3dB[3][3], a3x3dBInv[3][3];
      vBmatrix(oCell, a3x3dB);
      (void) dMatInv(a3x3dB, a3x3dBInv);
      bool bMismatch = false;
      for (size_t n = 0; n < voFiles.size(); n++)
        {
          tBatchFile&    rF = voFiles[n];
          const tOrient& rO = voOrient[n];
          double a3x3dM[3][3], a3x3dUB[3][3];
          vMatMul(rO.a3x3dA, a3x3dBInv, a3x3dM);
          rF.dOrthDev = dOrthonormalize(a3x3dM, rF.a3x3dU);
          double dDetU = 0.0, a3x3dTmp[3][3];
          dDetU = dMatInv(rF.a3x3dU, a3x3dTmp);
          if (0.9 > dDetU)
            {
              cout << "ERROR: " << rO.sFile << ": the orientation matrix is not a proper"
                   << " rotation (det " << dDetU << ")\n" << flush;
              return (4);
            }
          vMatMul(rF.a3x3dU, a3x3dB, a3x3dUB);
          double dSum = 0.0, dSumU = 0.0;
          size_t nNum = rF.vdLambda.size();
          for (size_t k = 0; k < nNum; k++)
            {
              double a3dH[3] = {(double) rF.vnHKL[3 * k], (double) rF.vnHKL[3 * k + 1],
                                (double) rF.vnHKL[3 * k + 2]};
              double a3dX[3], a3dXU[3];
              vMatVec(rO.a3x3dA, a3dH, a3dX);
              vMatVec(a3x3dUB, a3dH, a3dXU);
              double dRel  = dLaueLambda(a3dX, rO.a3dS0) / rF.vdLambda[k] - 1.0;
              double dRelU = dLaueLambda(a3dXU, rO.a3dS0) / rF.vdLambda[k] - 1.0;
              dSum  += dRel * dRel;
              dSumU += dRelU * dRelU;
            }
          rF.dLamRms  = sqrt(dSum / (double) std::max((size_t) 1, nNum));
          rF.dLamRmsU = sqrt(dSumU / (double) std::max((size_t) 1, nNum));
          if (0.01 < rF.dLamRms)
            {
              cout << "ERROR: " << rO.sFile << " does not belong to " << rF.sFile
                   << ": the wavelengths it predicts differ by "
                   << 100.0 * rF.dLamRms << " % rms.\n" << flush;
              bMismatch = true;
            }
          else if (0.001 < rF.dLamRms)
            cout << "WARNING: " << rO.sFile << " predicts the wavelengths of " << rF.sFile
                 << " only to " << 100.0 * rF.dLamRms << " % rms.\n" << flush;
        }
      if (bMismatch)
        {
          cout << "       Give the refined header of each run in the order of the afout files.\n"
               << flush;
          return (4);
        }
    }

  try
    {
      double dLambdaSum = 0.0;
      long   lNum = 0;
      for (size_t n = 0; n < voFiles.size(); n++)
        {
          dLambdaSum += voFiles[n].dLambdaSum;
          lNum       += voFiles[n].lNumUsed;
        }

      gemmi::Mtz oMtz(true);
      oMtz.title = ("" != sTitle) ? sTitle : "From dtexport (d*TREK TOF-Laue)";
      oMtz.set_spacegroup(poSG);
      oMtz.set_cell_for_all(oCell);
      gemmi::Mtz::Dataset& roDS = oMtz.add_dataset(sDataset);
      roDS.project_name = sProject;
      roDS.crystal_name = sCrystal;
      roDS.wavelength   = dLambdaSum / (double) lNum;
      int nDS = roDS.id;
      oMtz.add_column("M/ISYM", 'Y', nDS, -1, false);
      oMtz.add_column("BATCH",  'B', nDS, -1, false);
      oMtz.add_column("I",      'J', nDS, -1, false);
      oMtz.add_column("SIGI",   'Q', nDS, -1, false);
      oMtz.add_column("LAMBDA", 'R', nDS, -1, false);
      if (bRot)
        oMtz.add_column("ROT",  'R', nDS, -1, false);
      if ((int) oMtz.columns.size() != nCols)
        throw std::runtime_error("internal error: column count");
      oMtz.set_data(vfData.data(), vfData.size());

      for (size_t n = 0; n < voFiles.size(); n++)
        {
          const tBatchFile& rF = voFiles[n];
          gemmi::Mtz::Batch oBatch;
          oBatch.number = rF.nBatch;
          string sBase  = rF.sFile.substr(rF.sFile.find_last_of('/') + 1);
          oBatch.title  = sBase.substr(0, 70);
          oBatch.set_cell(oCell);
          oBatch.set_dataset_id(nDS);
          oBatch.ints[14] = 3;                          // LDTYPE: 3 = Laue
          oBatch.set_wavelength((float) (rF.dLambdaSum / (double) rF.lNumUsed));
          if (rF.bOrient)
            {
              // MTZ batch header layout as CCP4 MtzArrayToBatch
              const tOrient& rO = voOrient[n];
              oBatch.ints[12] = 1;                      // NCRYST
              oBatch.ints[15] = 1;                      // JSAXS: scan axis E1
              oBatch.ints[17] = 1;                      // NGONAX
              oBatch.ints[19] = 1;                      // NDET
              for (int j = 0; j < 3; j++)
                for (int i = 0; i < 3; i++)
                  oBatch.floats[6 + 3 * j + i] = (float) rF.a3x3dU[i][j];   // UMAT, by column
              oBatch.floats[36] = -0.5f;                // PHISTT
              oBatch.floats[37] =  0.5f;                // PHIEND
              oBatch.floats[47] =  1.0f;                // PHIRANGE
              for (int i = 0; i < 3; i++)
                {
                  oBatch.floats[38 + i] = (float) rO.a3dE1[i];   // SCANAX
                  oBatch.floats[59 + i] = (float) rO.a3dE1[i];   // E1
                  oBatch.floats[80 + i] = (float) rO.a3dS0[i];   // SOURCE
                  oBatch.floats[83 + i] = (float) rO.a3dS0[i];   // SO
                }
              oBatch.axes.push_back(rO.sAxis.substr(0, 8));
            }
          oMtz.batches.push_back(oBatch);
        }

      // Original indices -> asymmetric unit and M/ISYM, then sort as POINTLESS
      oMtz.indices_switched_to_original = true;
      if (!oMtz.switch_to_asu_hkl())
        throw std::runtime_error("could not reduce the indices to the asymmetric unit");
      oMtz.sort(5);
      oMtz.update_reso();

      time_t tNow = time(NULL);
      char a40cDate[40];
      strftime(a40cDate, sizeof(a40cDate), "%Y-%m-%d %H:%M:%S", localtime(&tNow));
      oMtz.history.push_back(string("From dtexport ") + a40cDate + ", "
                             + std::to_string(voFiles.size()) + " afout files"
                             + (voOrient.empty() ? "" : " with orientations"));
      oMtz.write_to_file(sOut);

      char a300cLine[300];
      cout << "\n  Batch  Read      Used   I=0    I<0   Lambda range      Mean  File\n";
      for (size_t n = 0; n < voFiles.size(); n++)
        {
          const tBatchFile& rF = voFiles[n];
          snprintf(a300cLine, sizeof(a300cLine),
                   "  %5d %6ld %8ld %6ld %6ld   %5.3f - %5.3f  %6.3f  ",
                   rF.nBatch, rF.lNumRead, rF.lNumUsed, rF.lNumZero, rF.lNumNeg,
                   rF.a2dLambda[0], rF.a2dLambda[1], rF.dLambdaSum / (double) rF.lNumUsed);
          cout << a300cLine << rF.sFile << '\n';
        }
      if (!voOrient.empty())
        {
          cout << "\n  Batch  Orientation header                  Axis      rms dlambda/lambda (%)"
               << "  Cell misfit\n"
               << "                                                         header UB  U, cell\n";
          for (size_t n = 0; n < voFiles.size(); n++)
            {
              const tBatchFile& rF = voFiles[n];
              string sShown = voOrient[n].sFile;
              if (34 < sShown.size())
                sShown = "..." + sShown.substr(sShown.size() - 31);
              snprintf(a300cLine, sizeof(a300cLine), "  %5d  %-34s  %-8s  %9.4f  %7.4f  %10.4f\n",
                       rF.nBatch, sShown.c_str(), voOrient[n].sAxis.c_str(),
                       100.0 * rF.dLamRms, 100.0 * rF.dLamRmsU, rF.dOrthDev);
              cout << a300cLine;
            }
          cout << "  (Cell misfit: max |M'M - I| of M = G (UB)_dtrek B^-1 before it is made a\n"
               << "  rotation; it measures how far the refined cell is from the MTZ cell.)\n";
        }
      snprintf(a300cLine, sizeof(a300cLine),
               "\n%s: %d reflections, %d batches%s, %s, cell %.4f %.4f %.4f %.3f %.3f %.3f,"
               " resolution %.2f - %.2f A\n",
               sOut.c_str(), oMtz.nreflections, (int) oMtz.batches.size(),
               voOrient.empty() ? "" : " with orientations",
               poSG->xhm().c_str(), oCell.a, oCell.b, oCell.c, oCell.alpha, oCell.beta,
               oCell.gamma, oMtz.resolution_low(), oMtz.resolution_high());
      cout << a300cLine << flush;

      long lNumZero = 0, lNumNeg = 0;
      for (size_t n = 0; n < voFiles.size(); n++)
        {
          lNumZero += voFiles[n].lNumZero;
          lNumNeg  += voFiles[n].lNumNeg;
        }
      if ( (0 == lNumNeg) && (lNumZero > lNum / 20) )
        cout << "WARNING: " << lNumZero << " intensities are exactly 0 and none is negative:\n"
             << "         the afout files were probably written without negative intensities.\n"
             << "         For AIMLESS, write them with dtref2afout -keepneg.\n" << flush;
    }
  catch (const std::exception& roErr)
    {
      cout << "ERROR: " << roErr.what() << '\n' << flush;
      return (4);
    }
  return (0);
}
