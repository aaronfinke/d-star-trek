//
// dtlauenorm.cpp
//    Wavelength normalization and run scaling of TOF-Laue reflection lists
//    from dtintegrate (the LAUENORM step of the original LANSCE PCS
//    workflow).  Writes one normalized, unmerged reflection list (sBatch =
//    run number) for merging with dtscaleaverage, and the normalization
//    curve.  See CLaueNorm.h for the model.
//
//    Usage: dtlauenorm (-spacegroup N | -head file.head) [options]
//                      list1.ref [list2.ref ...]
//

//+Include files

#include <stdlib.h>
#include <vector>
#include "Dtrek.h"
#include "dtrekdefs.h"
#include "dtreksys.h"
#include "Cimage_header.h"
#include "Creflnlist.h"
#include "Cspacegroup.h"
#include "CLaueNorm.h"

#if !defined(VC6)
using std::cout;
using std::flush;
#endif

static void
vUsage(void)
{
  cout << "\nUsage: dtlauenorm (-spacegroup N | -head file.head) [options]\n"
       << "                  list1.ref [list2.ref ...]\n\n"
       << "  Each list (a TOF-Laue dtintegrate output) is one run.  Fits\n"
       << "     I = k_run * N(lambda) * J_hkl\n"
       << "  with N piecewise linear and smooth, and writes I/(k_run N(lambda)).\n\n"
       << "  -spacegroup N      Space group number (or take it from -head).\n"
       << "  -head file         Header with CRYSTAL_SPACEGROUP.\n"
       << "  -wave lmin lmax    Wavelength range (default: that of the data).\n"
       << "  -knots n           Knots of N(lambda) (default 10).\n"
       << "  -smooth s          Smoothness restraint weight (default 0.005).\n"
       << "  -sigma s           Only I > s*sigma(I) are used for the fit (default 1).\n"
       << "  -reject z          Outlier limit in sigmas (default 5, 0: none).\n"
       << "  -errfac e          Error model: var = sigma^2 + (e I)^2 (default 0.05).\n"
       << "  -cycles n          Maximum cycles per rejection pass (default 50).\n"
       << "  -wavefield name    Wavelength field (default fCalwavelength).\n"
       << "  -anom              Keep Friedel mates separate.\n"
       << "  -keeprejects       Also write the rejected outliers.\n"
       << "  -out file          Normalized list (default dtlauenorm.ref).\n"
       << "  -curve file        Normalization curve table (default dtlauenorm_curve.txt).\n"
       << flush;
}

int
main(int argc, char *argv[])
{
  int     nSpacegroup  = 0;
  double  a2dWave[2]   = {0.0, 0.0};
  int     nKnots       = 10;
  double  dSmooth      = 0.005;
  double  dSigma       = 1.0;
  double  dReject      = 5.0;
  double  dErrFac      = 0.05;
  int     nCycles      = 50;
  bool    bAnom        = false;
  bool    bKeepRejects = false;
  Cstring sWaveField   = Creflnlist::ms_sfCalcWavelength;
  Cstring sOut         = "dtlauenorm.ref";
  Cstring sCurve       = "dtlauenorm_curve.txt";
  std::vector<Cstring> vsIn;

  vDtrekSetModuleName("dtlauenorm");
  cout << "\ndtlauenorm: TOF-Laue wavelength normalization\n" << flush;

  for (int i = 1; i < argc; i++)
    {
      Cstring sArg = argv[i];
      bool    bMore1 = (i + 1 < argc);
      if ( ("-spacegroup" == sArg) && bMore1 )
        nSpacegroup = atoi(argv[++i]);
      else if ( ("-head" == sArg) && bMore1 )
        {
          Cstring       sHead = argv[++i];
          Cimage_header oHeader(sHead);
          if (   !oHeader.bIsAvailable()
              || (0 != oHeader.nGetValue(Cstring("CRYSTAL_SPACEGROUP"), &nSpacegroup)) )
            {
              cout << "ERROR: no CRYSTAL_SPACEGROUP in " << argv[i] << '\n' << flush;
              return (1);
            }
        }
      else if ( ("-wave" == sArg) && (i + 2 < argc) )
        {
          a2dWave[0] = atof(argv[++i]);
          a2dWave[1] = atof(argv[++i]);
          if (a2dWave[0] > a2dWave[1])
            std::swap(a2dWave[0], a2dWave[1]);
        }
      else if ( ("-knots" == sArg) && bMore1 )
        nKnots = atoi(argv[++i]);
      else if ( ("-smooth" == sArg) && bMore1 )
        dSmooth = atof(argv[++i]);
      else if ( ("-sigma" == sArg) && bMore1 )
        dSigma = atof(argv[++i]);
      else if ( ("-reject" == sArg) && bMore1 )
        dReject = atof(argv[++i]);
      else if ( ("-errfac" == sArg) && bMore1 )
        dErrFac = atof(argv[++i]);
      else if ( ("-cycles" == sArg) && bMore1 )
        nCycles = atoi(argv[++i]);
      else if ( ("-wavefield" == sArg) && bMore1 )
        sWaveField = argv[++i];
      else if ("-anom" == sArg)
        bAnom = true;
      else if ("-keeprejects" == sArg)
        bKeepRejects = true;
      else if ( ("-out" == sArg) && bMore1 )
        sOut = argv[++i];
      else if ( ("-curve" == sArg) && bMore1 )
        sCurve = argv[++i];
      else if ( ("-help" == sArg) || ("-h" == sArg) )
        {
          vUsage();
          return (0);
        }
      else if ('-' == sArg.GetAt(0))
        {
          cout << "ERROR: unknown or incomplete option " << sArg << '\n' << flush;
          vUsage();
          return (1);
        }
      else
        vsIn.push_back(sArg);
    }

  if (vsIn.empty() || (0 >= nSpacegroup))
    {
      if (0 >= nSpacegroup)
        cout << "ERROR: give -spacegroup N or -head file.head\n" << flush;
      vUsage();
      return (1);
    }

  Cspacegroup oSpace(nSpacegroup);
  CLaueNorm   oNorm;
  oNorm.vSetLambdaRange(a2dWave[0], a2dWave[1]);
  oNorm.vSetNumKnots(nKnots);
  oNorm.vSetSmooth(dSmooth);
  oNorm.vSetSigmaCut(dSigma);
  oNorm.vSetRejectLimit(dReject);
  oNorm.vSetErrorFactor(dErrFac);
  oNorm.vSetMaxCycles(nCycles);

  cout << "\nSpace group " << nSpacegroup << ", " << vsIn.size() << " runs:\n" << flush;
  for (size_t r = 0; r < vsIn.size(); r++)
    {
      Creflnlist oList;
      if (0 != oList.nRead(vsIn[r]))
        {
          cout << "ERROR: cannot read " << vsIn[r] << '\n' << flush;
          return (2);
        }
      int nAdded = oNorm.nAddList(oList, oSpace, (int) r, sWaveField, bAnom);
      if (0 > nAdded)
        return (2);
      printf("  %3d  %-40s %6d observations\n", (int) r + 1, vsIn[r].string(), nAdded);
    }
  fflush(stdout);

  if (0 != oNorm.nFit())
    return (3);
  oNorm.vListResults();

  Creflnlist oOut;
  int nNumOut = oNorm.nMakeOutput(oOut, bKeepRejects);
  if (0 != oOut.nWrite(sOut))
    {
      cout << "ERROR: cannot write " << sOut << '\n' << flush;
      return (4);
    }
  cout << "\n" << nNumOut << " normalized reflections written to " << sOut << '\n';
  if (0 == oNorm.nWriteCurve(sCurve))
    cout << "Normalization curve written to " << sCurve << '\n';
  cout << flush;
  return (0);
}
