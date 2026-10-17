//
// dtpcsproject.cpp
//    Write one TOF channel, or the sum or maximum over a range of TOF
//    channels, of a LANSCE PCS NeXus file (HDF5) as a d*TREK image.
//    The maximum image is the "overlay" of the PCS paper (Langan & Greene,
//    J. Appl. Cryst. 2004, 37, 253-257), useful for checking found spots
//    and predicted reflections.
//
//    Usage: dtpcsproject in_file.h5 out_image.img
//                        [-chan N | -sum [c0 c1] | -max [c0 c1]]
//

//+Include files

#include <stdlib.h>
#include "Dtrek.h"
#include "dtrekdefs.h"
#include "dtreksys.h"
#include "Cimage.h"
#include "nexuspcs.h"

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif

static void
vUsage(void)
{
  cout << "\nUsage: dtpcsproject in_file.h5 out_image.img\n"
       << "                    [-chan N | -sum [c0 c1] | -max [c0 c1]]\n\n"
       << "  -chan N       Write TOF channel N (0-based).\n"
       << "  -sum [c0 c1]  Write the sum over channels c0 to c1 (default: all).\n"
       << "                This is a conventional Laue pattern.  Default mode.\n"
       << "  -max [c0 c1]  Write the maximum over channels c0 to c1 (default: all),\n"
       << "                the \"overlay\" image.\n"
       << flush;
}

int main(int   argc,      // Command line argument count
         char *argv[])    // Pointers to command line args
{
  int     nStat;
  int     nRequested;
  int     nChan0, nChan1;
  bool    bRange;
  Cstring sIn, sOut, sOpt;
  Cimage  oImage;

  vDtrekSetModuleName("dtpcsproject");
  cout << "\ndtpcsproject: LANSCE PCS TOF NeXus channel/sum/maximum images\n";
  cout << D_K_DTREKVersion << endl << flush;
  cout << "Command line:\n " << sGetCommandLine(argc, argv, 71) << endl << flush;

  argc--; argv++;
  if (2 > argc)
    {
      cout << "ERROR - not enough arguments!\n";
      vUsage();
      return (1);
    }
  sIn  = (const char*) argv[0];
  sOut = (const char*) argv[1];

  nRequested = nPCS_CHANNEL_SUM;
  nChan0     = 0;
  nChan1     = -1;
  bRange     = false;
  for (int i = 2; i < argc; i++)
    {
      sOpt = (const char*) argv[i];
      if ( ("-chan" == sOpt) && (i+1 < argc) )
        {
          nRequested = atoi(argv[++i]);
        }
      else if ( ("-sum" == sOpt) || ("-max" == sOpt) )
        {
          nRequested = ("-sum" == sOpt) ? nPCS_CHANNEL_SUM : nPCS_CHANNEL_MAX;
          if ( (i+2 < argc) && Cstring(argv[i+1]).bIsNumeric()
                            && Cstring(argv[i+2]).bIsNumeric() )
            {
              nChan0 = atoi(argv[++i]);
              nChan1 = atoi(argv[++i]);
              bRange = true;
            }
        }
      else
        {
          cout << "ERROR - unrecognized option " << sOpt << "\n";
          vUsage();
          return (1);
        }
    }

  nStat = oImage.nRead(sIn, nRequested);
  if ( (0 != nStat) || !oImage.bIsAvailable() )
    {
      cout << "ERROR - could not read " << sIn << " as a PCS NeXus file!\n" << flush;
      return (2);
    }

  if (bRange)
    {
      CPCSVolume *poVol = CPCSVolume::poGetInstance();
      if (nPCS_CHANNEL_SUM == nRequested)
        nStat = poVol->nGetSum(nChan0, nChan1, (float*) oImage.pvGetData());
      else
        nStat = poVol->nGetMax(nChan0, nChan1, (float*) oImage.pvGetData());
      oImage.m_oHeader.nReplaceValue(Cstring("PCS_CHANNEL_RANGE"),
                                     Cstring(nChan0) + ' ' + Cstring(nChan1));
    }

  if (0 == nStat)
    nStat = oImage.nWrite(sOut);
  if (0 != nStat)
    {
      cout << "ERROR - writing " << sOut << "!\n" << flush;
      return (3);
    }
  cout << "Wrote " << sOut << "\n" << flush;
  return (0);
}
