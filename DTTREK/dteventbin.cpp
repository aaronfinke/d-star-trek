//
// dteventbin.cpp
//    Histogram the neutron events of a NeXus event file (e.g. SNS MANDI)
//    into TOF-Laue panel histograms for d*TREK, written as an NXlauetof
//    HDF5 file.  See nxevents.h.
//
//    Usage: dteventbin events.nxs.h5 [options]
//

//+Include files

#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "nxevents.h"

using std::cout;
using std::flush;
using std::string;
using std::vector;

static void
vUsage(void)
{
  cout << "\nUsage: dteventbin events.nxs.h5 [options]\n\n"
       << "  Histograms the events of each rectangular detector panel into\n"
       << "  (y, x, TOF) arrays and writes them as an NXlauetof file for d*TREK.\n\n"
       << "  -out file         Output file (default: <run>_lauetof.h5 in the\n"
       << "                    current directory).\n"
       << "  -detcal file      ISAW/Mantid .DetCal calibration: L1, T0 and panel\n"
       << "                    positions (default: the instrument definition in\n"
       << "                    the event file).\n"
       << "  -goniocal file    Mantid goniometer calibration (chi-offset,\n"
       << "                    goniometer-tilt), recorded in the output only.\n"
       << "  -tof tmin tmax    TOF range in microseconds (default: that of the events).\n"
       << "  -logbin r         Logarithmic TOF bins, t[i+1]/t[i] = 1 + r (default 0.005).\n"
       << "  -linbin dt        Linear TOF bins of dt microseconds instead.\n"
       << "  -bin n            Sum n x n pixels (default 1).\n"
       << "  -gonio n1 n2 n3   DAS log names of the goniometer angles, outermost\n"
       << "                    first (default omega chi phi, axes y, z, y).\n"
       << "  -banks b1,b2,...  Only these bank numbers.\n"
       << flush;
}

int
main(int argc, char *argv[])
{
  if (2 > argc)
    {
      vUsage();
      return (1);
    }

  string sIn = argv[1];
  if ( ("-help" == sIn) || ("-h" == sIn) )
    {
      vUsage();
      return (0);
    }
  string sOut, sDetCal, sGonioCal;
  double a2dTOF[2] = {0.0, 0.0};
  double dLog = 0.005, dLin = 0.0;
  int    nBin = 1;
  vector<string> vsGonio;
  vector<int>    vnBanks;

  for (int i = 2; i < argc; i++)
    {
      string sArg = argv[i];
      bool   bMore = (i + 1 < argc);
      if ( ("-out" == sArg) && bMore )
        sOut = argv[++i];
      else if ( ("-detcal" == sArg) && bMore )
        sDetCal = argv[++i];
      else if ( ("-goniocal" == sArg) && bMore )
        sGonioCal = argv[++i];
      else if ( ("-tof" == sArg) && (i + 2 < argc) )
        {
          a2dTOF[0] = atof(argv[++i]);
          a2dTOF[1] = atof(argv[++i]);
        }
      else if ( ("-logbin" == sArg) && bMore )
        {
          dLog = atof(argv[++i]);
          dLin = 0.0;
        }
      else if ( ("-linbin" == sArg) && bMore )
        {
          dLin = atof(argv[++i]);
          dLog = 0.0;
        }
      else if ( ("-bin" == sArg) && bMore )
        nBin = atoi(argv[++i]);
      else if ( ("-gonio" == sArg) && bMore )
        {
          while ( (i + 1 < argc) && ('-' != argv[i + 1][0]) )
            vsGonio.push_back(argv[++i]);
        }
      else if ( ("-banks" == sArg) && bMore )
        {
          string sList = argv[++i];
          for (size_t k = 0; k < sList.size(); k++)
            if (',' == sList[k])
              sList[k] = ' ';
          std::istringstream oList(sList);
          int n;
          while (oList >> n)
            vnBanks.push_back(n);
        }
      else
        {
          cout << "ERROR: unknown or incomplete option " << sArg << '\n' << flush;
          vUsage();
          return (1);
        }
    }

  if ("" == sOut)
    {
      string sBase = sIn.substr(sIn.find_last_of('/') + 1);
      const char *apcSuffix[3] = {".nxs.h5", ".h5", ".nxs"};
      for (int k = 0; k < 3; k++)
        {
          size_t nLen = strlen(apcSuffix[k]);
          if ( (sBase.size() > nLen) && (0 == sBase.compare(sBase.size() - nLen, nLen, apcSuffix[k])) )
            {
              sBase = sBase.substr(0, sBase.size() - nLen);
              break;
            }
        }
      sOut = sBase + "_lauetof.h5";
    }

  CNXEventBinner oBinner;
  if (!vsGonio.empty())
    {
      // Axes for the given names: omega/phi about y, chi about z (NeXus frame)
      vector<double> vdAxes;
      for (size_t k = 0; k < vsGonio.size(); k++)
        {
          bool bChi = (string::npos != vsGonio[k].find("chi"));
          vdAxes.push_back(0.0);
          vdAxes.push_back(bChi ? 0.0 : 1.0);
          vdAxes.push_back(bChi ? 1.0 : 0.0);
        }
      oBinner.vSetGonio(vsGonio, vdAxes);
    }
  if (0 != oBinner.nOpen(sIn))
    return (2);
  if ( ("" != sDetCal) && (0 != oBinner.nApplyDetCal(sDetCal)) )
    return (3);
  if ( ("" != sGonioCal) && (0 != oBinner.nReadGonioCal(sGonioCal)) )
    return (3);
  if (0.0 < a2dTOF[1])
    oBinner.vSetTOFRange(a2dTOF[0], a2dTOF[1]);
  if (0.0 < dLin)
    oBinner.vSetLinearBins(dLin);
  else
    oBinner.vSetLogBins(dLog);
  oBinner.vSetPixelBinning(nBin);
  oBinner.vSetBanks(vnBanks);
  oBinner.vList();

  if (0 != oBinner.nBin())
    return (4);
  if (0 != oBinner.nWrite(sOut))
    return (5);
  cout << "Histograms written to " << sOut << '\n' << flush;
  return (0);
}
