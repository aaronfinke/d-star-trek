//
// dtref2afout.cpp
//    Convert d*TREK TOF-Laue reflection lists (dtintegrate output) to the
//    "afout" format read by LAUENORM, as the awk script negative.awk of the
//    original LANSCE PCS workflow did:
//
//       h k l lambda 0 scale*I scale*sigma     (%5d%5d%5d%10.5f%10.5f%10d%10d)
//
//    with negative intensities written as 0.  Unlike negative.awk, which
//    used fixed column numbers, the columns are found by field name, so the
//    lists of this d*TREK (which have an embedded header) work as well.
//    The numbers are converted from the text of the .ref file exactly as
//    awk does, so the output for a legacy list is identical to the legacy
//    .afout file.
//
//    Usage: dtref2afout [-scale f|auto] [-keepneg] [-wavefield name] [-out file]
//                       file1.ref [file2.ref ...]
//
//    The intensities of this d*TREK are about 1000 times smaller than those
//    of the original PCS d*TREK (whose -type 204 scale is not known), so the
//    default scale of 10 of negative.awk truncates them to a few digits; use
//    -scale auto (or e.g. -scale 10000) for them.  LAUENORM rescales its
//    input anyway, so only the number of digits matters.
//

//+Include files

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

using std::cout;
using std::flush;
using std::string;
using std::vector;

static void
vUsage(void)
{
  cout << "\nUsage: dtref2afout [-scale f|auto] [-keepneg] [-wavefield name] [-out file]\n"
       << "                   file1.ref [file2.ref ...]\n\n"
       << "  Writes the LAUENORM input file of each TOF-Laue reflection list:\n"
       << "     h k l lambda 0 scale*I scale*sigma\n"
       << "  as the negative.awk script of the original PCS workflow did.\n\n"
       << "  -scale f        Multiplier of I and sigma before they are truncated to\n"
       << "                  integers (default 10, as negative.awk).  auto: the power\n"
       << "                  of 10 that makes the largest scaled intensity 1e6-1e7\n"
       << "                  (use it for lists of this d*TREK, whose intensities are\n"
       << "                  about 1000x smaller than those of the original one).\n"
       << "  -keepneg        Write negative intensities (default: write them as 0).\n"
       << "  -wavefield name Field with the wavelength (default fCalwavelength).\n"
       << "  -out file       Output file; only with a single input list.  The default\n"
       << "                  is the input file name with .ref replaced by .afout, in\n"
       << "                  the current directory.\n"
       << flush;
}

static int
nFindField(const vector<string>& rvsNames, const string& rsName)
{
  for (size_t i = 0; i < rvsNames.size(); i++)
    if (rvsNames[i] == rsName)
      return ((int) i);
  return (-1);
}

// awk's %d of a non-integer number truncates towards zero
static long long
llAwkInt(const double dValue)
{
  return ((long long) dValue);
}

static int
nConvert(const string& rsIn, const string& rsOut, const double dScale,
         const bool bKeepNeg, const string& rsWaveField)
{
  std::ifstream oIn(rsIn.c_str());
  if (!oIn)
    {
      cout << "ERROR: cannot open " << rsIn << '\n' << flush;
      return (1);
    }

  // Header: "nInts nFloats nStrings [nHeaderLines]", the embedded header
  // lines (if any), then one field name per line.

  string sLine;
  if (!std::getline(oIn, sLine))
    {
      cout << "ERROR: " << rsIn << " is empty\n" << flush;
      return (2);
    }
  std::istringstream oCounts(sLine);
  int a4nCounts[4] = {0, 0, 0, 0};
  int nNumCounts   = 0;
  while ( (nNumCounts < 4) && (oCounts >> a4nCounts[nNumCounts]) )
    nNumCounts++;
  if (3 > nNumCounts)
    {
      cout << "ERROR: " << rsIn << " is not a d*TREK reflection list\n" << flush;
      return (2);
    }
  for (int i = 0; i < a4nCounts[3]; i++)
    std::getline(oIn, sLine);

  int nNumFields = a4nCounts[0] + a4nCounts[1] + a4nCounts[2];
  vector<string> vsNames;
  for (int i = 0; i < nNumFields; i++)
    {
      if (!std::getline(oIn, sLine))
        {
          cout << "ERROR: " << rsIn << ": truncated header\n" << flush;
          return (2);
        }
      std::istringstream oName(sLine);
      string sName;
      oName >> sName;
      vsNames.push_back(sName);
    }

  int nFH = nFindField(vsNames, "nH");
  int nFK = nFindField(vsNames, "nK");
  int nFL = nFindField(vsNames, "nL");
  int nFI = nFindField(vsNames, "fIntensity");
  int nFS = nFindField(vsNames, "fSigmaI");
  int nFW = nFindField(vsNames, rsWaveField);
  if ( (0 > nFH) || (0 > nFK) || (0 > nFL) || (0 > nFI) || (0 > nFS) )
    {
      cout << "ERROR: " << rsIn << " lacks nH, nK, nL, fIntensity or fSigmaI\n" << flush;
      return (3);
    }
  if (0 > nFW)
    {
      cout << "ERROR: " << rsIn << " has no field " << rsWaveField
           << " (not a TOF-Laue dtintegrate list?)\n" << flush;
      return (3);
    }

  FILE *pFOut = fopen(rsOut.c_str(), "w");
  if (NULL == pFOut)
    {
      cout << "ERROR: cannot write " << rsOut << '\n' << flush;
      return (4);
    }

  vector< vector<string> > vvsRecs;
  vector<string> vsTok;
  double dMaxI = 0.0;
  while (std::getline(oIn, sLine))
    {
      std::istringstream oTok(sLine);
      string sTok;
      vsTok.clear();
      while (oTok >> sTok)
        vsTok.push_back(sTok);
      if ((int) vsTok.size() < nNumFields)
        continue;       // blank or incomplete line
      double dI = fabs(strtod(vsTok[nFI].c_str(), NULL));
      if (dI > dMaxI)
        dMaxI = dI;
      vvsRecs.push_back(vsTok);
    }

  double dUseScale = dScale;
  if (0.0 >= dUseScale)                      // -scale auto
    {
      dUseScale = 1.0;
      if (0.0 < dMaxI)
        dUseScale = pow(10.0, floor(log10(9.99e6 / dMaxI)));
    }

  int nNumOut  = 0;
  int nNumTiny = 0;
  for (size_t n = 0; n < vvsRecs.size(); n++)
    {
      const vector<string>& rvsTok = vvsRecs[n];
      double dI = strtod(rvsTok[nFI].c_str(), NULL);
      double dS = strtod(rvsTok[nFS].c_str(), NULL);
      long long llI = 0;
      if ( (0.0 < dI) || bKeepNeg )
        llI = llAwkInt(dI * dUseScale);
      if (100 > llAwkInt(dS * dUseScale))
        nNumTiny++;
      fprintf(pFOut, "%5ld%5ld%5ld%10.5f%10.5f%10lld%10lld\n",
              strtol(rvsTok[nFH].c_str(), NULL, 10),
              strtol(rvsTok[nFK].c_str(), NULL, 10),
              strtol(rvsTok[nFL].c_str(), NULL, 10),
              strtod(rvsTok[nFW].c_str(), NULL), 0.0,
              llI, llAwkInt(dS * dUseScale));
      nNumOut++;
    }
  fclose(pFOut);
  cout << rsIn << " -> " << rsOut << ": " << nNumOut << " reflections, scale "
       << dUseScale << '\n' << flush;
  if (nNumTiny > nNumOut / 2)
    cout << "WARNING: most scaled sigmas are below 100, so the integer intensities\n"
         << "         keep few digits.  Use -scale auto.\n" << flush;
  return (0);
}

int
main(int argc, char *argv[])
{
  double         dScale     = 10.0;
  bool           bKeepNeg   = false;
  string         sWaveField = "fCalwavelength";
  string         sOut       = "";
  vector<string> vsIn;

  for (int i = 1; i < argc; i++)
    {
      string sArg = argv[i];
      if ( ("-scale" == sArg) && (i + 1 < argc) )
        {
          i++;
          dScale = ("auto" == string(argv[i])) ? 0.0 : atof(argv[i]);
          if ( (0.0 >= dScale) && ("auto" != string(argv[i])) )
            {
              cout << "ERROR: -scale needs a positive number or auto\n" << flush;
              return (1);
            }
        }
      else if ("-keepneg" == sArg)
        bKeepNeg = true;
      else if ( ("-wavefield" == sArg) && (i + 1 < argc) )
        sWaveField = argv[++i];
      else if ( ("-out" == sArg) && (i + 1 < argc) )
        sOut = argv[++i];
      else if ( ("-help" == sArg) || ("-h" == sArg) )
        {
          vUsage();
          return (0);
        }
      else if ('-' == sArg[0])
        {
          cout << "ERROR: unknown option " << sArg << '\n' << flush;
          vUsage();
          return (1);
        }
      else
        vsIn.push_back(sArg);
    }

  if (vsIn.empty())
    {
      vUsage();
      return (1);
    }
  if ( ("" != sOut) && (1 != vsIn.size()) )
    {
      cout << "ERROR: -out needs exactly one input list\n" << flush;
      return (1);
    }

  int nStat = 0;
  for (size_t i = 0; i < vsIn.size(); i++)
    {
      string sThisOut = sOut;
      if ("" == sThisOut)
        {
          string sBase = vsIn[i];
          size_t nSlash = sBase.find_last_of('/');
          if (string::npos != nSlash)
            sBase = sBase.substr(nSlash + 1);
          size_t nLen = sBase.length();
          if ( (4 < nLen) && (".ref" == sBase.substr(nLen - 4)) )
            sBase = sBase.substr(0, nLen - 4);
          sThisOut = sBase + ".afout";
        }
      if (0 != nConvert(vsIn[i], sThisOut, dScale, bKeepNeg, sWaveField))
        nStat = 1;
    }
  return (nStat);
}
