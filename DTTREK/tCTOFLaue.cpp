//
// tCTOFLaue.cpp
//    Tests the TOF-Laue physics of CTOFLaue and Cdetector:
//      1. channel -> wavelength -> channel round trip;
//      2. the TOF of integer channels equals the bin edges of the header;
//      3. given a header and a spot list of the original PCS d*TREK (dtfind
//         output with fObs_rot_mid and fDetwavelength), the observed
//         wavelengths are reproduced with the flight path to the detector
//         cylinder (as the original dtfind did); the difference that a
//         flight path to a flat detector (L2 = D / (S.N)) would make is
//         reported for comparison.
//
//    Usage: tCTOFLaue file.head [legacy_find.ref]
//    Exit status 0 if all tests pass.
//

//+Include files

#include <math.h>
#include <stdio.h>
#include <vector>
#include "Dtrek.h"
#include "dtrekvec.h"
#include "Cimage_header.h"
#include "Cdetector.h"
#include "Creflnlist.h"
#include "CTOFLaue.h"

#if !defined(VC6)
using std::cout;
using std::flush;
#endif

int
main(int argc, char *argv[])
{
  if (2 > argc)
    {
      cout << "Usage: tCTOFLaue file.head [legacy_find.ref]\n" << flush;
      return (1);
    }

  int     nFail = 0;
  Cstring sHead = argv[1];
  Cimage_header oHeader(sHead);
  if (!oHeader.bIsAvailable())
    {
      cout << "FAIL: cannot read " << sHead << '\n' << flush;
      return (1);
    }
  Cstring sDetName = "Detector1_";
  (void) oHeader.nGetValue(Cstring("DETECTOR_NAMES"), &sDetName);

  CTOFLaue oTOF(oHeader, sDetName);
  if (!oTOF.bIsAvailable())
    {
      cout << "FAIL: the header does not describe TOF-Laue data\n" << flush;
      return (1);
    }
  int nNumChan = oTOF.nGetNumChannels();

  // 1. Round trip

  double dMaxErr = 0.0;
  for (double dL2 = 700.0; dL2 <= 820.0; dL2 += 20.0)
    for (double dChan = 0.0; dChan <= (double) nNumChan; dChan += 0.125)
      {
        double dLambda = oTOF.dChanToLambda(dChan, dL2);
        dMaxErr = std::max(dMaxErr, fabs(oTOF.dLambdaToChan(dLambda, dL2) - dChan));
      }
  printf("1. channel -> lambda -> channel: max error %.2e channel  %s\n",
         dMaxErr, (1.0e-6 > dMaxErr) ? "PASS" : "FAIL");
  if (1.0e-6 <= dMaxErr)
    nFail++;

  // 2. Bin edges

  std::vector<double> vdEdges(nNumChan + 1, 0.0);
  double dBin = 1.0;
  (void) oHeader.nGetValue(sDetName + "ROTATION_BIN", &dBin);
  if (0 != oHeader.nGetValue(sDetName + "ROTATION_ARRAY", nNumChan + 1, &vdEdges[0]))
    {
      printf("2. bin edges: FAIL (cannot read %sROTATION_ARRAY)\n", sDetName.string());
      nFail++;
    }
  else
    {
      dMaxErr = 0.0;
      for (int c = 0; c <= nNumChan; c++)
        dMaxErr = std::max(dMaxErr, fabs(oTOF.dChanToTOF((double) c) - vdEdges[c] * dBin));
      printf("2. TOF of integer channels vs bin edges: max error %.2e us  %s\n",
             dMaxErr, (1.0e-6 > dMaxErr) ? "PASS" : "FAIL");
      if (1.0e-6 <= dMaxErr)
        nFail++;
    }

  // 3. Legacy observed wavelengths

  if (3 <= argc)
    {
      Creflnlist oList;
      Cstring    sRef = argv[2];
      if (0 != oList.nRead(sRef))
        {
          cout << "FAIL: cannot read " << sRef << '\n' << flush;
          return (1);
        }
      int nFIWave = oList.nGetFieldIndex(Creflnlist::ms_sfDetWavelength);
      int nFIRot  = oList.m_nFI_fObsRotMid;
      if ( (0 > nFIWave) || (0 > nFIRot) )
        {
          cout << "FAIL: " << sRef << " lacks fDetwavelength or fObs_rot_mid\n" << flush;
          return (1);
        }
      Cdetector oDet(oHeader, sDetName);
      double a3x3dDD[3][3], a3dDN[3];
      (void) oDet.nCalcGetDDDN(&a3x3dDD[0][0], a3dDN);
      (void) fNormVec3D(a3dDN);
      double dDist = fabs(oDet.fGetDistance());

      double dSumFlat = 0.0, dSumCyl = 0.0, dMaxCyl = 0.0, dSumRel = 0.0;
      int    nNum = 0;
      for (int i = 0; i < oList.nGetNumReflns(); i++)
        {
          Crefln *poRefln = oList.poGetRefln(i);
          double dLambdaLeg = poRefln->fGetField(nFIWave);
          double dChan      = poRefln->fGetField(nFIRot);
          double dLambdaCyl, a3dS[3], dL2Cyl;
          float  a3fMM[3];
          if (0 != oDet.nCalcTOFLaueObs(oTOF, poRefln->fGetField(oList.m_nFI_fObsPx0),
                                        poRefln->fGetField(oList.m_nFI_fObsPx1), dChan,
                                        &dLambdaCyl, a3fMM, a3dS, &dL2Cyl))
            continue;
          double dCos = fabs(fDot3D(a3dS, a3dDN));
          if (1.0e-3 > dCos)
            continue;
          double dLambdaFlat = oTOF.dChanToLambda(dChan, dDist / dCos);
          dSumFlat += (dLambdaFlat - dLambdaLeg) * (dLambdaFlat - dLambdaLeg);
          dSumCyl  += (dLambdaCyl  - dLambdaLeg) * (dLambdaCyl  - dLambdaLeg);
          dSumRel  += (dLambdaCyl - dLambdaLeg) / dLambdaLeg;
          dMaxCyl   = std::max(dMaxCyl, fabs(dLambdaCyl - dLambdaLeg) / dLambdaLeg);
          nNum++;
        }
      if (0 == nNum)
        {
          printf("3. legacy wavelengths: FAIL (no usable spots)\n");
          nFail++;
        }
      else
        {
          double dRmsCyl = sqrt(dSumCyl / nNum);
          printf("3. legacy fDetwavelength of %d spots, cylinder flight path:\n"
                 "   rms difference %.5f A, mean %+.3f%%, max %.3f%%  %s\n", nNum, dRmsCyl,
                 100.0 * dSumRel / nNum, 100.0 * dMaxCyl,
                 (0.001 > dRmsCyl) ? "PASS" : "FAIL");
          printf("   (a flat-detector flight path would give rms %.5f A)\n",
                 sqrt(dSumFlat / nNum));
          if (0.001 <= dRmsCyl)
            nFail++;
        }
    }

  printf("%s\n", (0 == nFail) ? "All tests passed." : "Some tests FAILED.");
  return ( (0 == nFail) ? 0 : 1 );
}
