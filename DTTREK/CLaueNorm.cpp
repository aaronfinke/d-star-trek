//
// CLaueNorm.cpp
//    Wavelength normalization and run scaling of TOF-Laue intensities.
//    See CLaueNorm.h.
//

//+Include files

#include <math.h>
#include <stdio.h>
#include <algorithm>
#include "CLaueNorm.h"
#include "Creflnlist.h"
#include "Cspacegroup.h"
#include "dtsvd.h"

#if !defined(VC6)
using std::cout;
using std::flush;
#endif

//+Code begin

CLaueNorm::CLaueNorm()
{
  m_a2dLambda[0] = 0.0;
  m_a2dLambda[1] = 0.0;      // 0, 0: use the range of the data
  m_nKnots       = 10;
  m_dSmooth      = 0.005;
  m_dSigmaCut    = 1.0;
  m_dReject      = 5.0;
  m_dErrFac      = 0.05;
  m_nMaxCycles   = 50;
  m_nCyclesDone  = 0;
  m_nNumRejected = 0;
  m_dSigmaScale  = 1.0;
}

CLaueNorm::~CLaueNorm()
{
}

int
CLaueNorm::nAddList(Creflnlist& roList, Cspacegroup& roSpace, const int nRun,
                    const Cstring& rsWaveField, const bool bAnom)
{
  int nFIWave = roList.nGetFieldIndex(rsWaveField);
  if (0 > nFIWave)
    nFIWave = roList.nGetFieldIndex(Creflnlist::ms_sfDetWavelength);
  if (0 > nFIWave)
    {
      cout << "ERROR in CLaueNorm::nAddList: no wavelength field "
           << rsWaveField << " in the reflection list.\n" << flush;
      return (-1);
    }
  if (0 != roList.nReduce(roSpace, bAnom ? 1 : 0))
    return (-2);

  if ((int) m_vdRunScale.size() <= nRun)
    m_vdRunScale.resize(nRun + 1, 1.0);

  int nNum = roList.nGetNumReflns();
  int nAdded = 0;
  for (int i = 0; i < nNum; i++)
    {
      Crefln *poRefln = roList.poGetRefln(i);
      tObs oObs;
      int  nPacked = poRefln->nGetField(roList.m_nFI_nPackedHKL);
      std::map<int, int>::iterator oIt = m_mPackedToUnique.find(nPacked);
      if (oIt == m_mPackedToUnique.end())
        {
          oObs.nUnique = (int) m_vdJ.size();
          m_mPackedToUnique[nPacked] = oObs.nUnique;
          m_vdJ.push_back(0.0);
          m_vnUsedCount.push_back(0);
        }
      else
        oObs.nUnique = oIt->second;
      oObs.nRun      = nRun;
      oObs.a3nHKL[0] = poRefln->nGetH();
      oObs.a3nHKL[1] = poRefln->nGetK();
      oObs.a3nHKL[2] = poRefln->nGetL();
      oObs.fLambda   = poRefln->fGetField(nFIWave);
      oObs.fI        = poRefln->fGetIntensity();
      oObs.fSig      = poRefln->fGetSigmaI();
      oObs.bUse      = false;
      oObs.bReject   = false;
      if ( (0.0f >= oObs.fLambda) || (0.0f >= oObs.fSig) )
        continue;
      m_vObs.push_back(oObs);
      nAdded++;
    }
  return (nAdded);
}

void
CLaueNorm::vBasis(const double dLambda, int *pnKnot, double *pdFrac) const
{
  // The knots are equally spaced in ln(lambda), like the TOF channels:
  // N changes fastest at short wavelengths.

  double dL0   = log(m_a2dLambda[0]);
  double dStep = (log(m_a2dLambda[1]) - dL0) / (double) (m_nKnots - 1);
  double dT    = (log(std::max(dLambda, 1.0e-3)) - dL0) / dStep;
  int    nK    = (int) floor(dT);
  nK = std::max(0, std::min(m_nKnots - 2, nK));
  *pnKnot = nK;
  *pdFrac = std::max(0.0, std::min(1.0, dT - (double) nK));
}

double
CLaueNorm::dNormAt(const double dLambda) const
{
  int    nK;
  double dF;
  vBasis(dLambda, &nK, &dF);
  return ((1.0 - dF) * m_vdKnot[nK] + dF * m_vdKnot[nK + 1]);
}

double
CLaueNorm::dGetNorm(const double dLambda) const
{
  if ((int) m_vdKnot.size() != m_nKnots)
    return (1.0);
  return (dNormAt(dLambda));
}

double
CLaueNorm::dWeight(const tObs& roObs) const
{
  double dVar = (double) roObs.fSig * roObs.fSig
                + m_dErrFac * m_dErrFac * (double) roObs.fI * roObs.fI;
  return ( (0.0 < dVar) ? 1.0 / dVar : 0.0 );
}

void
CLaueNorm::vUpdateJ(void)
{
  // Least-squares J_u for the model I = k N J, from the used observations.

  std::vector<double> vdNum(m_vdJ.size(), 0.0), vdDen(m_vdJ.size(), 0.0);
  std::fill(m_vnUsedCount.begin(), m_vnUsedCount.end(), 0);
  for (size_t j = 0; j < m_vObs.size(); j++)
    {
      const tObs& roObs = m_vObs[j];
      if (!roObs.bUse || roObs.bReject)
        continue;
      double dA = m_vdRunScale[roObs.nRun] * dNormAt(roObs.fLambda);
      double dW = dWeight(roObs);
      vdNum[roObs.nUnique] += dW * dA * roObs.fI;
      vdDen[roObs.nUnique] += dW * dA * dA;
      m_vnUsedCount[roObs.nUnique]++;
    }
  for (size_t u = 0; u < m_vdJ.size(); u++)
    m_vdJ[u] = (0.0 < vdDen[u]) ? vdNum[u] / vdDen[u] : 0.0;
}

int
CLaueNorm::nUpdateKnots(void)
{
  // Linear least squares for the knot values, with a second-difference
  // smoothness restraint.  Only multiply measured reflections carry
  // information on N.

  int nM = m_nKnots;
  std::vector<double> vdA(nM * nM, 0.0), vdB(nM, 0.0), vdX(nM, 0.0);
  for (size_t j = 0; j < m_vObs.size(); j++)
    {
      const tObs& roObs = m_vObs[j];
      if (!roObs.bUse || roObs.bReject || (2 > m_vnUsedCount[roObs.nUnique]))
        continue;
      int    nK;
      double dF;
      vBasis(roObs.fLambda, &nK, &dF);
      double dC  = m_vdRunScale[roObs.nRun] * m_vdJ[roObs.nUnique];
      double a2dD[2] = { (1.0 - dF) * dC, dF * dC };
      double dW  = dWeight(roObs);
      for (int p = 0; p < 2; p++)
        {
          vdB[nK + p] += dW * a2dD[p] * roObs.fI;
          for (int q = 0; q < 2; q++)
            vdA[(nK + p) * nM + nK + q] += dW * a2dD[p] * a2dD[q];
        }
    }

  double dTrace = 0.0;
  for (int m = 0; m < nM; m++)
    dTrace += vdA[m * nM + m];
  if (0.0 >= dTrace)
    return (1);

  // Knots with no data are tied to their neighbours by the restraint only;
  // the restraint weight is relative to the mean diagonal element.
  double dMu = std::max(m_dSmooth, 1.0e-6) * dTrace / (double) nM;
  for (int m = 1; m < nM - 1; m++)
    {
      double a3dD[3] = {1.0, -2.0, 1.0};
      for (int p = 0; p < 3; p++)
        for (int q = 0; q < 3; q++)
          vdA[(m - 1 + p) * nM + m - 1 + q] += dMu * a3dD[p] * a3dD[q];
    }

  if (0 != nSolveAxB_svd(nM, &vdA[0], &vdX[0], &vdB[0]))
    return (2);

  double dMax = 0.0;
  for (int m = 0; m < nM; m++)
    dMax = std::max(dMax, vdX[m]);
  if (0.0 >= dMax)
    return (3);
  for (int m = 0; m < nM; m++)
    m_vdKnot[m] = std::max(vdX[m], 1.0e-3 * dMax);
  return (0);
}

void
CLaueNorm::vUpdateRunScales(void)
{
  std::vector<double> vdNum(m_vdRunScale.size(), 0.0), vdDen(m_vdRunScale.size(), 0.0);
  for (size_t j = 0; j < m_vObs.size(); j++)
    {
      const tObs& roObs = m_vObs[j];
      if (!roObs.bUse || roObs.bReject || (2 > m_vnUsedCount[roObs.nUnique]))
        continue;
      double dA = dNormAt(roObs.fLambda) * m_vdJ[roObs.nUnique];
      double dW = dWeight(roObs);
      vdNum[roObs.nRun] += dW * dA * roObs.fI;
      vdDen[roObs.nRun] += dW * dA * dA;
    }
  for (size_t r = 0; r < m_vdRunScale.size(); r++)
    if ( (0.0 < vdDen[r]) && (0.0 < vdNum[r]) )
      m_vdRunScale[r] = vdNum[r] / vdDen[r];
}

void
CLaueNorm::vNormalize(void)
{
  // Mean run scale 1 and maximum of N 1 (J takes up the rest).

  double dSum = 0.0;
  for (size_t r = 0; r < m_vdRunScale.size(); r++)
    dSum += m_vdRunScale[r];
  if (0.0 < dSum)
    {
      double dF = (double) m_vdRunScale.size() / dSum;
      for (size_t r = 0; r < m_vdRunScale.size(); r++)
        m_vdRunScale[r] *= dF;
    }
  double dMax = *std::max_element(m_vdKnot.begin(), m_vdKnot.end());
  if (0.0 < dMax)
    for (size_t m = 0; m < m_vdKnot.size(); m++)
      m_vdKnot[m] /= dMax;
}

int
CLaueNorm::nRejectOutliers(void)
{
  // Among reflections measured at least three times, reject the worst
  // observation of a reflection if it is more than m_dReject sigmas from
  // the mean of the others.  The sigmas are rescaled by the robust spread
  // of all deviations, in case the input sigmas are mis-scaled (as those
  // of the original PCS d*TREK are).

  std::vector<double> vdSumW(m_vdJ.size(), 0.0), vdSumWI(m_vdJ.size(), 0.0);
  for (size_t j = 0; j < m_vObs.size(); j++)
    {
      const tObs& roObs = m_vObs[j];
      if (!roObs.bUse || roObs.bReject)
        continue;
      double dA = m_vdRunScale[roObs.nRun] * dNormAt(roObs.fLambda);
      double dW = dWeight(roObs) * dA * dA;
      vdSumW[roObs.nUnique]  += dW;
      vdSumWI[roObs.nUnique] += dW * roObs.fI / dA;
    }

  std::vector<double> vdZ(m_vObs.size(), 0.0), vdAbsZ;
  for (size_t j = 0; j < m_vObs.size(); j++)
    {
      const tObs& roObs = m_vObs[j];
      if (!roObs.bUse || roObs.bReject || (3 > m_vnUsedCount[roObs.nUnique]))
        continue;
      double dA     = m_vdRunScale[roObs.nRun] * dNormAt(roObs.fLambda);
      double dW     = dWeight(roObs) * dA * dA;
      double dWRest = vdSumW[roObs.nUnique] - dW;
      if (0.0 >= dWRest)
        continue;
      double dMeanRest = (vdSumWI[roObs.nUnique] - dW * roObs.fI / dA) / dWRest;
      // variance of (I/A - mean of the others) = 1/w + 1/wrest
      vdZ[j] = (roObs.fI / dA - dMeanRest) / sqrt(1.0 / dW + 1.0 / dWRest);
      vdAbsZ.push_back(fabs(vdZ[j]));
    }
  if (vdAbsZ.empty())
    return (0);
  std::nth_element(vdAbsZ.begin(), vdAbsZ.begin() + vdAbsZ.size() / 2, vdAbsZ.end());
  double dSpread = std::max(1.0, 1.4826 * vdAbsZ[vdAbsZ.size() / 2]);
  double dLimit  = m_dReject * dSpread;

  // Worst observation of each reflection
  std::vector<int> vnWorst(m_vdJ.size(), -1);
  for (size_t j = 0; j < m_vObs.size(); j++)
    {
      if (fabs(vdZ[j]) <= dLimit)
        continue;
      int nU = m_vObs[j].nUnique;
      if ( (0 > vnWorst[nU]) || (fabs(vdZ[j]) > fabs(vdZ[vnWorst[nU]])) )
        vnWorst[nU] = (int) j;
    }
  int nNew = 0;
  for (size_t u = 0; u < vnWorst.size(); u++)
    if (0 <= vnWorst[u])
      {
        m_vObs[vnWorst[u]].bReject = true;
        nNew++;
      }
  m_nNumRejected += nNew;
  m_dSigmaScale   = dSpread;
  return (nNew);
}

int
CLaueNorm::nFit(void)
{
  if (m_vObs.empty())
    {
      cout << "ERROR in CLaueNorm::nFit: no observations.\n" << flush;
      return (1);
    }

  // Wavelength range: default from the data

  if (m_a2dLambda[1] <= m_a2dLambda[0])
    {
      m_a2dLambda[0] = 1.0e10;
      m_a2dLambda[1] = 0.0;
      for (size_t j = 0; j < m_vObs.size(); j++)
        {
          m_a2dLambda[0] = std::min(m_a2dLambda[0], (double) m_vObs[j].fLambda);
          m_a2dLambda[1] = std::max(m_a2dLambda[1], (double) m_vObs[j].fLambda);
        }
    }

  int nNumUsed = 0;
  for (size_t j = 0; j < m_vObs.size(); j++)
    {
      tObs& roObs = m_vObs[j];
      roObs.bReject = false;
      roObs.bUse    =    (roObs.fLambda >= m_a2dLambda[0])
                      && (roObs.fLambda <= m_a2dLambda[1])
                      && (0.0f < roObs.fI)
                      && (roObs.fI >= m_dSigmaCut * roObs.fSig);
      if (roObs.bUse)
        nNumUsed++;
    }

  m_vdKnot.assign(m_nKnots, 1.0);
  std::fill(m_vdRunScale.begin(), m_vdRunScale.end(), 1.0);
  m_nNumRejected = 0;
  m_nCyclesDone  = 0;

  cout << "\nTOF-Laue wavelength normalization\n"
       << "  Observations:              " << m_vObs.size() << '\n'
       << "  Unique reflections:        " << m_vdJ.size() << '\n'
       << "  Runs:                      " << m_vdRunScale.size() << '\n'
       << "  Used for the fit:          " << nNumUsed
       << "  (I > " << m_dSigmaCut << " sigma, "
       << m_a2dLambda[0] << " <= lambda <= " << m_a2dLambda[1] << ")\n"
       << "  Knots:                     " << m_nKnots << '\n'
       << "  Smoothness restraint:      " << m_dSmooth << '\n'
       << "  Error model factor:        " << m_dErrFac << '\n'
       << "  Outlier rejection (sigma): " << m_dReject << '\n' << flush;

  int nPass;
  for (nPass = 0; nPass < 10; nPass++)
    {
      int nCycle;
      for (nCycle = 0; nCycle < m_nMaxCycles; nCycle++)
        {
          std::vector<double> vdOldKnot  = m_vdKnot;
          std::vector<double> vdOldScale = m_vdRunScale;
          vUpdateJ();
          if (0 != nUpdateKnots())
            {
              cout << "ERROR in CLaueNorm::nFit: no multiply measured reflections"
                   << " to determine the wavelength curve.\n" << flush;
              return (2);
            }
          vUpdateJ();
          vUpdateRunScales();
          vNormalize();
          m_nCyclesDone++;

          double dChange = 0.0;
          for (size_t m = 0; m < m_vdKnot.size(); m++)
            dChange = std::max(dChange, fabs(m_vdKnot[m] - vdOldKnot[m])
                                        / std::max(1.0e-6, m_vdKnot[m]));
          for (size_t r = 0; r < m_vdRunScale.size(); r++)
            dChange = std::max(dChange, fabs(m_vdRunScale[r] - vdOldScale[r])
                                        / std::max(1.0e-6, m_vdRunScale[r]));
          if (1.0e-5 > dChange)
            break;
        }
      vUpdateJ();
      if ( (0.0 >= m_dReject) || (0 == nRejectOutliers()) )
        break;
    }
  vUpdateJ();
  return (0);
}

void
CLaueNorm::vStats(const int nNumBins, std::vector<double>& rvdNum,
                  std::vector<double>& rvdDen, std::vector<int>& rvnCount,
                  const bool bByRun)
{
  // Laue merging R = sum |I' - <I'>| / sum I' of the normalized intensities
  // of multiply measured reflections, by wavelength bin or by run.

  rvdNum.assign(nNumBins, 0.0);
  rvdDen.assign(nNumBins, 0.0);
  rvnCount.assign(nNumBins, 0);
  for (size_t j = 0; j < m_vObs.size(); j++)
    {
      const tObs& roObs = m_vObs[j];
      if (!roObs.bUse || roObs.bReject || (2 > m_vnUsedCount[roObs.nUnique]))
        continue;
      double dScale = m_vdRunScale[roObs.nRun] * dNormAt(roObs.fLambda);
      double dI     = roObs.fI / dScale;
      int    nBin;
      if (bByRun)
        nBin = roObs.nRun;
      else
        nBin = (int) ((roObs.fLambda - m_a2dLambda[0])
                      / (m_a2dLambda[1] - m_a2dLambda[0]) * nNumBins);
      nBin = std::max(0, std::min(nNumBins - 1, nBin));
      rvdNum[nBin] += fabs(dI - m_vdJ[roObs.nUnique]);
      rvdDen[nBin] += dI;
      rvnCount[nBin]++;
    }
}

void
CLaueNorm::vListResults(void)
{
  int nM = m_nKnots;
  double dLStep = (log(m_a2dLambda[1]) - log(m_a2dLambda[0])) / (double) (nM - 1);

  printf("\nNormalization curve (%d cycles, %d outliers rejected,"
         " sigma scale %.2f)\n", m_nCyclesDone, m_nNumRejected, m_dSigmaScale);
  printf("  N(lambda) has a maximum of 1 and corresponds to LAUENORM's\n"
         "  \"normalised inverse lambda curve\"; the intensities are divided by\n"
         "  k_run * N(lambda).\n\n");
  printf("   Lambda       N\n");
  for (int m = 0; m < nM; m++)
    printf("  %7.3f  %7.4f\n", m_a2dLambda[0] * exp(m * dLStep), m_vdKnot[m]);

  std::vector<double> vdNum, vdDen;
  std::vector<int>    vnCount;
  int nBins = 10;
  double dStep = (m_a2dLambda[1] - m_a2dLambda[0]) / (double) nBins;
  vStats(nBins, vdNum, vdDen, vnCount, false);
  double dNumAll = 0.0, dDenAll = 0.0;
  int    nCountAll = 0;
  printf("\nLaue merging R-factors, R = sum|I - <I>| / sum I of the normalized\n"
         "intensities of multiply measured reflections (I+ and I- merged unless\n"
         "-anom; LAUENORM's RFACT1 is the same quantity)\n\n");
  printf("   Lambda range       Nobs       R\n");
  for (int b = 0; b < nBins; b++)
    {
      printf("  %6.3f - %6.3f  %7d  %7.3f\n", m_a2dLambda[0] + b * dStep,
             m_a2dLambda[0] + (b + 1) * dStep, vnCount[b],
             (0.0 < vdDen[b]) ? vdNum[b] / vdDen[b] : 0.0);
      dNumAll += vdNum[b];
      dDenAll += vdDen[b];
      nCountAll += vnCount[b];
    }
  printf("  %-15s  %7d  %7.3f\n", "All", nCountAll,
         (0.0 < dDenAll) ? dNumAll / dDenAll : 0.0);

  int nRuns = (int) m_vdRunScale.size();
  vStats(nRuns, vdNum, vdDen, vnCount, true);
  printf("\nRun scale factors\n\n");
  printf("    Run    Scale     Nobs       R\n");
  for (int r = 0; r < nRuns; r++)
    printf("  %5d  %7.4f  %7d  %7.3f\n", r + 1, m_vdRunScale[r], vnCount[r],
           (0.0 < vdDen[r]) ? vdNum[r] / vdDen[r] : 0.0);
  fflush(stdout);
}

int
CLaueNorm::nWriteCurve(const Cstring& rsFilename)
{
  FILE *pF = fopen(rsFilename.string(), "w");
  if (NULL == pF)
    return (1);
  fprintf(pF, "# TOF-Laue wavelength normalization curve (dtlauenorm)\n"
              "# lambda  N(lambda) (max 1; LAUENORM's normalised inverse lambda curve)\n");
  int    nSteps = 100;
  for (int i = 0; i <= nSteps; i++)
    {
      double dL = m_a2dLambda[0] + (m_a2dLambda[1] - m_a2dLambda[0]) * i / (double) nSteps;
      fprintf(pF, "%8.4f %10.6f\n", dL, dNormAt(dL));
    }
  fclose(pF);
  return (0);
}

int
CLaueNorm::nMakeOutput(Creflnlist& roOut, const bool bKeepRejected)
{
  int nFIWave  = roOut.nExpandGetField(Creflnlist::ms_sfCalcWavelength);
  int nFIBatch = roOut.nExpandGetField(Creflnlist::ms_ssBatch);
  int nFINorm  = roOut.nExpandGetField(Cstring("fLaueNorm"));
  int nNum     = 0;
  char acBatch[16];

  for (size_t j = 0; j < m_vObs.size(); j++)
    {
      const tObs& roObs = m_vObs[j];
      if ( (roObs.fLambda < m_a2dLambda[0]) || (roObs.fLambda > m_a2dLambda[1]) )
        continue;
      if (roObs.bReject && !bKeepRejected)
        continue;
      double dScale = m_vdRunScale[roObs.nRun] * dNormAt(roObs.fLambda);
      Crefln oRefln(&roOut, roObs.a3nHKL[0], roObs.a3nHKL[1], roObs.a3nHKL[2]);
      oRefln.vSetIntensity((float) (roObs.fI / dScale));
      oRefln.vSetSigmaI((float) (roObs.fSig / dScale));
      oRefln.vSetField(nFIWave, roObs.fLambda);
      oRefln.vSetField(nFINorm, (float) dScale);
      snprintf(acBatch, sizeof(acBatch), "%03d", roObs.nRun + 1);
      oRefln.vSetField(nFIBatch, Cstring(acBatch));
      roOut.nInsert(&oRefln);
      nNum++;
    }
  return (nNum);
}
