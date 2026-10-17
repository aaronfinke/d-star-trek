#ifndef DT_CLAUENORM_H
#define DT_CLAUENORM_H
//
// CLaueNorm.h
//    Wavelength normalization and run scaling of TOF-Laue intensities, the
//    step done by LAUENORM in the original LANSCE PCS workflow (Langan &
//    Greene, J. Appl. Cryst. 2004, 37, 253-257).
//
//    The integrated intensities of dtintegrate are already divided by the
//    geometric part 1/sin^2(theta) of the Laue Lorentz factor, so an
//    observation j of the unique reflection u in run r is modelled as
//
//       I_j = k_r * N(lambda_j) * J_u
//
//    where N(lambda) (the incident spectrum times the detector efficiency
//    times lambda^4) is piecewise linear between knots equally spaced in
//    ln(lambda), with a second-difference smoothness restraint, k_r is a
//    scale factor per run and J_u is the normalized intensity.  The fit alternates between
//    J (weighted means of symmetry equivalents), N and k, and rejects
//    outliers among multiply measured reflections.
//

//+Include files

#include <vector>
#include <map>
#include "Dtrek.h"
#include "Cstring.h"

class Creflnlist;
class Cspacegroup;

class DTREK_EXPORT CLaueNorm
{
public:

CLaueNorm();
~CLaueNorm();

// Add the observations of one reflection list as run nRun (0-based).
// The wavelength is taken from rsWaveField (fCalwavelength by default,
// else fDetwavelength).  Returns the number of observations added or < 0.
int    nAddList(Creflnlist& roList, Cspacegroup& roSpace, const int nRun,
                const Cstring& rsWaveField, const bool bAnom = false);

inline void vSetLambdaRange(const double dMin, const double dMax)
  { m_a2dLambda[0] = dMin; m_a2dLambda[1] = dMax; }
inline void vSetNumKnots(const int nKnots)     { m_nKnots     = (nKnots < 2) ? 2 : nKnots; }
inline void vSetSmooth(const double dSmooth)   { m_dSmooth    = dSmooth; }
inline void vSetSigmaCut(const double dSig)    { m_dSigmaCut  = dSig; }
inline void vSetRejectLimit(const double dRej) { m_dReject    = dRej; }
inline void vSetErrorFactor(const double dE)   { m_dErrFac    = dE; }
inline void vSetMaxCycles(const int nCycles)   { m_nMaxCycles = nCycles; }

int    nFit(void);

double dGetNorm(const double dLambda) const;     // N(lambda), max of 1
inline double dGetRunScale(const int nRun) const
  { return ( (0 <= nRun) && (nRun < (int) m_vdRunScale.size()) ) ? m_vdRunScale[nRun] : 1.0; }
inline int    nGetNumRuns(void) const   { return ((int) m_vdRunScale.size()); }

// Normalized, unmerged output list (I / (k_r N(lambda)), sBatch = run
// number + 1); outliers are left out unless bKeepRejected.
int    nMakeOutput(Creflnlist& roOut, const bool bKeepRejected = false);
int    nWriteCurve(const Cstring& rsFilename);
void   vListResults(void);

private:

struct tObs
{
  int    nRun;
  int    nUnique;
  int    a3nHKL[3];
  float  fLambda;
  float  fI;
  float  fSig;
  bool   bUse;       // I/sigma and wavelength cuts passed
  bool   bReject;    // outlier
};

std::vector<tObs>    m_vObs;
std::map<int, int>   m_mPackedToUnique;
std::vector<double>  m_vdJ;           // normalized intensity of each unique
std::vector<int>     m_vnUsedCount;   // used, non-rejected observations of each unique
std::vector<double>  m_vdKnot;        // N at the knots
std::vector<double>  m_vdRunScale;
double               m_a2dLambda[2];
int                  m_nKnots;
double               m_dSmooth;
double               m_dSigmaCut;
double               m_dReject;
double               m_dErrFac;
int                  m_nMaxCycles;
int                  m_nCyclesDone;
int                  m_nNumRejected;
double               m_dSigmaScale;   // robust spread of the normalized deviations

void   vBasis(const double dLambda, int *pnKnot, double *pdFrac) const;
double dNormAt(const double dLambda) const;
double dWeight(const tObs& roObs) const;
void   vUpdateJ(void);
int    nUpdateKnots(void);
void   vUpdateRunScales(void);
void   vNormalize(void);
int    nRejectOutliers(void);
void   vStats(const int nNumBins, std::vector<double>& rvdNum,
              std::vector<double>& rvdDen, std::vector<int>& rvnCount, const bool bByRun);
};

#endif   // DT_CLAUENORM_H
