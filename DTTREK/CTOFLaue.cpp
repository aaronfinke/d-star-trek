//
// CTOFLaue.cpp
//    This file contains the member functions of class CTOFLaue which maps
//    time-of-flight channels to TOF and wavelength for wavelength-resolved
//    neutron Laue data (LANSCE PCS).  See CTOFLaue.h for the conventions.
//

//+Include files

#include <algorithm>
#include <math.h>
#include <stdio.h>
#include "CTOFLaue.h"
#include "dtrekdefs.h"

#if !defined(VC6)
using std::cout;
using std::flush;
#endif

//+Definitions, constants, and initialization of static member variables

// m_n / h = 1.67492749804e-27 kg / 6.62607015e-34 J s  in us/(m*Angstrom)
const double CTOFLaue::ms_dTOFConst = 252.7784;

Cstring CTOFLaue::ms_sRotationArray   = "ROTATION_ARRAY";
Cstring CTOFLaue::ms_sRotationBin     = "ROTATION_BIN";
Cstring CTOFLaue::ms_sRotationDelay   = "ROTATION_DELAY";
Cstring CTOFLaue::ms_sRotationL1      = "ROTATION_L1";
Cstring CTOFLaue::ms_sWavelengthRange = "TOF_WAVELENGTH_RANGE";
Cstring CTOFLaue::ms_sWidthModel      = "TOF_WIDTH_MODEL";

//+Code begin

CTOFLaue::CTOFLaue()
{
  m_bAvailable        = false;
  m_sPrefix           = "";
  m_dDelay            = 0.0;
  m_dL1               = 0.0;
  m_a2dLambdaRange[0] = 0.0;
  m_a2dLambdaRange[1] = 1000.0;
  m_a3dWidthModel[0]  = 2.0;
  m_a3dWidthModel[1]  = 0.0;
  m_a3dWidthModel[2]  = 0.5;
}

CTOFLaue::CTOFLaue(Cimage_header& roHeader, const Cstring& rsPrefix)
{
  m_bAvailable        = false;
  m_sPrefix           = "";
  m_dDelay            = 0.0;
  m_dL1               = 0.0;
  m_a2dLambdaRange[0] = 0.0;
  m_a2dLambdaRange[1] = 1000.0;
  m_a3dWidthModel[0]  = 2.0;
  m_a3dWidthModel[1]  = 0.0;
  m_a3dWidthModel[2]  = 0.5;
  (void) nInitValues(roHeader, rsPrefix);
}

CTOFLaue::~CTOFLaue()
{
  m_vdEdges.clear();
}

Cstring
CTOFLaue::sGetDetectorPrefix(Cimage_header& roHeader)
{
  Cstring sPrefix = "";
  if (0 != roHeader.nGetValue(Cstring(D_K_DetectorNames), 1, &sPrefix))
    sPrefix = "";
  return (sPrefix);
}

bool
CTOFLaue::bIsTOFLaue(Cimage_header& roHeader)
{
  Cstring sPrefix;
  Cstring sTemp;

  if (!roHeader.bIsAvailable())
    return (false);
  sPrefix = sGetDetectorPrefix(roHeader);
  return (   (0 == roHeader.nGetValue(sPrefix + ms_sRotationArray, &sTemp))
          && (0 == roHeader.nGetValue(sPrefix + ms_sRotationL1, &sTemp)) );
}

int
CTOFLaue::nNormalizeHeader(Cimage_header& roHeader)
{
  // Headers written by the original PCS version of d*TREK store
  //   DETECTOR_DIMENSIONS=256 2048 96;
  // without the detector prefix, while Cdetector requires the prefixed
  // two-value keyword.  Add it if missing.

  int     a3nDim[3];
  Cstring sPrefix;

  if (!bIsTOFLaue(roHeader))
    return (0);

  sPrefix = sGetDetectorPrefix(roHeader);
  if (   ("" != sPrefix)
      && (0 != roHeader.nGetValue(sPrefix + D_K_DetectorDimensions, 2, a3nDim))
      && (0 == roHeader.nGetValue(Cstring(D_K_DetectorDimensions), 2, a3nDim)) )
    {
      (void) roHeader.nReplaceValue(sPrefix + D_K_DetectorDimensions, 2, a3nDim);
    }

  // Neutron counts of 0 are valid (d*TREK otherwise treats pixels below 1
  // as bad); negative values mark bad pixels.

  float fTemp;
  if (0 != roHeader.nGetValue(Cstring("DTREK_NONUNF_OKVALUE"), &fTemp))
    (void) roHeader.nReplaceValue(Cstring("DTREK_NONUNF_OKVALUE"), 0.0f, 1);
  return (0);
}

int
CTOFLaue::nInitValues(Cimage_header& roHeader, const Cstring& rsPrefix)
{
  int     nStat;
  double  dBin;
  double  a2dTemp[2];
  Cstring sTemp;

  m_bAvailable = false;
  m_vdEdges.clear();

  m_sPrefix = rsPrefix;
  if ("" == m_sPrefix)
    m_sPrefix = sGetDetectorPrefix(roHeader);

  nStat = roHeader.nGetValue(m_sPrefix + ms_sRotationArray, &sTemp);
  if (0 != nStat)
    {
      cout << "ERROR in CTOFLaue: missing " << m_sPrefix << ms_sRotationArray
           << " keyword.\n" << flush;
      return (1);
    }

  // Parse the bin edges from the string so there is no limit on their number

  {
    const char *pcStart = sTemp.string();
    char       *pcEnd   = NULL;
    double      dValue;
    while (NULL != pcStart)
      {
        dValue = strtod(pcStart, &pcEnd);
        if (pcEnd == pcStart)
          break;
        m_vdEdges.push_back(dValue);
        pcStart = pcEnd;
      }
  }
  if (2 > m_vdEdges.size())
    {
      cout << "ERROR in CTOFLaue: " << m_sPrefix << ms_sRotationArray
           << " needs at least 2 values.\n" << flush;
      m_vdEdges.clear();
      return (2);
    }

  dBin = 1.0;
  (void) roHeader.nGetValue(m_sPrefix + ms_sRotationBin, &dBin);
  for (size_t i = 0; i < m_vdEdges.size(); i++)
    {
      m_vdEdges[i] *= dBin;
      if ( (0 < i) && (m_vdEdges[i] <= m_vdEdges[i-1]) )
        {
          cout << "ERROR in CTOFLaue: TOF bin edges are not increasing.\n" << flush;
          m_vdEdges.clear();
          return (3);
        }
    }

  m_dDelay = 0.0;
  (void) roHeader.nGetValue(m_sPrefix + ms_sRotationDelay, &m_dDelay);

  nStat = roHeader.nGetValue(m_sPrefix + ms_sRotationL1, &m_dL1);
  if ( (0 != nStat) || (0.0 >= m_dL1) )
    {
      cout << "ERROR in CTOFLaue: missing or invalid " << m_sPrefix
           << ms_sRotationL1 << " keyword.\n" << flush;
      m_vdEdges.clear();
      return (4);
    }

  if (0 == roHeader.nGetValue(ms_sWavelengthRange, 2, a2dTemp))
    vSetLambdaRange(a2dTemp[0], a2dTemp[1]);

  double a3dWidth[3];
  if (0 == roHeader.nGetValue(ms_sWidthModel, 3, a3dWidth))
    vSetWidthModel(a3dWidth[0], a3dWidth[1], a3dWidth[2]);
  else if (0 == roHeader.nGetValue(ms_sWidthModel, 2, a3dWidth))
    vSetWidthModel(a3dWidth[0], a3dWidth[1], 0.0);

  m_bAvailable = true;
  return (0);
}

int
CTOFLaue::nUpdateHeader(Cimage_header* poHeader)
{
  // Only the refinable quantities are written back; the bin edges
  // are a property of the data file and are never changed.

  if (!m_bAvailable)
    return (1);
  (void) poHeader->nReplaceValue(m_sPrefix + ms_sRotationL1, m_dL1, 3);
  (void) poHeader->nReplaceValue(m_sPrefix + ms_sRotationDelay, m_dDelay, 3);
  (void) poHeader->nReplaceValue(ms_sWavelengthRange, 2, m_a2dLambdaRange, 4);
  (void) poHeader->nReplaceValue(ms_sWidthModel, 3, m_a3dWidthModel, 4);
  return (0);
}

void
CTOFLaue::vSetLambdaRange(const double dLambda1, const double dLambda2)
{
  m_a2dLambdaRange[0] = std::min(dLambda1, dLambda2);
  m_a2dLambdaRange[1] = std::max(dLambda1, dLambda2);
}

double
CTOFLaue::dChanToTOF(const double dChan) const
{
  // Channel coordinate c runs from 0 at the start of the first bin to N
  // at the end of the last bin.  Within a bin the TOF is linear in c.
  // Outside [0, N] extrapolate with the width of the end bins.

  int nNum = nGetNumChannels();
  int nBin;

  if (0 >= nNum)
    return (0.0);
  nBin = (int) floor(dChan);
  if (0 > nBin)
    nBin = 0;
  else if (nNum <= nBin)
    nBin = nNum - 1;
  return (m_vdEdges[nBin] + (dChan - (double) nBin)
                            * (m_vdEdges[nBin+1] - m_vdEdges[nBin]));
}

double
CTOFLaue::dTOFToChan(const double dTOF) const
{
  int nNum = nGetNumChannels();
  int nBin;

  if (0 >= nNum)
    return (0.0);

  // First edge greater than dTOF, then step back to the bin containing it

  nBin = (int) (std::upper_bound(m_vdEdges.begin(), m_vdEdges.end(), dTOF)
                - m_vdEdges.begin()) - 1;
  if (0 > nBin)
    nBin = 0;
  else if (nNum <= nBin)
    nBin = nNum - 1;
  return ((double) nBin + (dTOF - m_vdEdges[nBin])
                          / (m_vdEdges[nBin+1] - m_vdEdges[nBin]));
}

double
CTOFLaue::dGetChannelWidthTOF(const double dChan) const
{
  int nNum = nGetNumChannels();
  int nBin;

  if (0 >= nNum)
    return (0.0);
  nBin = (int) floor(dChan);
  if (0 > nBin)
    nBin = 0;
  else if (nNum <= nBin)
    nBin = nNum - 1;
  return (m_vdEdges[nBin+1] - m_vdEdges[nBin]);
}

double
CTOFLaue::dChanToLambda(const double dChan, const double dL2) const
{
  return ( (dChanToTOF(dChan) + m_dDelay)
           / (ms_dTOFConst * (m_dL1 + dL2) * 0.001) );
}

double
CTOFLaue::dLambdaToChan(const double dLambda, const double dL2) const
{
  return (dTOFToChan(dLambda * ms_dTOFConst * (m_dL1 + dL2) * 0.001 - m_dDelay));
}

double
CTOFLaue::dDChanDLambda(const double dLambda, const double dL2) const
{
  double dChan = dLambdaToChan(dLambda, dL2);
  return (ms_dTOFConst * (m_dL1 + dL2) * 0.001 / dGetChannelWidthTOF(dChan));
}

double
CTOFLaue::dWidthChan(const double dLambda, const double dSinTheta) const
{
  double dWidth = m_a3dWidthModel[0] + m_a3dWidthModel[1] * dLambda;
  if (0.0 < dSinTheta)
    dWidth += m_a3dWidthModel[2] / dSinTheta;
  return (std::max(dWidth, 1.0));
}

double
CTOFLaue::dLorentzLegacy(const double dSinTheta)
{
  if (0.0 >= dSinTheta)
    return (0.0);
  return (1.0 / (dSinTheta * dSinTheta));
}

double
CTOFLaue::dLorentzLaue(const double dLambda, const double dSinTheta)
{
  if (0.0 >= dSinTheta)
    return (0.0);
  return (dLambda * dLambda * dLambda * dLambda / (2.0 * dSinTheta * dSinTheta));
}

void
CTOFLaue::vList(void)
{
  int nNum = nGetNumChannels();

  cout << "\nTOF-Laue geometry:\n"
       << "  Number of TOF channels:   " << nNum << '\n';
  if (0 < nNum)
    cout << "  TOF range (us):           " << m_vdEdges[0] << " to "
         << m_vdEdges[nNum] << '\n';
  cout << "  TOF delay (us):           " << m_dDelay << '\n'
       << "  L1 (mm):                  " << m_dL1 << '\n'
       << "  Wavelength range (A):     " << m_a2dLambdaRange[0] << " to "
       << m_a2dLambdaRange[1] << '\n'
       << "  Width model (channels):   " << m_a3dWidthModel[0] << " + "
       << m_a3dWidthModel[1] << " * lambda + " << m_a3dWidthModel[2]
       << " / sin(theta)\n" << flush;
}
