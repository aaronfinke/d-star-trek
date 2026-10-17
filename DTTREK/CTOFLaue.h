#ifndef DT_CTOFLAUE_H
#define DT_CTOFLAUE_H
//
// CTOFLaue.h
//    This file is the header file for class CTOFLaue, which encapsulates the
//    time-of-flight (TOF) wavelength-resolved Laue geometry used with neutron
//    data from the LANSCE Protein Crystallography Station (PCS).
//
//    In TOF-Laue mode the third image dimension is a TOF channel, not a
//    rotation angle.  As in the original PCS version of d*TREK, the channel
//    index is carried in the "rotation" fields of scans and reflection lists
//    (SCAN_ROTATION=0 96 1, fObs_rot_mid, fCalc_rot_mid, ...), while the
//    crystal orientation is fixed by CRYSTAL_GONIO_VALUES.  This class maps
//    channels to TOF and wavelength:
//
//       t(c)   = linear interpolation of the TOF bin edges at channel c  (us)
//       lambda = (t(c) + delay) / (K * (L1 + L2))
//
//    with K = m_n/h = 252.7784 us/(m*Angstrom), L1 the moderator-to-sample
//    distance and L2 the sample-to-pixel distance.
//
//    Header keywords (the ROTATION_* keywords are the ones written by the
//    original PCS d*TREK and carry the detector prefix, e.g. Detector1_):
//       <prefix>ROTATION_ARRAY  TOF bin edges (N+1 values, units of ROTATION_BIN us)
//       <prefix>ROTATION_BIN    multiplier from ROTATION_ARRAY units to us
//       <prefix>ROTATION_DELAY  TOF offset added to t(c) (us)
//       <prefix>ROTATION_L1     moderator-to-sample distance (mm)
//       TOF_WAVELENGTH_RANGE    lambda_min lambda_max (Angstrom)
//       TOF_WIDTH_MODEL         a b c: predicted spot width in channels
//                               = a + b*lambda + c/sin(theta).  The default
//                               2 0 0.5 is what the original PCS d*TREK wrote
//                               as fCalc_rot_width (verified on 2010-2013 runs).
//

//+Include files

#include <vector>
#include "Cstring.h"
#include "Cimage_header.h"

//+Code begin

class DTREK_EXPORT CTOFLaue {

public:

static const double ms_dTOFConst;     // m_n/h in us/(m*Angstrom)

static Cstring ms_sRotationArray;     // "ROTATION_ARRAY"  (detector-prefixed)
static Cstring ms_sRotationBin;       // "ROTATION_BIN"    (detector-prefixed)
static Cstring ms_sRotationDelay;     // "ROTATION_DELAY"  (detector-prefixed)
static Cstring ms_sRotationL1;        // "ROTATION_L1"     (detector-prefixed)
static Cstring ms_sWavelengthRange;   // "TOF_WAVELENGTH_RANGE"
static Cstring ms_sWidthModel;        // "TOF_WIDTH_MODEL"

private:

bool                m_bAvailable;
Cstring             m_sPrefix;        // Detector prefix, e.g. "Detector1_"
std::vector<double> m_vdEdges;        // TOF bin edges in us (N+1 values)
double              m_dDelay;         // us
double              m_dL1;            // mm
double              m_a2dLambdaRange[2];  // Angstrom, [0] min, [1] max
double              m_a3dWidthModel[3];   // channels: [0] + [1]*lambda + [2]/sin(theta)

public:

CTOFLaue();
CTOFLaue(Cimage_header& roHeader, const Cstring& rsPrefix = "");
~CTOFLaue();

// True if the header describes TOF-Laue data (has ROTATION_ARRAY and ROTATION_L1)
static bool    bIsTOFLaue(Cimage_header& roHeader);
static Cstring sGetDetectorPrefix(Cimage_header& roHeader);

// Add keywords that the rest of d*TREK needs but that headers written by
// the original PCS d*TREK lack (e.g. <prefix>DETECTOR_DIMENSIONS).
static int     nNormalizeHeader(Cimage_header& roHeader);

int  nInitValues(Cimage_header& roHeader, const Cstring& rsPrefix = "");
int  nUpdateHeader(Cimage_header* poHeader);

inline bool   bIsAvailable(void) const      { return (m_bAvailable); }
inline int    nGetNumChannels(void) const   { return ((int)m_vdEdges.size() - 1); }
inline double dGetL1(void) const            { return (m_dL1); }
inline double dGetDelay(void) const         { return (m_dDelay); }
inline void   vSetL1(const double dL1)      { m_dL1 = dL1; }
inline void   vSetDelay(const double dDel)  { m_dDelay = dDel; }
inline double dGetLambdaMin(void) const     { return (m_a2dLambdaRange[0]); }
inline double dGetLambdaMax(void) const     { return (m_a2dLambdaRange[1]); }
void          vSetLambdaRange(const double dLambda1, const double dLambda2);
inline bool   bLambdaInRange(const double dLambda) const
  { return ( (dLambda >= m_a2dLambdaRange[0]) && (dLambda <= m_a2dLambdaRange[1]) ); }
inline void   vSetWidthModel(const double dA, const double dB, const double dC)
  { m_a3dWidthModel[0] = dA; m_a3dWidthModel[1] = dB; m_a3dWidthModel[2] = dC; }
inline double dGetWidthModel(const int nIdx) const { return (m_a3dWidthModel[nIdx]); }

// Channel <-> TOF (us, before the delay is added)
double dChanToTOF(const double dChan) const;
double dTOFToChan(const double dTOF) const;
double dGetChannelWidthTOF(const double dChan) const;   // us per channel at dChan

// Channel <-> wavelength for a sample-to-pixel distance dL2 (mm)
double dChanToLambda(const double dChan, const double dL2) const;
double dLambdaToChan(const double dLambda, const double dL2) const;
double dDChanDLambda(const double dLambda, const double dL2) const;

// Predicted spot width along TOF in channels
double dWidthChan(const double dLambda, const double dSinTheta) const;

// Lorentz factor as written by the original PCS d*TREK: 1/sin^2(theta)
static double dLorentzLegacy(const double dSinTheta);
// Wavelength-resolved Laue Lorentz factor: lambda^4 / (2 sin^2(theta))
static double dLorentzLaue(const double dLambda, const double dSinTheta);

void vList(void);
};

#endif   // DT_CTOFLAUE_H
