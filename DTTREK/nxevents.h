#ifndef DT_NXEVENTS_H
#define DT_NXEVENTS_H
//
// nxevents.h
//    Histogram neutron event data (NeXus NXevent_data, as written at SNS,
//    e.g. MANDI) into TOF-Laue histograms, one (y, x, TOF) array per
//    rectangular detector panel, and write them as an NXlauetof file (the
//    layout used for the ESS NMX reduced data), which d*TREK reads as a
//    TOF-Laue "scan" (see nexuspcs.h).
//
//    Geometry comes from the Mantid instrument definition embedded in the
//    NeXus file (entry/instrument/instrument_xml), optionally replaced by
//    an ISAW/Mantid .DetCal calibration (L1, T0 and panel positions).  The
//    goniometer angles are the time averages of DAS logs.  All positions
//    are written in the NeXus frame (z along the beam, y up), in metres.
//

//+Include files

#include <string>
#include <vector>
#include <stdint.h>

struct tEventPanel
{
  std::string sName;           // "bank12"
  int         nBank;           // 12
  bool        bHasEvents;

  // Pixel ID layout of a Mantid rectangular_detector
  long        lIdStart;
  bool        bFillByY;        // idfillbyfirst="y"
  long        lIdStepByRow;
  long        lIdStep;
  int         nXPixels;        // Unbinned
  int         nYPixels;

  // Geometry (NeXus frame, metres): centre of pixel (0,0), unit vectors
  // along increasing x and y index, pixel pitch
  double      a3dOrigin[3];
  double      a3dXAxis[3];
  double      a3dYAxis[3];
  double      dXSize;
  double      dYSize;
  bool        bCalibrated;     // From a DetCal file

  // Histogram after binning: [y][x][tof]
  int         nXBinned;
  int         nYBinned;
  std::vector<uint32_t> vnCounts;
  uint64_t    llEvents;        // Events histogrammed
};

class CNXEventBinner
{
public:

CNXEventBinner();
~CNXEventBinner();

// Read the instrument, logs and metadata of an event NeXus file.
int  nOpen(const std::string& rsFile);

// Replace L1, T0 and the panel geometry by an ISAW/Mantid .DetCal file.
int  nApplyDetCal(const std::string& rsFile);

// Record a Mantid goniometer calibration (chi-offset, goniometer-tilt).
// The values are written to the output for reference; they are not applied.
int  nReadGonioCal(const std::string& rsFile);

// Goniometer: DAS log names and rotation axes (NeXus frame), outermost first.
// Default: omega (0 1 0), chi (0 0 1), phi (0 1 0) (Mantid "Universal").
void vSetGonio(const std::vector<std::string>& rvsNames,
               const std::vector<double>& rvdAxes);

// TOF binning (microseconds).  dRatio > 0: logarithmic bins with
// t[i+1]/t[i] = 1 + dRatio; dWidth > 0: linear bins.  If no range is set
// it is the range of the event TOFs.
inline void vSetTOFRange(const double dMin, const double dMax)
  { m_a2dTOFRange[0] = dMin; m_a2dTOFRange[1] = dMax; }
inline void vSetLogBins(const double dRatio)   { m_dLogRatio = dRatio; m_dLinWidth = 0.0; }
inline void vSetLinearBins(const double dWidth) { m_dLinWidth = dWidth; m_dLogRatio = 0.0; }
inline void vSetPixelBinning(const int nBin)    { m_nPixBin = (1 > nBin) ? 1 : nBin; }
void vSetBanks(const std::vector<int>& rvnBanks);

int  nBin(void);
int  nWrite(const std::string& rsFile);
void vList(void) const;

inline int nGetNumPanels(void) const { return ((int) m_vPanels.size()); }

private:

std::string               m_sFile;
std::vector<tEventPanel>  m_vPanels;
std::string               m_sInstrument;
double                    m_a3dSource[3];
double                    m_a3dSample[3];
double                    m_dL1;            // m
double                    m_dT0;            // us, added to the event TOF
std::string               m_sDetCal;
std::vector<int>          m_vnBanks;        // Restrict to these (empty: all)

std::vector<std::string>  m_vsGonioNames;
std::vector<double>       m_vdGonioAxes;    // 3 per axis
std::vector<double>       m_vdGonioValues;
std::string               m_sGonioCal;
double                    m_dChiOffset;
std::vector<double>       m_vdGonioTilt;

std::string               m_sTitle, m_sRunNumber, m_sStartTime, m_sEndTime;
double                    m_dDuration;
double                    m_dProtonCharge;

double                    m_a2dTOFRange[2];
double                    m_dLogRatio;
double                    m_dLinWidth;
int                       m_nPixBin;
std::vector<double>       m_vdEdges;        // us

int  nParseIDF(const std::string& rsXML);
int  nReadGonio(long lFile);
int  nTOFRangeFromEvents(long lFile);
int  nMakeEdges(void);
};

#endif   // DT_NXEVENTS_H
