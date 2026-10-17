#ifndef DT_NEXUSPCS_H
#define DT_NEXUSPCS_H
//
// nexuspcs.h
//    Reader for the NeXus time-of-flight histogram files written at the
//    LANSCE Protein Crystallography Station (PCS) by ldaq2nx.
//
//    The original files are HDF4 (*.nx.hdf); they must be converted to HDF5
//    first (h4toh5, see scripts/pcs_h4toh5.sh).  Each file holds one crystal
//    setting as eight detector segments
//       /run_N/tof_detector_K/data   [TOF][arc][vertical], K = 1..8
//    which are stitched into one 256 x 2048 x 96 volume:
//       px0 (fast, vertical) = vertical index
//       px1 (slow, arc)      = K-1 * stride + (arc - cutoff)
//    where stride = segment width - overlap and cutoff is the number of
//    border pixels dropped on each side of a segment (8 binned pixels).
//
//    Each TOF channel is presented to the rest of d*TREK as one image of a
//    scan whose "rotation" axis is the channel index (see CTOFLaue.h).
//
//    The same class reads NXlauetof histogram files (ESS NMX reduced data,
//    or event data histogrammed by dteventbin, e.g. SNS MANDI): one
//    (y, x, TOF) array per flat rectangular detector panel.  The panels are
//    tiled into one mosaic image per channel, panel after panel along px1
//    with a gap of masked pixels between them:
//       px0 = x index,  px1 = panel start + y index
//    and the header describes each panel (Panels_spatial, see Cspatial.h)
//    in the d*TREK laboratory frame (NeXus x -> -X, y -> Y, z -> -Z; mm).
//

//+Include files

#include <vector>
#include "Cstring.h"
#include "Cimage_header.h"

//+Definitions and constants

// Special channel numbers for nReadPCSNexusData() and the requested channel
const int nPCS_CHANNEL_SUM = -1;     // Sum over all TOF channels (a Laue pattern)
const int nPCS_CHANNEL_MAX = -2;     // Maximum over all TOF channels (an "overlay")

// A flat detector panel of a multi-panel (NXlauetof) detector, in the
// d*TREK laboratory frame (mm)
struct tTOFPanel
{
  int    nBank;          // Bank number of the instrument
  int    nStart;         // First mosaic row (px1) of the panel
  int    nNx;            // Pixels along px0 (x index)
  int    nNy;            // Pixels along px1 (y index)
  double a3dOrigin[3];   // Centre of pixel (0, 0)
  double a3dU[3];        // Unit vector of increasing px0
  double a3dV[3];        // Unit vector of increasing px1
  double dDu;            // Pixel pitch along px0
  double dDv;            // Pixel pitch along px1
};

const int nTOF_PANEL_GAP = 8;        // Masked pixels between panels in the mosaic

// Detector response correction modes (see CPCSVolume::nSetResponse)
const int nPCS_RESPONSE_NONE    = 0;
const int nPCS_RESPONSE_PIXEL   = 1;
const int nPCS_RESPONSE_SEGMENT = 2;
const int nPCS_RESPONSE_FULL    = 3;
const int nPCS_RESPONSE_SPECTRUM = 4;

//+Code begin

class DTREK_EXPORT CPCSVolume {

public:

static CPCSVolume* poGetInstance(void);
static bool        bIsHDFSignature(const char* pcBuffer);   // HDF5, or HDF4 (converted copy used)

// Load and stitch a file (no-op if it is already loaded).
int  nLoad(const Cstring& rsFilename);

// Detector response correction from an isotropic incoherent scatterer (V/Nb)
// run in the same format, applied when data are loaded (-cor file -type N):
//   nPCS_RESPONSE_PIXEL:    per-pixel efficiency P(x,y) (summed over TOF)
//   nPCS_RESPONSE_SEGMENT:  P(x,y) times the TOF shape of each detector
//                           segment relative to the whole detector
//   nPCS_RESPONSE_FULL:     per-pixel, per-channel response, smoothed
//   nPCS_RESPONSE_SPECTRUM: the V/Nb TOF spectrum of the whole detector,
//                           normalized to a mean of 1 over the channels.
//                           This REMOVES the incident spectrum, as the
//                           original PCS d*TREK did with -type 3.
// The pixel, segment and full responses have a mean of 1 over the detector at
// each TOF channel, so they do NOT remove the incident spectrum (that is left
// to the wavelength normalization); with them, pixels with almost no V/Nb
// counts are set to -1 (bad).
int  nSetResponse(const Cstring& rsFilename, const int nMode);

// Build a complete d*TREK header for channel nChan (or nPCS_CHANNEL_SUM/MAX).
int  nBuildHeader(Cimage_header* poHeader, const Cstring& rsNameForTemplate,
                  const int nChan);

// Copy one image into pfData (m_nDim0 * m_nDim1 floats, dim0 fastest).
// nChan may be a channel index or nPCS_CHANNEL_SUM / nPCS_CHANNEL_MAX.
int  nGetImage(const int nChan, float* pfData);
int  nGetSum(const int nChan0, const int nChan1, float* pfData);
int  nGetMax(const int nChan0, const int nChan1, float* pfData);

// Direct read-only access to one channel ([dim1][dim0], dim0 fastest)
inline const float* pfGetChannel(const int nChan) const
  { return ( (0 <= nChan) && (nChan < m_nNumChan) )
           ? &m_vfData[(size_t) nChan * m_nDim0 * m_nDim1] : NULL; }

inline int nGetDim0(void)        const { return (m_nDim0); }
inline int nGetDim1(void)        const { return (m_nDim1); }
inline int nGetNumChannels(void) const { return (m_nNumChan); }
inline const Cstring& sGetFilename(void) const { return (m_sFilename); }

private:

CPCSVolume();
~CPCSVolume();

void vClear(void);
int  nApplyResponse(void);
int  nLoadLauetof(long long llFile, const Cstring& rsFilename);
int  nBuildPanelHeader(Cimage_header* poHeader, const Cstring& rsNameForTemplate,
                       const int nChan);

Cstring             m_sFilename;
long                m_lFileTime;      // Modification time of the loaded file
Cstring             m_sResponseFile;  // V/Nb response file ("" for none)
int                 m_nResponseMode;
int                 m_nSegStride;     // Arc pixels between segment starts
int                 m_nNumSeg;

int                 m_nDim0;          // Vertical pixels (fast)
int                 m_nDim1;          // Arc pixels (slow)
int                 m_nNumChan;       // TOF channels
std::vector<float>  m_vfData;         // [chan][dim1][dim0]

// Metadata needed for the header

Cstring             m_sRunName;
std::vector<double> m_vdTOFEdges;     // In units of m_dTimeBin us
double              m_dTimeBin;
double              m_dTimeDelay;
double              m_dL1;
double              m_a2dBeamCenter[2];
double              m_a3dDetSize[3];
double              m_a2dPixelSize[2];
double              m_a3dFastVec[3];
double              m_a3dSlowVec[3];
std::vector<Cstring> m_vsDetGonioNames;
std::vector<double>  m_vdDetGonioValues;
std::vector<double>  m_vdDetGonioVectors;
std::vector<Cstring> m_vsCrysGonioNames;
std::vector<double>  m_vdCrysGonioValues;
std::vector<double>  m_vdCrysGonioVectors;
Cstring              m_sCrysGonioType;
std::vector<double>  m_vdUnitCell;
std::vector<double>  m_vdOrientAngles;
int                  m_nSpacegroup;
double               m_dMosaicity;
Cstring              m_sDetType;
Cstring              m_sDetModel;
Cstring              m_sStartTime;
Cstring              m_sEndTime;
Cstring              m_sComment;
double               m_dMicroAmpHours;

// Multi-panel (NXlauetof) detector
bool                   m_bPanels;
std::vector<tTOFPanel> m_vPanels;
double                 m_dTOFBinRatio;   // Median t[i+1]/t[i] - 1 of the channels
Cstring                m_sInstrument;

static CPCSVolume*   ms_poInstance;
};

// Map a -type value to a response correction mode: "none", "pixel",
// "segment", "full" or "spectrum".  The original PCS d*TREK used -type 3 (dtfind) and
// -type 204 (dtintegrate), whose meaning is not documented; they are
// accepted: 3 removes the incident spectrum ("spectrum", verified against
// legacy dtfind intensities) and 204 is mapped to "full" (legacy dtintegrate
// intensities have no incident-spectrum correction; with pixel, segment and
// full they correlate equally well, CC 0.98, full marginally best).
// Returns -1 if unknown.
int  nPCSResponseMode(const Cstring& rsType);

// Functions called from Cimage_header::nRead and Cimage::nRead

void vSetPCSRequestedChannel(const int nChan);
int  nGetPCSRequestedChannel(void);

int  nReadPCSNexusHeader(const Cstring& rsFilename,
                         const Cstring& rsNameForTemplate,
                         Cimage_header* poHeader);
int  nReadPCSNexusData(const Cstring& rsFilename, float* pfData,
                       Cimage_header* poHeader);

#endif   // DT_NEXUSPCS_H
