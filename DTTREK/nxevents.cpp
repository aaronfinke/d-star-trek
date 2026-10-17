//
// nxevents.cpp
//    Histogram NeXus event data into TOF-Laue panel histograms and write
//    them as an NXlauetof file.  See nxevents.h.
//

//+Include files

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include "hdf5.h"
#include "h5util.h"
#include "nxevents.h"

using std::cout;
using std::flush;
using std::string;
using std::vector;

//+Local helpers

namespace {

const double c_dDegToRad = 3.14159265358979323846 / 180.0;

// --- 3-vectors and rotation matrices

void
vRotMat(const double dAngleDeg, const double a3dAxisIn[3], double a3x3dR[3][3])
{
  // Right-handed rotation by dAngleDeg about the axis (Rodrigues).
  double a3dAxis[3] = {a3dAxisIn[0], a3dAxisIn[1], a3dAxisIn[2]};
  double dLen = sqrt(a3dAxis[0] * a3dAxis[0] + a3dAxis[1] * a3dAxis[1] + a3dAxis[2] * a3dAxis[2]);
  if (0.0 >= dLen)
    {
      a3dAxis[0] = 0.0; a3dAxis[1] = 0.0; a3dAxis[2] = 1.0;
      dLen = 1.0;
    }
  double x = a3dAxis[0] / dLen, y = a3dAxis[1] / dLen, z = a3dAxis[2] / dLen;
  double c = cos(dAngleDeg * c_dDegToRad), s = sin(dAngleDeg * c_dDegToRad), t = 1.0 - c;
  a3x3dR[0][0] = t*x*x + c;   a3x3dR[0][1] = t*x*y - s*z; a3x3dR[0][2] = t*x*z + s*y;
  a3x3dR[1][0] = t*x*y + s*z; a3x3dR[1][1] = t*y*y + c;   a3x3dR[1][2] = t*y*z - s*x;
  a3x3dR[2][0] = t*x*z - s*y; a3x3dR[2][1] = t*y*z + s*x; a3x3dR[2][2] = t*z*z + c;
}

void
vMatMul(const double a[3][3], const double b[3][3], double c[3][3])
{
  double t[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      t[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
  memcpy(c, t, sizeof(t));
}

void
vMatVec(const double a[3][3], const double v[3], double w[3])
{
  double t[3];
  for (int i = 0; i < 3; i++)
    t[i] = a[i][0] * v[0] + a[i][1] * v[1] + a[i][2] * v[2];
  w[0] = t[0]; w[1] = t[1]; w[2] = t[2];
}

double
dNorm(const double v[3])
{
  return (sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]));
}

// --- Minimal XML tag scanner for Mantid instrument definitions

struct tTag
{
  string                    sName;     // tag name without '/'
  bool                      bClose;    // </name>
  bool                      bSelfClose;// <name ... />
  std::map<string, string>  mAttr;
};

bool
bNextTag(const string& rsXML, size_t& rnPos, tTag& rTag)
{
  while (true)
    {
      size_t nStart = rsXML.find('<', rnPos);
      if (string::npos == nStart)
        return (false);
      if (0 == rsXML.compare(nStart, 4, "<!--"))
        {
          size_t nEnd = rsXML.find("-->", nStart + 4);
          rnPos = (string::npos == nEnd) ? rsXML.size() : nEnd + 3;
          continue;
        }
      if ( (0 == rsXML.compare(nStart, 2, "<?")) || (0 == rsXML.compare(nStart, 2, "<!")) )
        {
          size_t nEnd = rsXML.find('>', nStart);
          rnPos = (string::npos == nEnd) ? rsXML.size() : nEnd + 1;
          continue;
        }
      size_t nEnd = rsXML.find('>', nStart);
      if (string::npos == nEnd)
        return (false);
      rnPos = nEnd + 1;

      string sBody = rsXML.substr(nStart + 1, nEnd - nStart - 1);
      rTag.mAttr.clear();
      rTag.bClose     = (!sBody.empty() && ('/' == sBody[0]));
      rTag.bSelfClose = (!sBody.empty() && ('/' == sBody[sBody.size() - 1]));
      if (rTag.bClose)
        sBody = sBody.substr(1);
      if (rTag.bSelfClose)
        sBody = sBody.substr(0, sBody.size() - 1);
      size_t i = 0;
      while ( (i < sBody.size()) && !isspace((unsigned char) sBody[i]) )
        i++;
      rTag.sName = sBody.substr(0, i);
      // key="value" pairs
      while (i < sBody.size())
        {
          while ( (i < sBody.size()) && isspace((unsigned char) sBody[i]) )
            i++;
          size_t nEq = sBody.find('=', i);
          if (string::npos == nEq)
            break;
          string sKey = sBody.substr(i, nEq - i);
          sKey.erase(sKey.find_last_not_of(" \t\r\n") + 1);
          size_t nQ1 = sBody.find_first_of("\"'", nEq + 1);
          if (string::npos == nQ1)
            break;
          size_t nQ2 = sBody.find(sBody[nQ1], nQ1 + 1);
          if (string::npos == nQ2)
            break;
          rTag.mAttr[sKey] = sBody.substr(nQ1 + 1, nQ2 - nQ1 - 1);
          i = nQ2 + 1;
        }
      return (true);
    }
}

double
dAttr(const tTag& rTag, const string& rsKey, const double dDefault)
{
  std::map<string, string>::const_iterator it = rTag.mAttr.find(rsKey);
  return ( (rTag.mAttr.end() == it) ? dDefault : atof(it->second.c_str()) );
}

string
sAttr(const tTag& rTag, const string& rsKey)
{
  std::map<string, string>::const_iterator it = rTag.mAttr.find(rsKey);
  return ( (rTag.mAttr.end() == it) ? string("") : it->second );
}

// Position of a <location> tag (x y z, or r t p spherical), metres
void
vLocation(const tTag& rTag, double a3dPos[3])
{
  if (rTag.mAttr.count("r"))
    {
      double r = dAttr(rTag, "r", 0.0);
      double t = dAttr(rTag, "t", 0.0) * c_dDegToRad;
      double p = dAttr(rTag, "p", 0.0) * c_dDegToRad;
      a3dPos[0] = r * sin(t) * cos(p);
      a3dPos[1] = r * sin(t) * sin(p);
      a3dPos[2] = r * cos(t);
    }
  else
    {
      a3dPos[0] = dAttr(rTag, "x", 0.0);
      a3dPos[1] = dAttr(rTag, "y", 0.0);
      a3dPos[2] = dAttr(rTag, "z", 0.0);
    }
}

int
nBankNumber(const string& rsName)
{
  size_t i = rsName.find_first_of("0123456789");
  return ( (string::npos == i) ? -1 : atoi(rsName.c_str() + i) );
}

}  // end of anonymous namespace

//+Code begin

CNXEventBinner::CNXEventBinner()
{
  for (int i = 0; i < 3; i++)
    {
      m_a3dSource[i] = 0.0;
      m_a3dSample[i] = 0.0;
    }
  m_dL1           = 0.0;
  m_dT0           = 0.0;
  m_dChiOffset    = 0.0;
  m_dDuration     = 0.0;
  m_dProtonCharge = 0.0;
  m_a2dTOFRange[0] = m_a2dTOFRange[1] = 0.0;
  m_dLogRatio     = 0.005;
  m_dLinWidth     = 0.0;
  m_nPixBin       = 1;
  const char *apcNames[3] = {"omega", "chi", "phi"};
  const double a9dAxes[9] = {0, 1, 0,  0, 0, 1,  0, 1, 0};
  for (int i = 0; i < 3; i++)
    m_vsGonioNames.push_back(apcNames[i]);
  m_vdGonioAxes.assign(a9dAxes, a9dAxes + 9);
}

CNXEventBinner::~CNXEventBinner()
{
}

void
CNXEventBinner::vSetGonio(const vector<string>& rvsNames, const vector<double>& rvdAxes)
{
  m_vsGonioNames = rvsNames;
  if (rvdAxes.size() == 3 * rvsNames.size())
    m_vdGonioAxes = rvdAxes;
}

void
CNXEventBinner::vSetBanks(const vector<int>& rvnBanks)
{
  m_vnBanks = rvnBanks;
}

int
CNXEventBinner::nParseIDF(const string& rsXML)
{
  // Rectangular detector types, the components using them (pixel ID layout
  // and location with nested rotations), the source and sample positions.
  // Components are taken at the top level (no assemblies), which covers
  // the SNS single-crystal instruments.

  struct tRectType { double xstart, xstep, ystart, ystep; int nx, ny; };
  std::map<string, tRectType> mRect;
  std::map<string, string>    mIs;      // type name -> "Source", "SamplePos"
  size_t nPos = 0;
  tTag   oTag;

  // Pass 1: types
  while (bNextTag(rsXML, nPos, oTag))
    {
      if ( ("type" != oTag.sName) || oTag.bClose )
        continue;
      string sName = sAttr(oTag, "name");
      string sIs   = sAttr(oTag, "is");
      std::transform(sIs.begin(), sIs.end(), sIs.begin(), ::tolower);
      if ("rectangular_detector" == sIs)
        {
          tRectType oR;
          oR.xstart = dAttr(oTag, "xstart", 0.0);
          oR.xstep  = dAttr(oTag, "xstep", 0.0);
          oR.ystart = dAttr(oTag, "ystart", 0.0);
          oR.ystep  = dAttr(oTag, "ystep", 0.0);
          oR.nx     = (int) dAttr(oTag, "xpixels", 0.0);
          oR.ny     = (int) dAttr(oTag, "ypixels", 0.0);
          mRect[sName] = oR;
        }
      else if ( ("source" == sIs) || ("samplepos" == sIs) )
        mIs[sName] = sIs;
    }

  // Pass 2: components and their locations
  nPos = 0;
  string sCompType;
  tTag   oComp;
  bool   bInComp = false;
  while (bNextTag(rsXML, nPos, oTag))
    {
      if ("component" == oTag.sName)
        {
          if (oTag.bClose)
            {
              bInComp = false;
              continue;
            }
          oComp     = oTag;
          sCompType = sAttr(oTag, "type");
          bInComp   = !oTag.bSelfClose;
          continue;
        }
      if ( !bInComp || ("location" != oTag.sName) || oTag.bClose )
        continue;

      // Location, then nested <rot> tags (outermost first) up to </location>
      double a3dPos[3];
      double a3x3dR[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
      vLocation(oTag, a3dPos);
      tTag oLoc = oTag;
      if (oLoc.mAttr.count("rot"))
        {
          double a3dAxis[3] = {dAttr(oLoc, "axis-x", 0.0), dAttr(oLoc, "axis-y", 0.0),
                               dAttr(oLoc, "axis-z", 1.0)};
          vRotMat(dAttr(oLoc, "rot", 0.0), a3dAxis, a3x3dR);
        }
      if (!oLoc.bSelfClose)
        {
          tTag oInner;
          while (bNextTag(rsXML, nPos, oInner))
            {
              if ( ("location" == oInner.sName) && oInner.bClose )
                break;
              if ( ("rot" == oInner.sName) && !oInner.bClose )
                {
                  double a3dAxis[3] = {dAttr(oInner, "axis-x", 0.0), dAttr(oInner, "axis-y", 0.0),
                                       dAttr(oInner, "axis-z", 1.0)};
                  double a3x3dRot[3][3];
                  vRotMat(dAttr(oInner, "val", 0.0), a3dAxis, a3x3dRot);
                  vMatMul(a3x3dR, a3x3dRot, a3x3dR);
                }
            }
        }

      std::map<string, string>::const_iterator itIs = mIs.find(sCompType);
      if (mIs.end() != itIs)
        {
          double *pdDest = ("source" == itIs->second) ? m_a3dSource : m_a3dSample;
          for (int i = 0; i < 3; i++)
            pdDest[i] = a3dPos[i];
          continue;
        }
      std::map<string, tRectType>::const_iterator itR = mRect.find(sCompType);
      if (mRect.end() == itR)
        continue;

      const tRectType& rR = itR->second;
      tEventPanel oP;
      oP.sName        = sAttr(oLoc, "name");
      if ("" == oP.sName)
        oP.sName = sAttr(oComp, "name");
      if ("" == oP.sName)
        oP.sName = sCompType;
      oP.nBank        = nBankNumber(oP.sName);
      oP.bHasEvents   = false;
      oP.lIdStart     = (long) dAttr(oComp, "idstart", 0.0);
      string sFill    = sAttr(oComp, "idfillbyfirst");
      oP.bFillByY     = ("x" != sFill);
      oP.lIdStepByRow = (long) dAttr(oComp, "idstepbyrow", oP.bFillByY ? rR.ny : rR.nx);
      oP.lIdStep      = (long) dAttr(oComp, "idstep", 1.0);
      oP.nXPixels     = rR.nx;
      oP.nYPixels     = rR.ny;
      double a3dLocal[3] = {rR.xstart, rR.ystart, 0.0};
      double a3dEx[3]    = {(0.0 <= rR.xstep) ? 1.0 : -1.0, 0.0, 0.0};
      double a3dEy[3]    = {0.0, (0.0 <= rR.ystep) ? 1.0 : -1.0, 0.0};
      vMatVec(a3x3dR, a3dLocal, oP.a3dOrigin);
      for (int i = 0; i < 3; i++)
        oP.a3dOrigin[i] += a3dPos[i];
      vMatVec(a3x3dR, a3dEx, oP.a3dXAxis);
      vMatVec(a3x3dR, a3dEy, oP.a3dYAxis);
      oP.dXSize       = fabs(rR.xstep);
      oP.dYSize       = fabs(rR.ystep);
      oP.bCalibrated  = false;
      oP.nXBinned     = 0;
      oP.nYBinned     = 0;
      oP.llEvents     = 0;
      m_vPanels.push_back(oP);
    }

  double a3dL1[3] = {m_a3dSample[0] - m_a3dSource[0], m_a3dSample[1] - m_a3dSource[1],
                     m_a3dSample[2] - m_a3dSource[2]};
  m_dL1 = dNorm(a3dL1);
  return (m_vPanels.empty() ? 1 : 0);
}

int
CNXEventBinner::nOpen(const string& rsFile)
{
  h5u::CErrorSilencer oSilencer;
  m_sFile = rsFile;
  m_vPanels.clear();

  hid_t hFile = H5Fopen(rsFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (0 > hFile)
    {
      cout << "ERROR: cannot open " << rsFile << " as an HDF5 file.\n" << flush;
      return (1);
    }

  string sXML;
  if (   (0 != h5u::nReadString(hFile, "/entry/instrument/instrument_xml/data", &sXML))
      || ("" == sXML) )
    {
      cout << "ERROR: " << rsFile << " has no entry/instrument/instrument_xml/data\n"
           << "       (the Mantid instrument definition is needed for the pixel IDs).\n" << flush;
      H5Fclose(hFile);
      return (2);
    }
  if (0 != nParseIDF(sXML))
    {
      cout << "ERROR: no rectangular detector panels in the instrument definition.\n" << flush;
      H5Fclose(hFile);
      return (3);
    }

  (void) h5u::nReadString(hFile, "/entry/instrument/name", &m_sInstrument);
  (void) h5u::nReadString(hFile, "/entry/title", &m_sTitle);
  (void) h5u::nReadString(hFile, "/entry/run_number", &m_sRunNumber);
  (void) h5u::nReadString(hFile, "/entry/start_time", &m_sStartTime);
  (void) h5u::nReadString(hFile, "/entry/end_time", &m_sEndTime);
  (void) h5u::nReadDouble(hFile, "/entry/duration", &m_dDuration);
  (void) h5u::nReadDouble(hFile, "/entry/proton_charge", &m_dProtonCharge);

  for (size_t i = 0; i < m_vPanels.size(); i++)
    m_vPanels[i].bHasEvents = h5u::bPathExists(hFile, "/entry/" + m_vPanels[i].sName
                                                       + "_events/event_id");
  int nStat = nReadGonio((long) hFile);
  H5Fclose(hFile);
  return (nStat);
}

int
CNXEventBinner::nReadGonio(long lFile)
{
  hid_t hFile = (hid_t) lFile;
  m_vdGonioValues.clear();
  for (size_t i = 0; i < m_vsGonioNames.size(); i++)
    {
      vector<double> vdValues, vdTimes;
      string sLog = "/entry/DASlogs/" + m_vsGonioNames[i];
      if (0 != h5u::nReadDoubles(hFile, sLog + "/value", vdValues))
        {
          cout << "WARNING: no DAS log " << m_vsGonioNames[i] << "; using 0.\n" << flush;
          m_vdGonioValues.push_back(0.0);
          continue;
        }
      // Time-weighted average if times are available, else the mean
      double dSum = 0.0, dW = 0.0;
      if (   (0 == h5u::nReadDoubles(hFile, sLog + "/time", vdTimes))
          && (vdTimes.size() == vdValues.size()) && (1 < vdTimes.size())
          && (0.0 < m_dDuration) )
        {
          for (size_t j = 0; j < vdValues.size(); j++)
            {
              double dEnd = (j + 1 < vdTimes.size()) ? vdTimes[j + 1] : m_dDuration;
              double dLen = std::max(0.0, std::min(dEnd, m_dDuration) - std::max(0.0, vdTimes[j]));
              dSum += dLen * vdValues[j];
              dW   += dLen;
            }
        }
      if (0.0 >= dW)
        {
          dSum = 0.0;
          for (size_t j = 0; j < vdValues.size(); j++)
            dSum += vdValues[j];
          dW = (double) vdValues.size();
        }
      m_vdGonioValues.push_back(dSum / dW);
    }
  return (0);
}

int
CNXEventBinner::nApplyDetCal(const string& rsFile)
{
  std::ifstream oIn(rsFile.c_str());
  if (!oIn)
    {
      cout << "ERROR: cannot open " << rsFile << '\n' << flush;
      return (1);
    }
  string sLine;
  int    nUsed = 0;
  while (std::getline(oIn, sLine))
    {
      std::istringstream oLine(sLine);
      int nCode;
      if (!(oLine >> nCode))
        continue;
      if (7 == nCode)
        {
          double dL1cm, dT0;
          if (oLine >> dL1cm >> dT0)
            {
              m_dL1 = dL1cm * 0.01;
              m_dT0 = dT0;
            }
        }
      else if (5 == nCode)
        {
          int    nDet, nRows, nCols;
          double dW, dH, dDepth, dDist, a3dC[3], a3dBase[3], a3dUp[3];
          if (!(oLine >> nDet >> nRows >> nCols >> dW >> dH >> dDepth >> dDist
                      >> a3dC[0] >> a3dC[1] >> a3dC[2]
                      >> a3dBase[0] >> a3dBase[1] >> a3dBase[2]
                      >> a3dUp[0] >> a3dUp[1] >> a3dUp[2]))
            continue;
          for (size_t i = 0; i < m_vPanels.size(); i++)
            {
              tEventPanel& rP = m_vPanels[i];
              if (rP.nBank != nDet)
                continue;
              if ( (rP.nXPixels != nCols) || (rP.nYPixels != nRows) )
                cout << "WARNING: " << rP.sName << " is " << rP.nXPixels << " x " << rP.nYPixels
                     << " pixels in the instrument definition but " << nCols << " x " << nRows
                     << " in " << rsFile << ".\n" << flush;
              // Pixel (col, row) centre = C + (col + 0.5 - ncols/2) w/ncols Base
              //                             + (row + 0.5 - nrows/2) h/nrows Up
              rP.dXSize = 0.01 * dW / nCols;
              rP.dYSize = 0.01 * dH / nRows;
              double dLb = dNorm(a3dBase), dLu = dNorm(a3dUp);
              for (int k = 0; k < 3; k++)
                {
                  rP.a3dXAxis[k]  = a3dBase[k] / dLb;
                  rP.a3dYAxis[k]  = a3dUp[k] / dLu;
                  rP.a3dOrigin[k] = 0.01 * a3dC[k]
                                    + (0.5 - 0.5 * nCols) * rP.dXSize * rP.a3dXAxis[k]
                                    + (0.5 - 0.5 * nRows) * rP.dYSize * rP.a3dYAxis[k];
                }
              rP.bCalibrated = true;
              nUsed++;
            }
        }
    }
  m_sDetCal = rsFile;
  int nUncal = 0;
  for (size_t i = 0; i < m_vPanels.size(); i++)
    if (!m_vPanels[i].bCalibrated)
      nUncal++;
  cout << "Calibration " << rsFile << ": L1 " << m_dL1 << " m, T0 " << m_dT0 << " us, "
       << nUsed << " panels";
  if (0 < nUncal)
    cout << " (" << nUncal << " panels keep the instrument definition geometry)";
  cout << ".\n" << flush;
  return (0);
}

int
CNXEventBinner::nReadGonioCal(const string& rsFile)
{
  // Mantid parameter file with "chi-offset" and "goniometer-tilt"
  std::ifstream oIn(rsFile.c_str());
  if (!oIn)
    {
      cout << "ERROR: cannot open " << rsFile << '\n' << flush;
      return (1);
    }
  std::stringstream oBuf;
  oBuf << oIn.rdbuf();
  string sXML = oBuf.str();
  size_t nPos = 0;
  tTag   oTag;
  string sParam;
  while (bNextTag(sXML, nPos, oTag))
    {
      if ( ("parameter" == oTag.sName) && !oTag.bClose )
        sParam = sAttr(oTag, "name");
      else if ( ("value" == oTag.sName) && ("" != sParam) )
        {
          string sVal = sAttr(oTag, "val");
          if ("chi-offset" == sParam)
            m_dChiOffset = atof(sVal.c_str());
          else if ("goniometer-tilt" == sParam)
            {
              m_vdGonioTilt.clear();
              std::replace(sVal.begin(), sVal.end(), ',', ' ');
              std::istringstream oVals(sVal);
              double dV;
              while (oVals >> dV)
                m_vdGonioTilt.push_back(dV);
            }
        }
    }
  m_sGonioCal = rsFile;
  return (0);
}

int
CNXEventBinner::nTOFRangeFromEvents(long lFile)
{
  hid_t  hFile = (hid_t) lFile;
  double dMin = 1.0e30, dMax = -1.0e30;
  const hsize_t nChunk = 1 << 22;
  vector<float> vfTOF;

  for (size_t i = 0; i < m_vPanels.size(); i++)
    {
      if (!m_vPanels[i].bHasEvents)
        continue;
      string sPath = "/entry/" + m_vPanels[i].sName + "_events/event_time_offset";
      hid_t hData = H5Dopen2(hFile, sPath.c_str(), H5P_DEFAULT);
      if (0 > hData)
        continue;
      string sUnits;
      (void) h5u::nReadAttrString(hFile, sPath, "units", &sUnits);
      double dScale = ( ("ns" == sUnits) || ("nanosecond" == sUnits) ) ? 0.001 : 1.0;
      hid_t   hSpace = H5Dget_space(hData);
      hsize_t nTotal = (hsize_t) H5Sget_simple_extent_npoints(hSpace);
      H5Sclose(hSpace);
      for (hsize_t nOff = 0; nOff < nTotal; nOff += nChunk)
        {
          hsize_t nCount = std::min(nChunk, nTotal - nOff);
          vfTOF.resize((size_t) nCount);
          if (0 != h5u::nReadSlab1D(hData, H5T_NATIVE_FLOAT, nOff, nCount, &vfTOF[0]))
            break;
          for (size_t j = 0; j < vfTOF.size(); j++)
            {
              double t = vfTOF[j] * dScale + m_dT0;
              dMin = std::min(dMin, t);
              dMax = std::max(dMax, t);
            }
        }
      H5Dclose(hData);
    }
  if (dMin >= dMax)
    return (1);
  m_a2dTOFRange[0] = dMin;
  m_a2dTOFRange[1] = dMax * (1.0 + 1.0e-7);
  return (0);
}

int
CNXEventBinner::nMakeEdges(void)
{
  m_vdEdges.clear();
  double t0 = m_a2dTOFRange[0], t1 = m_a2dTOFRange[1];
  if ( (0.0 >= t0) || (t1 <= t0) )
    return (1);
  if (0.0 < m_dLogRatio)
    {
      int n = (int) ceil(log(t1 / t0) / log(1.0 + m_dLogRatio));
      for (int i = 0; i <= n; i++)
        m_vdEdges.push_back(t0 * pow(1.0 + m_dLogRatio, (double) i));
    }
  else if (0.0 < m_dLinWidth)
    {
      int n = (int) ceil((t1 - t0) / m_dLinWidth);
      for (int i = 0; i <= n; i++)
        m_vdEdges.push_back(t0 + i * m_dLinWidth);
    }
  else
    return (2);
  return (0);
}

int
CNXEventBinner::nBin(void)
{
  h5u::CErrorSilencer oSilencer;

  // Restrict to the requested banks and to panels with geometry and events

  vector<tEventPanel> vKeep;
  for (size_t i = 0; i < m_vPanels.size(); i++)
    {
      if (!m_vPanels[i].bHasEvents)
        continue;
      if (   !m_vnBanks.empty()
          && (m_vnBanks.end() == std::find(m_vnBanks.begin(), m_vnBanks.end(), m_vPanels[i].nBank)) )
        continue;
      vKeep.push_back(m_vPanels[i]);
    }
  m_vPanels.swap(vKeep);
  if (m_vPanels.empty())
    {
      cout << "ERROR: no detector panels with events to histogram.\n" << flush;
      return (1);
    }

  hid_t hFile = H5Fopen(m_sFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (0 > hFile)
    return (2);

  if ( (0.0 >= m_a2dTOFRange[0]) || (m_a2dTOFRange[1] <= m_a2dTOFRange[0]) )
    {
      cout << "Finding the TOF range of the events...\n" << flush;
      if (0 != nTOFRangeFromEvents((long) hFile))
        {
          cout << "ERROR: no events with a valid TOF.\n" << flush;
          H5Fclose(hFile);
          return (3);
        }
    }
  if (0 != nMakeEdges())
    {
      cout << "ERROR: invalid TOF binning.\n" << flush;
      H5Fclose(hFile);
      return (4);
    }
  const int    nChan  = (int) m_vdEdges.size() - 1;
  const double dT0    = m_vdEdges[0];
  const double dTEnd  = m_vdEdges[nChan];
  const double dLogR  = (0.0 < m_dLogRatio) ? log(1.0 + m_dLogRatio) : 0.0;
  const int    nB     = m_nPixBin;

  cout << "Histogramming " << m_vPanels.size() << " panels: TOF " << dT0 << " to " << dTEnd
       << " us in " << nChan << " channels, " << nB << " x " << nB << " pixel binning.\n" << flush;

  const hsize_t nChunk = 1 << 22;
  vector<uint32_t> vnIds;
  vector<float>    vfTOF;
  uint64_t llTotal = 0, llOutside = 0;

  for (size_t p = 0; p < m_vPanels.size(); p++)
    {
      tEventPanel& rP = m_vPanels[p];
      rP.nXBinned = rP.nXPixels / nB;
      rP.nYBinned = rP.nYPixels / nB;
      rP.vnCounts.assign((size_t) rP.nYBinned * rP.nXBinned * nChan, 0);
      rP.llEvents = 0;

      string sGroup = "/entry/" + rP.sName + "_events";
      hid_t hIds  = H5Dopen2(hFile, (sGroup + "/event_id").c_str(), H5P_DEFAULT);
      hid_t hTOF  = H5Dopen2(hFile, (sGroup + "/event_time_offset").c_str(), H5P_DEFAULT);
      if ( (0 > hIds) || (0 > hTOF) )
        {
          if (0 <= hIds) H5Dclose(hIds);
          if (0 <= hTOF) H5Dclose(hTOF);
          continue;
        }
      string sUnits;
      (void) h5u::nReadAttrString(hFile, sGroup + "/event_time_offset", "units", &sUnits);
      const double dScale = ( ("ns" == sUnits) || ("nanosecond" == sUnits) ) ? 0.001 : 1.0;
      hid_t   hSpace = H5Dget_space(hIds);
      hsize_t nTotal = (hsize_t) H5Sget_simple_extent_npoints(hSpace);
      H5Sclose(hSpace);

      for (hsize_t nOff = 0; nOff < nTotal; nOff += nChunk)
        {
          hsize_t nCount = std::min(nChunk, nTotal - nOff);
          vnIds.resize((size_t) nCount);
          vfTOF.resize((size_t) nCount);
          if (   (0 != h5u::nReadSlab1D(hIds, H5T_NATIVE_UINT32, nOff, nCount, &vnIds[0]))
              || (0 != h5u::nReadSlab1D(hTOF, H5T_NATIVE_FLOAT, nOff, nCount, &vfTOF[0])) )
            {
              cout << "ERROR reading the events of " << rP.sName << '\n' << flush;
              break;
            }
          for (size_t j = 0; j < (size_t) nCount; j++)
            {
              long lLocal = (long) vnIds[j] - rP.lIdStart;
              if (0 > lLocal)
                continue;
              long ix, iy;
              if (rP.bFillByY)
                {
                  ix = lLocal / rP.lIdStepByRow;
                  iy = (lLocal % rP.lIdStepByRow) / rP.lIdStep;
                }
              else
                {
                  iy = lLocal / rP.lIdStepByRow;
                  ix = (lLocal % rP.lIdStepByRow) / rP.lIdStep;
                }
              ix /= nB;
              iy /= nB;
              if ( (ix >= rP.nXBinned) || (iy >= rP.nYBinned) )
                continue;
              double t = vfTOF[j] * dScale + m_dT0;
              if ( (t < dT0) || (t >= dTEnd) )
                {
                  llOutside++;
                  continue;
                }
              int c = (0.0 < dLogR) ? (int) (log(t / dT0) / dLogR)
                                    : (int) ((t - dT0) / m_dLinWidth);
              if (c >= nChan)
                c = nChan - 1;
              rP.vnCounts[((size_t) iy * rP.nXBinned + ix) * nChan + c]++;
              rP.llEvents++;
            }
        }
      H5Dclose(hIds);
      H5Dclose(hTOF);
      llTotal += rP.llEvents;

      // Geometry of the binned pixels
      if (1 < nB)
        {
          for (int k = 0; k < 3; k++)
            rP.a3dOrigin[k] += 0.5 * (nB - 1) * (rP.dXSize * rP.a3dXAxis[k]
                                                + rP.dYSize * rP.a3dYAxis[k]);
          rP.dXSize *= nB;
          rP.dYSize *= nB;
        }
    }
  H5Fclose(hFile);
  cout << llTotal << " events histogrammed";
  if (0 < llOutside)
    cout << " (" << llOutside << " outside the TOF range)";
  cout << ".\n" << flush;
  return (0);
}

int
CNXEventBinner::nWrite(const string& rsFile)
{
  hid_t hFile = H5Fcreate(rsFile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (0 > hFile)
    {
      cout << "ERROR: cannot create " << rsFile << '\n' << flush;
      return (1);
    }
  const int nChan = (int) m_vdEdges.size() - 1;
  vector<double> vdCentres(nChan);
  for (int c = 0; c < nChan; c++)
    vdCentres[c] = 0.5 * (m_vdEdges[c] + m_vdEdges[c + 1]);

  hid_t hEntry = h5u::hCreateGroup(hFile, "entry", "NXentry");
  (void) h5u::nWriteString(hEntry, "definition", "NXlauetof");
  if ("" != m_sTitle)     (void) h5u::nWriteString(hEntry, "title", m_sTitle);
  if ("" != m_sRunNumber) (void) h5u::nWriteString(hEntry, "run_number", m_sRunNumber);
  if ("" != m_sStartTime) (void) h5u::nWriteString(hEntry, "start_time", m_sStartTime);
  if ("" != m_sEndTime)   (void) h5u::nWriteString(hEntry, "end_time", m_sEndTime);
  if (0.0 < m_dDuration)  (void) h5u::nWriteDouble(hEntry, "duration", m_dDuration, "second");
  if (0.0 < m_dProtonCharge)
    (void) h5u::nWriteDouble(hEntry, "proton_charge", m_dProtonCharge, "picoCoulomb");

  // Instrument: source and panels
  hid_t hInst = h5u::hCreateGroup(hEntry, "instrument", "NXinstrument");
  (void) h5u::nWriteString(hInst, "name", ("" != m_sInstrument) ? m_sInstrument : string("unknown"));
  {
    hid_t hSrc = h5u::hCreateGroup(hInst, "source", "NXsource");
    (void) h5u::nWriteDouble(hSrc, "distance", m_dL1, "m");
    hid_t hD = H5Dopen2(hSrc, "distance", H5P_DEFAULT);
    vector<double> vdSrcPos(3, 0.0);
    vdSrcPos[2] = -m_dL1;
    (void) h5u::nWriteAttrDoubles(hD, "position", vdSrcPos);
    H5Dclose(hD);
    (void) h5u::nWriteString(hSrc, "probe", "neutron");
    (void) h5u::nWriteString(hSrc, "type", "Spallation Neutron Source");
    H5Gclose(hSrc);
  }

  bool bDeflate = (0 < H5Zfilter_avail(H5Z_FILTER_DEFLATE));
  for (size_t p = 0; p < m_vPanels.size(); p++)
    {
      const tEventPanel& rP = m_vPanels[p];
      if (rP.vnCounts.empty())
        continue;
      hid_t hDet = h5u::hCreateGroup(hInst, rP.sName, "NXdetector");

      hsize_t a3nDims[3] = {(hsize_t) rP.nYBinned, (hsize_t) rP.nXBinned, (hsize_t) nChan};
      hid_t hSpace = H5Screate_simple(3, a3nDims, NULL);
      hid_t hProp  = H5Pcreate(H5P_DATASET_CREATE);
      hsize_t a3nChunk[3] = {1, (hsize_t) rP.nXBinned, (hsize_t) nChan};
      H5Pset_chunk(hProp, 3, a3nChunk);
      if (bDeflate)
        {
          H5Pset_shuffle(hProp);
          H5Pset_deflate(hProp, 4);
        }
      hid_t hData = H5Dcreate2(hDet, "data", H5T_STD_U32LE, hSpace, H5P_DEFAULT, hProp, H5P_DEFAULT);
      H5Dwrite(hData, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rP.vnCounts[0]);
      vector<string> vsAxes;
      vsAxes.push_back("y_pixel_offset");
      vsAxes.push_back("x_pixel_offset");
      vsAxes.push_back("tof");
      (void) h5u::nWriteAttrStrings(hData, "axes", vsAxes);
      (void) h5u::nWriteAttrDoubles(hData, "signal", vector<double>(1, 1.0));
      (void) h5u::nWriteAttrString(hData, "units", "counts");
      H5Dclose(hData);
      H5Pclose(hProp);
      H5Sclose(hSpace);

      vector<double> vdOrigin(rP.a3dOrigin, rP.a3dOrigin + 3);
      vector<double> vdY(rP.a3dYAxis, rP.a3dYAxis + 3);
      vector<double> vdX(rP.a3dXAxis, rP.a3dXAxis + 3);
      (void) h5u::nWriteDoubles(hDet, "origin", vdOrigin, "m");
      hid_t hO = H5Dopen2(hDet, "origin", H5P_DEFAULT);
      (void) h5u::nWriteAttrDoubles(hO, "first_pixel_position", vdOrigin);
      H5Dclose(hO);
      // As in the ESS NMX files: fast_axis is along the y index (data
      // dimension 0), slow_axis along the x index (data dimension 1).
      (void) h5u::nWriteDoubles(hDet, "fast_axis", vdY);
      hid_t hF = H5Dopen2(hDet, "fast_axis", H5P_DEFAULT);
      (void) h5u::nWriteAttrString(hF, "dim", "y_pixel_offset");
      H5Dclose(hF);
      (void) h5u::nWriteDoubles(hDet, "slow_axis", vdX);
      hid_t hS = H5Dopen2(hDet, "slow_axis", H5P_DEFAULT);
      (void) h5u::nWriteAttrString(hS, "dim", "x_pixel_offset");
      H5Dclose(hS);
      (void) h5u::nWriteDouble(hDet, "x_pixel_size", rP.dXSize, "m");
      (void) h5u::nWriteDouble(hDet, "y_pixel_size", rP.dYSize, "m");
      (void) h5u::nWriteDoubles(hDet, "time_of_flight", vdCentres, "us");
      (void) h5u::nWriteDoubles(hDet, "original_time_edges", m_vdEdges, "us");

      double a3dCentre[3];
      for (int k = 0; k < 3; k++)
        a3dCentre[k] = rP.a3dOrigin[k] + 0.5 * (rP.nXBinned - 1) * rP.dXSize * rP.a3dXAxis[k]
                                       + 0.5 * (rP.nYBinned - 1) * rP.dYSize * rP.a3dYAxis[k];
      double dDist = dNorm(a3dCentre);
      (void) h5u::nWriteDouble(hDet, "distance", dDist, "m");
      (void) h5u::nWriteDouble(hDet, "polar_angle",
                               acos(std::max(-1.0, std::min(1.0, a3dCentre[2] / dDist))) / c_dDegToRad,
                               "deg");
      (void) h5u::nWriteDouble(hDet, "azimuthal_angle",
                               atan2(a3dCentre[1], a3dCentre[0]) / c_dDegToRad, "deg");
      (void) h5u::nWriteAttrDoubles(hDet, "bank_number", vector<double>(1, (double) rP.nBank));
      (void) h5u::nWriteAttrString(hDet, "geometry",
                                   rP.bCalibrated ? "calibration file" : "instrument definition");
      H5Gclose(hDet);
    }
  H5Gclose(hInst);

  // Sample and goniometer (NXtransformations, outermost axis last in the chain)
  hid_t hSample = h5u::hCreateGroup(hEntry, "sample", "NXsample");
  (void) h5u::nWriteDoubles(hSample, "position", vector<double>(3, 0.0), "m");
  if (!m_vsGonioNames.empty())
    {
      hid_t hTrans = h5u::hCreateGroup(hSample, "transformations", "NXtransformations");
      for (size_t i = 0; i < m_vsGonioNames.size(); i++)
        {
          (void) h5u::nWriteDouble(hTrans, m_vsGonioNames[i], m_vdGonioValues[i], "deg");
          hid_t hA = H5Dopen2(hTrans, m_vsGonioNames[i].c_str(), H5P_DEFAULT);
          (void) h5u::nWriteAttrString(hA, "transformation_type", "rotation");
          (void) h5u::nWriteAttrDoubles(hA, "vector",
                                        vector<double>(m_vdGonioAxes.begin() + 3 * i,
                                                       m_vdGonioAxes.begin() + 3 * i + 3));
          (void) h5u::nWriteAttrString(hA, "depends_on",
                                       (0 == i) ? string(".")
                                                : "/entry/sample/transformations/" + m_vsGonioNames[i - 1]);
          H5Dclose(hA);
        }
      H5Gclose(hTrans);
      (void) h5u::nWriteString(hSample, "depends_on",
                               "/entry/sample/transformations/" + m_vsGonioNames.back());
    }
  if ("" != m_sGonioCal)
    {
      hid_t hCal = h5u::hCreateGroup(hSample, "goniometer_calibration", "NXcollection");
      (void) h5u::nWriteString(hCal, "file", m_sGonioCal);
      (void) h5u::nWriteDouble(hCal, "chi_offset", m_dChiOffset, "deg");
      if (9 == m_vdGonioTilt.size())
        (void) h5u::nWriteDoubles(hCal, "goniometer_tilt", m_vdGonioTilt);
      (void) h5u::nWriteString(hCal, "note", "Recorded, not applied to the goniometer angles");
      H5Gclose(hCal);
    }
  H5Gclose(hSample);

  // How the file was made
  hid_t hProg = h5u::hCreateGroup(hEntry, "reducer", "NXprogram");
  (void) h5u::nWriteString(hProg, "program", "d*TREK dteventbin");
  (void) h5u::nWriteString(hProg, "event_file", m_sFile);
  if ("" != m_sDetCal)
    (void) h5u::nWriteString(hProg, "calibration_file", m_sDetCal);
  (void) h5u::nWriteDouble(hProg, "t0_shift", m_dT0, "us");
  (void) h5u::nWriteDouble(hProg, "pixel_binning", (double) m_nPixBin);
  if (0.0 < m_dLogRatio)
    (void) h5u::nWriteDouble(hProg, "tof_log_bin_ratio", m_dLogRatio);
  else
    (void) h5u::nWriteDouble(hProg, "tof_linear_bin_width", m_dLinWidth, "us");
  H5Gclose(hProg);

  H5Gclose(hEntry);
  H5Fclose(hFile);
  return (0);
}

void
CNXEventBinner::vList(void) const
{
  int nCal = 0, nEv = 0;
  for (size_t i = 0; i < m_vPanels.size(); i++)
    {
      nCal += m_vPanels[i].bCalibrated ? 1 : 0;
      nEv  += m_vPanels[i].bHasEvents ? 1 : 0;
    }
  cout << "Instrument " << ("" != m_sInstrument ? m_sInstrument : string("?"))
       << ", run " << m_sRunNumber << " (" << m_sTitle << ")\n"
       << "  Panels: " << m_vPanels.size() << " (" << nEv << " with events, "
       << nCal << " calibrated)\n"
       << "  L1: " << m_dL1 << " m, T0: " << m_dT0 << " us\n"
       << "  Goniometer:";
  for (size_t i = 0; i < m_vsGonioNames.size(); i++)
    cout << ' ' << m_vsGonioNames[i] << ' '
         << ( (i < m_vdGonioValues.size()) ? m_vdGonioValues[i] : 0.0 );
  cout << '\n' << flush;
}
