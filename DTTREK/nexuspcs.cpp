//
// nexuspcs.cpp
//    Reader for LANSCE PCS time-of-flight NeXus files (HDF5 version).
//    See nexuspcs.h for the file layout and the stitching rule.
//

//+Include files

#include <algorithm>
#include <string>
#include <string.h>
#include <sys/stat.h>
#include "hdf5.h"
#include "Dtrek.h"
#include "dtrekdefs.h"
#include "dtreksys.h"
#include "nexuspcs.h"
#include "CTOFLaue.h"

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
using std::min;
using std::max;
#endif

//+Definitions and constants

CPCSVolume* CPCSVolume::ms_poInstance = NULL;

static int s_nPCSRequestedChannel = nPCS_CHANNEL_SUM;

//+Local HDF5 helpers

namespace {

// Silence the HDF5 error stack while probing for optional datasets.

class CH5ErrorSilencer
{
  H5E_auto2_t  m_pfFunc;
  void        *m_pvData;
public:
  CH5ErrorSilencer()
    {
      H5Eget_auto2(H5E_DEFAULT, &m_pfFunc, &m_pvData);
      H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
    }
  ~CH5ErrorSilencer()
    {
      H5Eset_auto2(H5E_DEFAULT, m_pfFunc, m_pvData);
    }
};

bool
bH5PathExists(hid_t hLoc, const std::string& rsPath)
{
  // H5Lexists fails if an intermediate group is missing, so test
  // each component of the path in turn.

  size_t nPos = 0;
  while (std::string::npos != nPos)
    {
      nPos = rsPath.find('/', nPos + 1);
      std::string sPartial = rsPath.substr(0, nPos);
      if ( ("" == sPartial) || ("/" == sPartial) )
        continue;
      if (0 >= H5Lexists(hLoc, sPartial.c_str(), H5P_DEFAULT))
        return (false);
    }
  return (true);
}

int
nH5ReadDoubles(hid_t hLoc, const std::string& rsPath, std::vector<double>& rvdValues)
{
  hid_t    hData, hSpace;
  hssize_t nPoints;
  herr_t   nErr;

  rvdValues.clear();
  if (!bH5PathExists(hLoc, rsPath))
    return (1);
  hData = H5Dopen2(hLoc, rsPath.c_str(), H5P_DEFAULT);
  if (0 > hData)
    return (2);
  hSpace  = H5Dget_space(hData);
  nPoints = H5Sget_simple_extent_npoints(hSpace);
  H5Sclose(hSpace);
  if (0 >= nPoints)
    {
      H5Dclose(hData);
      return (3);
    }
  rvdValues.resize((size_t) nPoints);
  nErr = H5Dread(hData, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 &rvdValues[0]);
  H5Dclose(hData);
  if (0 > nErr)
    {
      rvdValues.clear();
      return (4);
    }
  return (0);
}

int
nH5ReadDouble(hid_t hLoc, const std::string& rsPath, double* pdValue)
{
  std::vector<double> vdTemp;
  if ( (0 != nH5ReadDoubles(hLoc, rsPath, vdTemp)) || (1 > vdTemp.size()) )
    return (1);
  *pdValue = vdTemp[0];
  return (0);
}

int
nH5ReadString(hid_t hLoc, const std::string& rsPath, Cstring* psValue)
{
  // Strings converted by h4toh5 are int8 arrays; also accept real strings.

  hid_t    hData, hType, hSpace, hMemType;
  hssize_t nPoints;
  size_t   nSize;
  herr_t   nErr = -1;
  std::vector<char> vcBuf;

  *psValue = "";
  if (!bH5PathExists(hLoc, rsPath))
    return (1);
  hData = H5Dopen2(hLoc, rsPath.c_str(), H5P_DEFAULT);
  if (0 > hData)
    return (2);
  hType   = H5Dget_type(hData);
  hSpace  = H5Dget_space(hData);
  nPoints = H5Sget_simple_extent_npoints(hSpace);
  H5Sclose(hSpace);

  if ( (H5T_STRING == H5Tget_class(hType)) && !H5Tis_variable_str(hType) )
    {
      nSize = H5Tget_size(hType) * (size_t) nPoints;
      vcBuf.assign(nSize + 1, '\0');
      hMemType = H5Tcopy(hType);
      nErr = H5Dread(hData, hMemType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vcBuf[0]);
      H5Tclose(hMemType);
    }
  else if (H5T_INTEGER == H5Tget_class(hType))
    {
      vcBuf.assign((size_t) nPoints + 1, '\0');
      nErr = H5Dread(hData, H5T_NATIVE_SCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vcBuf[0]);
    }
  H5Tclose(hType);
  H5Dclose(hData);
  if (0 > nErr)
    return (3);

  // Trim trailing NULs and blanks

  std::string sTemp(&vcBuf[0]);
  size_t nEnd = sTemp.find_last_not_of(" \t\r\n");
  sTemp = (std::string::npos == nEnd) ? "" : sTemp.substr(0, nEnd + 1);
  *psValue = sTemp.c_str();
  return (0);
}

// The original PCS files are HDF4 (*.nx.hdf).  If an HDF4 file is named
// (e.g. by a header from the original PCS d*TREK), use the HDF5 file of the
// same name converted by h4toh5 (*.h5) next to it.

Cstring
sPCSHDF5Name(const Cstring& rsFilename)
{
  unsigned char a4ucMagic[4] = {0, 0, 0, 0};
  FILE *pFile = fopen(rsFilename.string(), "rb");
  if (NULL == pFile)
    return (rsFilename);
  size_t nRead = fread(a4ucMagic, 1, 4, pFile);
  fclose(pFile);
  if (   (4 != nRead) || (0x0e != a4ucMagic[0]) || (0x03 != a4ucMagic[1])
      || (0x13 != a4ucMagic[2]) || (0x01 != a4ucMagic[3]) )
    return (rsFilename);   // Not HDF4

  std::string sName(rsFilename.string());
  const char *apcSuffix[3] = {".nx.hdf", ".hdf", ".nxs"};
  for (int i = 0; i < 3; i++)
    {
      size_t nLen = strlen(apcSuffix[i]);
      if (   (sName.size() > nLen)
          && (0 == sName.compare(sName.size() - nLen, nLen, apcSuffix[i])) )
        {
          std::string sH5 = sName.substr(0, sName.size() - nLen) + ".h5";
          if (bFileExists(Cstring(sH5.c_str())))
            return (Cstring(sH5.c_str()));
        }
    }
  cout << "ERROR: " << rsFilename << " is an HDF4 file.  Convert it to HDF5 with\n"
       << "       scripts/pcs_h4toh5.sh (h4toh5) first.\n" << flush;
  return (rsFilename);
}

long
lFileModTime(const Cstring& rsFilename)
{
  struct stat stBuf;
  if (0 != stat(rsFilename.string(), &stBuf))
    return (-1);
  return ((long) stBuf.st_mtime);
}

// Read the goniometer axes NXaxis1, NXaxis2, ... below a group.

void
vReadGonioAxes(hid_t hFile, const std::string& rsGroup,
               std::vector<Cstring>& rvsNames,
               std::vector<double>& rvdValues,
               std::vector<double>& rvdVectors)
{
  int nAxis;
  char acName[32];
  Cstring sName;
  double  dValue;
  std::vector<double> vdVec;

  rvsNames.clear();
  rvdValues.clear();
  rvdVectors.clear();
  for (nAxis = 1; nAxis <= 20; nAxis++)
    {
      snprintf(acName, sizeof(acName), "/NXaxis%d", nAxis);
      std::string sAxis = rsGroup + acName;
      if (!bH5PathExists(hFile, sAxis))
        break;
      if (0 != nH5ReadString(hFile, sAxis + "/name", &sName))
        sName = Cstring("axis") + Cstring(nAxis);
      if (0 != nH5ReadDouble(hFile, sAxis + "/current value", &dValue))
        dValue = 0.0;
      if ( (0 != nH5ReadDoubles(hFile, sAxis + "/vector", vdVec)) || (3 != vdVec.size()) )
        {
          vdVec.assign(3, 0.0);
          vdVec[0] = 1.0;
        }
      rvsNames.push_back(sName);
      rvdValues.push_back(dValue);
      rvdVectors.insert(rvdVectors.end(), vdVec.begin(), vdVec.end());
    }
}

Cstring
sJoin(const std::vector<Cstring>& rvsValues)
{
  Cstring sOut = "";
  for (size_t i = 0; i < rvsValues.size(); i++)
    {
      if (0 < i)
        sOut += ' ';
      sOut += rvsValues[i];
    }
  return (sOut);
}

Cstring
sJoin(const std::vector<double>& rvdValues, const char* pcFormat)
{
  Cstring sOut = "";
  char    acBuf[64];
  for (size_t i = 0; i < rvdValues.size(); i++)
    {
      snprintf(acBuf, sizeof(acBuf), pcFormat, rvdValues[i]);
      if (0 < i)
        sOut += ' ';
      sOut += acBuf;
    }
  return (sOut);
}

}  // end of anonymous namespace

//+Code begin

int
nPCSResponseMode(const Cstring& rsType)
{
  if ("none" == rsType)
    return (nPCS_RESPONSE_NONE);
  else if ("pixel" == rsType)
    return (nPCS_RESPONSE_PIXEL);
  else if ( ("spectrum" == rsType) || ("3" == rsType) ) // Legacy dtfind -type 3
    return (nPCS_RESPONSE_SPECTRUM);
  else if ("segment" == rsType)
    return (nPCS_RESPONSE_SEGMENT);
  else if ( ("full" == rsType) || ("204" == rsType) ) // Legacy dtintegrate -type 204
    return (nPCS_RESPONSE_FULL);
  return (-1);
}

void
vSetPCSRequestedChannel(const int nChan)
{
  s_nPCSRequestedChannel = nChan;
}

int
nGetPCSRequestedChannel(void)
{
  return (s_nPCSRequestedChannel);
}

CPCSVolume::CPCSVolume()
{
  m_sResponseFile = "";
  m_nResponseMode = nPCS_RESPONSE_NONE;
  vClear();
}

int
CPCSVolume::nSetResponse(const Cstring& rsFilename, const int nMode)
{
  if ( (nPCS_RESPONSE_NONE > nMode) || (nPCS_RESPONSE_SPECTRUM < nMode) )
    {
      cout << "ERROR in CPCSVolume: unknown response correction mode " << nMode << ".\n" << flush;
      return (1);
    }
  if ( (rsFilename != m_sResponseFile) || (nMode != m_nResponseMode) )
    {
      m_sResponseFile = rsFilename;
      m_nResponseMode = nMode;
      m_sFilename     = "";   // Force a reload with the new correction
    }
  return (0);
}

int
CPCSVolume::nApplyResponse(void)
{
  // Compute the detector response from the V/Nb file and divide the data by it.

  size_t nPix, p;
  int    nChan;
  double dMedian, dMean;
  CPCSVolume *poFlat;

  if ( ("" == m_sResponseFile) || (nPCS_RESPONSE_NONE == m_nResponseMode) )
    return (0);

  poFlat = new CPCSVolume();
  if (0 != poFlat->nLoad(m_sResponseFile))
    {
      delete poFlat;
      return (1);
    }
  if (   (poFlat->m_nDim0 != m_nDim0) || (poFlat->m_nDim1 != m_nDim1)
      || (poFlat->m_nNumChan != m_nNumChan) )
    {
      cout << "ERROR in CPCSVolume: response file " << m_sResponseFile
           << " has different dimensions from the data.\n" << flush;
      delete poFlat;
      return (2);
    }

  nPix = (size_t) m_nDim0 * m_nDim1;

  // Per-pixel sum over TOF and the pixels that are usable

  std::vector<double> vdPixSum(nPix, 0.0);
  for (nChan = 0; nChan < m_nNumChan; nChan++)
    {
      const float *pfV = &poFlat->m_vfData[nChan * nPix];
      for (p = 0; p < nPix; p++)
        vdPixSum[p] += pfV[p];
    }
  {
    std::vector<double> vdTemp;
    for (p = 0; p < nPix; p++)
      if (0.0 < vdPixSum[p])
        vdTemp.push_back(vdPixSum[p]);
    if (0 == vdTemp.size())
      {
        cout << "ERROR in CPCSVolume: response file has no counts.\n" << flush;
        delete poFlat;
        return (3);
      }
    std::nth_element(vdTemp.begin(), vdTemp.begin() + vdTemp.size() / 2, vdTemp.end());
    dMedian = vdTemp[vdTemp.size() / 2];
  }
  std::vector<char> vbValid(nPix, 0);
  dMean = 0.0;
  int nValid = 0;
  for (p = 0; p < nPix; p++)
    if (vdPixSum[p] > 0.1 * dMedian)
      {
        vbValid[p] = 1;
        dMean += vdPixSum[p];
        nValid++;
      }
  dMean /= (double) nValid;

  if (nPCS_RESPONSE_SPECTRUM == m_nResponseMode)
    {
      // Divide each channel by the V/Nb spectrum of the whole detector,
      // normalized to a mean of 1 over the channels (no pixel masking).

      std::vector<double> vdSpec(m_nNumChan, 0.0);
      double dSpecMean = 0.0;
      for (nChan = 0; nChan < m_nNumChan; nChan++)
        {
          const float *pfV = &poFlat->m_vfData[nChan * nPix];
          for (p = 0; p < nPix; p++)
            vdSpec[nChan] += pfV[p];
          dSpecMean += vdSpec[nChan];
        }
      dSpecMean /= (double) m_nNumChan;
      delete poFlat;
      for (nChan = 0; nChan < m_nNumChan; nChan++)
        {
          double dR     = (0.0 < dSpecMean) ? vdSpec[nChan] / dSpecMean : 0.0;
          float *pfData = &m_vfData[nChan * nPix];
          for (p = 0; p < nPix; p++)
            pfData[p] = (0.0 < dR) ? (float) (pfData[p] / dR) : 0.0f;
        }
      cout << "Detector response correction (incident spectrum) from "
           << m_sResponseFile << " applied.\n" << flush;
      return (0);
    }

  // Response of each pixel (and channel) normalized to a mean of 1

  std::vector<double> vdChanFactor;   // For the segment mode: [seg][chan]
  if (nPCS_RESPONSE_SEGMENT == m_nResponseMode)
    {
      int nSeg;
      int nNumSeg = max(1, m_nNumSeg);
      std::vector<double> vdSegSum((size_t) nNumSeg * m_nNumChan, 0.0);
      std::vector<double> vdAllSum(m_nNumChan, 0.0);
      for (nChan = 0; nChan < m_nNumChan; nChan++)
        {
          const float *pfV = &poFlat->m_vfData[nChan * nPix];
          for (p = 0; p < nPix; p++)
            {
              if (!vbValid[p])
                continue;
              nSeg = min(nNumSeg - 1, (int) (p / m_nDim0) / max(1, m_nSegStride));
              vdSegSum[(size_t) nSeg * m_nNumChan + nChan] += pfV[p];
              vdAllSum[nChan] += pfV[p];
            }
        }
      // Factor = (segment TOF shape) / (whole detector TOF shape)
      vdChanFactor.assign((size_t) nNumSeg * m_nNumChan, 1.0);
      double dAll = 0.0;
      for (nChan = 0; nChan < m_nNumChan; nChan++)
        dAll += vdAllSum[nChan];
      for (nSeg = 0; nSeg < nNumSeg; nSeg++)
        {
          double dSeg = 0.0;
          for (nChan = 0; nChan < m_nNumChan; nChan++)
            dSeg += vdSegSum[(size_t) nSeg * m_nNumChan + nChan];
          for (nChan = 0; nChan < m_nNumChan; nChan++)
            if ( (0.0 < dSeg) && (0.0 < vdAllSum[nChan]) )
              vdChanFactor[(size_t) nSeg * m_nNumChan + nChan]
                = (vdSegSum[(size_t) nSeg * m_nNumChan + nChan] / dSeg)
                  / (vdAllSum[nChan] / dAll);
        }
    }

  std::vector<float> vfSmooth;
  if (nPCS_RESPONSE_FULL == m_nResponseMode)
    {
      // Smooth over +-2 channels and 3 x 3 pixels, then normalize each channel
      int nC, nC0, nC1, i0, i1, j0, j1;
      vfSmooth.assign((size_t) m_nNumChan * nPix, 0.0f);
      for (nChan = 0; nChan < m_nNumChan; nChan++)
        {
          nC0 = max(0, nChan - 2);
          nC1 = min(m_nNumChan - 1, nChan + 2);
          double dChanMean = 0.0;
          int    nChanValid = 0;
          for (i1 = 0; i1 < m_nDim1; i1++)
            for (i0 = 0; i0 < m_nDim0; i0++)
              {
                p = (size_t) i1 * m_nDim0 + i0;
                if (!vbValid[p])
                  continue;
                double dSum = 0.0;
                int    nSum = 0;
                for (nC = nC0; nC <= nC1; nC++)
                  for (j1 = max(0, i1 - 1); j1 <= min(m_nDim1 - 1, i1 + 1); j1++)
                    for (j0 = max(0, i0 - 1); j0 <= min(m_nDim0 - 1, i0 + 1); j0++)
                      if (vbValid[(size_t) j1 * m_nDim0 + j0])
                        {
                          dSum += poFlat->m_vfData[nC * nPix + (size_t) j1 * m_nDim0 + j0];
                          nSum++;
                        }
                vfSmooth[nChan * nPix + p] = (float) (dSum / max(1, nSum));
                dChanMean += vfSmooth[nChan * nPix + p];
                nChanValid++;
              }
          dChanMean /= (double) max(1, nChanValid);
          if (0.0 < dChanMean)
            for (p = 0; p < nPix; p++)
              vfSmooth[nChan * nPix + p] = (float) (vfSmooth[nChan * nPix + p] / dChanMean);
        }
    }

  delete poFlat;

  // Divide the data by the response

  for (nChan = 0; nChan < m_nNumChan; nChan++)
    {
      float *pfData = &m_vfData[nChan * nPix];
      for (p = 0; p < nPix; p++)
        {
          double dR;
          if (!vbValid[p])
            {
              pfData[p] = -1.0f;   // Bad pixel
              continue;
            }
          if (nPCS_RESPONSE_FULL == m_nResponseMode)
            dR = vfSmooth[nChan * nPix + p];
          else
            {
              dR = vdPixSum[p] / dMean;
              if (nPCS_RESPONSE_SEGMENT == m_nResponseMode)
                {
                  int nSeg = min(max(1, m_nNumSeg) - 1,
                                 (int) (p / m_nDim0) / max(1, m_nSegStride));
                  dR *= vdChanFactor[(size_t) nSeg * m_nNumChan + nChan];
                }
            }
          pfData[p] = (0.0 < dR) ? (float) (pfData[p] / dR) : 0.0f;
        }
    }
  cout << "Detector response correction (mode " << m_nResponseMode << ") from "
       << m_sResponseFile << " applied; " << nPix - nValid << " pixels masked.\n" << flush;
  return (0);
}

CPCSVolume::~CPCSVolume()
{
  vClear();
}

CPCSVolume*
CPCSVolume::poGetInstance(void)
{
  if (NULL == ms_poInstance)
    ms_poInstance = new CPCSVolume();
  return (ms_poInstance);
}

bool
CPCSVolume::bIsHDFSignature(const char* pcBuffer)
{
  return (   (0 == memcmp(pcBuffer, "\211HDF\r\n\032\n", 8))     // HDF5
          || (0 == memcmp(pcBuffer, "\016\003\023\001", 4)) );  // HDF4
}

void
CPCSVolume::vClear(void)
{
  m_sFilename  = "";
  m_lFileTime  = -1;
  m_nSegStride = 0;
  m_nNumSeg    = 0;
  m_nDim0      = 0;
  m_nDim1      = 0;
  m_nNumChan   = 0;
  m_vfData.clear();
  m_sRunName   = "";
  m_vdTOFEdges.clear();
  m_dTimeBin   = 0.1;
  m_dTimeDelay = 0.0;
  m_dL1        = 0.0;
  m_a2dBeamCenter[0] = m_a2dBeamCenter[1] = 0.0;
  m_a3dDetSize[0] = m_a3dDetSize[1] = m_a3dDetSize[2] = 0.0;
  m_a2dPixelSize[0] = m_a2dPixelSize[1] = 0.0;
  m_a3dFastVec[0] = -1.0;  m_a3dFastVec[1] =  0.0; m_a3dFastVec[2] = 0.0;
  m_a3dSlowVec[0] =  0.0;  m_a3dSlowVec[1] = -1.0; m_a3dSlowVec[2] = 0.0;
  m_vsDetGonioNames.clear();
  m_vdDetGonioValues.clear();
  m_vdDetGonioVectors.clear();
  m_vsCrysGonioNames.clear();
  m_vdCrysGonioValues.clear();
  m_vdCrysGonioVectors.clear();
  m_sCrysGonioType = "Kappa";
  m_vdUnitCell.clear();
  m_vdOrientAngles.clear();
  m_nSpacegroup    = 1;
  m_dMosaicity     = 0.0;
  m_sDetType       = "cylinder";
  m_sDetModel      = "";
  m_sStartTime     = "";
  m_sEndTime       = "";
  m_sComment       = "";
  m_dMicroAmpHours = 0.0;
}

int
CPCSVolume::nLoad(const Cstring& rsFilename)
{
  hid_t   hFile, hData, hSpace;
  hsize_t a3nDims[3];
  int     nSeg, nNumSeg;
  int     nArc, nVert, nTOF;
  int     nOverlap, nCutoff, nStride;
  double  dTemp, dRatio;
  long    lModTime;
  char    acName[256];
  std::vector<double> vdTemp;
  std::vector<float>  vfSeg;

  Cstring sHDF5File = sPCSHDF5Name(rsFilename);   // HDF4 -> converted HDF5
  lModTime = lFileModTime(sHDF5File);
  if ( (rsFilename == m_sFilename) && (lModTime == m_lFileTime) && (0 < m_vfData.size()) )
    return (0);   // Already loaded

  vClear();

  CH5ErrorSilencer oSilencer;

  hFile = H5Fopen(sHDF5File.string(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (0 > hFile)
    {
      cout << "ERROR in CPCSVolume: cannot open " << rsFilename
           << " as an HDF5 file.\n" << flush;
      return (1);
    }

  // Find the run_N group

  {
    H5G_info_t stInfo;
    hsize_t    nIdx;
    ssize_t    nLen;
    if (0 <= H5Gget_info(hFile, &stInfo))
      {
        for (nIdx = 0; nIdx < stInfo.nlinks; nIdx++)
          {
            nLen = H5Lget_name_by_idx(hFile, "/", H5_INDEX_NAME, H5_ITER_INC, nIdx,
                                      acName, sizeof(acName), H5P_DEFAULT);
            if ( (0 < nLen) && (0 == strncmp(acName, "run_", 4)) )
              {
                m_sRunName = acName;
                break;
              }
          }
      }
  }
  if ("" == m_sRunName)
    {
      cout << "ERROR in CPCSVolume: no /run_N group in " << rsFilename << ".\n" << flush;
      H5Fclose(hFile);
      return (2);
    }

  std::string sRun = std::string("/") + m_sRunName.string();
  std::string sDet = sRun + "/NXinstrument/NXdetector1";

  // Count the detector segments

  nNumSeg = 0;
  for (nSeg = 1; nSeg <= 64; nSeg++)
    {
      snprintf(acName, sizeof(acName), "%s/tof_detector_%d/data", sRun.c_str(), nSeg);
      if (!bH5PathExists(hFile, acName))
        break;
      nNumSeg = nSeg;
    }
  if (0 == nNumSeg)
    {
      cout << "ERROR in CPCSVolume: no tof_detector_N/data datasets in "
           << rsFilename << ".\n" << flush;
      H5Fclose(hFile);
      return (3);
    }

  // Dimensions of one segment: [TOF][arc][vertical]

  snprintf(acName, sizeof(acName), "%s/tof_detector_1/data", sRun.c_str());
  hData  = H5Dopen2(hFile, acName, H5P_DEFAULT);
  hSpace = H5Dget_space(hData);
  if (3 != H5Sget_simple_extent_ndims(hSpace))
    {
      cout << "ERROR in CPCSVolume: segment data is not 3-dimensional.\n" << flush;
      H5Sclose(hSpace);
      H5Dclose(hData);
      H5Fclose(hFile);
      return (4);
    }
  H5Sget_simple_extent_dims(hSpace, a3nDims, NULL);
  H5Sclose(hSpace);
  H5Dclose(hData);
  nTOF  = (int) a3nDims[0];
  nArc  = (int) a3nDims[1];
  nVert = (int) a3nDims[2];

  // Segment overlap and border cutoff.  Files from 2011 give these in
  // unbinned pixels (segment size 512), later files in binned pixels,
  // and earlier files not at all (defaults: 16 and 8 binned pixels).

  dRatio = 1.0;
  if ( (0 == nH5ReadDouble(hFile, sDet + "/horizontal_segment_size_in_pixels", &dTemp))
      && (0.0 < dTemp) )
    dRatio = dTemp / (double) nArc;
  nOverlap = (16 * nArc) / 256;
  if (0 == nH5ReadDouble(hFile, sDet + "/horizontal_segment_overlap_in_pixels", &dTemp))
    nOverlap = (int) (dTemp / dRatio + 0.5);
  nCutoff  = (8 * nArc) / 256;
  if (0 == nH5ReadDouble(hFile, sDet + "/horizontal_border_cutoff_in_pixels", &dTemp))
    nCutoff = (int) (dTemp / dRatio + 0.5);
  nStride = nArc - nOverlap;
  m_nSegStride = nStride;
  m_nNumSeg    = nNumSeg;

  m_nDim0    = nVert;
  m_nDim1    = nStride * nNumSeg + nOverlap;   // 2048 for the PCS detector
  m_nNumChan = nTOF;
  if ( (0 == nH5ReadDoubles(hFile, sDet + "/detector dimensions", vdTemp))
      && (3 == vdTemp.size()) )
    {
      if ( ((int) vdTemp[0] != m_nDim0) || ((int) vdTemp[2] != m_nNumChan) )
        cout << "WARNING in CPCSVolume: detector dimensions " << vdTemp[0] << ' '
             << vdTemp[1] << ' ' << vdTemp[2] << " do not match the data.\n" << flush;
      if ((int) vdTemp[1] > m_nDim1)
        m_nDim1 = (int) vdTemp[1];
    }

  // Read and stitch the segments

  m_vfData.assign((size_t) m_nNumChan * m_nDim1 * m_nDim0, 0.0f);
  vfSeg.resize((size_t) nTOF * nArc * nVert);
  for (nSeg = 0; nSeg < nNumSeg; nSeg++)
    {
      snprintf(acName, sizeof(acName), "%s/tof_detector_%d/data", sRun.c_str(), nSeg + 1);
      hData  = H5Dopen2(hFile, acName, H5P_DEFAULT);
      hSpace = H5Dget_space(hData);
      H5Sget_simple_extent_dims(hSpace, a3nDims, NULL);
      H5Sclose(hSpace);
      if (   ((int) a3nDims[0] != nTOF) || ((int) a3nDims[1] != nArc)
          || ((int) a3nDims[2] != nVert) )
        {
          cout << "ERROR in CPCSVolume: segment " << nSeg + 1
               << " has different dimensions.\n" << flush;
          H5Dclose(hData);
          H5Fclose(hFile);
          vClear();
          return (5);
        }
      if (0 > H5Dread(hData, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vfSeg[0]))
        {
          cout << "ERROR in CPCSVolume: reading " << acName << ".\n" << flush;
          H5Dclose(hData);
          H5Fclose(hFile);
          vClear();
          return (6);
        }
      H5Dclose(hData);

      for (int nChan = 0; nChan < nTOF; nChan++)
        {
          for (int nA = nCutoff; nA < nArc - nCutoff; nA++)
            {
              int nPx1 = nSeg * nStride + (nA - nCutoff);
              if ( (0 > nPx1) || (m_nDim1 <= nPx1) )
                continue;
              memcpy(&m_vfData[((size_t) nChan * m_nDim1 + nPx1) * m_nDim0],
                     &vfSeg[((size_t) nChan * nArc + nA) * nVert],
                     nVert * sizeof(float));
            }
        }
    }

  // TOF bin edges

  snprintf(acName, sizeof(acName), "%s/tof_detector_1/tof", sRun.c_str());
  if ( (0 != nH5ReadDoubles(hFile, acName, m_vdTOFEdges))
      || ((int) m_vdTOFEdges.size() != m_nNumChan + 1) )
    {
      cout << "ERROR in CPCSVolume: missing or wrong-sized TOF bin edges.\n" << flush;
      H5Fclose(hFile);
      vClear();
      return (7);
    }

  // Instrument metadata

  (void) nH5ReadDouble(hFile, sDet + "/Time bin", &m_dTimeBin);
  (void) nH5ReadDouble(hFile, sDet + "/Time delay", &m_dTimeDelay);
  (void) nH5ReadDouble(hFile, sRun + "/NXinstrument/NXsource/L", &m_dL1);
  if ( (0 == nH5ReadDoubles(hFile, sDet + "/beam center", vdTemp)) && (2 == vdTemp.size()) )
    {
      m_a2dBeamCenter[0] = vdTemp[0];
      m_a2dBeamCenter[1] = vdTemp[1];
    }
  if ( (0 == nH5ReadDoubles(hFile, sDet + "/detector size", vdTemp)) && (3 == vdTemp.size()) )
    {
      m_a3dDetSize[0] = vdTemp[0];
      m_a3dDetSize[1] = vdTemp[1];
      m_a3dDetSize[2] = vdTemp[2];
    }
  if ( (0 == nH5ReadDoubles(hFile, sDet + "/pixel dimensions", vdTemp)) && (2 == vdTemp.size()) )
    {
      m_a2dPixelSize[0] = vdTemp[0];
      m_a2dPixelSize[1] = vdTemp[1];
    }
  if ( (0 == nH5ReadDoubles(hFile, sDet + "/fast pixel vector", vdTemp)) && (3 == vdTemp.size()) )
    for (int i = 0; i < 3; i++) m_a3dFastVec[i] = vdTemp[i];
  if ( (0 == nH5ReadDoubles(hFile, sDet + "/slow pixel vector", vdTemp)) && (3 == vdTemp.size()) )
    for (int i = 0; i < 3; i++) m_a3dSlowVec[i] = vdTemp[i];
  (void) nH5ReadString(hFile, sDet + "/type", &m_sDetType);
  (void) nH5ReadString(hFile, sDet + "/model", &m_sDetModel);

  vReadGonioAxes(hFile, sDet + "/NXDetector1Goniometer",
                 m_vsDetGonioNames, m_vdDetGonioValues, m_vdDetGonioVectors);
  vReadGonioAxes(hFile, sRun + "/NXsample/NXgoniometer",
                 m_vsCrysGonioNames, m_vdCrysGonioValues, m_vdCrysGonioVectors);
  (void) nH5ReadString(hFile, sRun + "/NXsample/NXgoniometer/type", &m_sCrysGonioType);
  (void) nH5ReadDoubles(hFile, sRun + "/NXsample/unit cell", m_vdUnitCell);
  (void) nH5ReadDoubles(hFile, sRun + "/NXsample/orientation angles", m_vdOrientAngles);
  if (0 == nH5ReadDouble(hFile, sRun + "/NXsample/space group", &dTemp))
    m_nSpacegroup = (int) dTemp;
  (void) nH5ReadDouble(hFile, sRun + "/NXsample/mosaicity", &m_dMosaicity);
  (void) nH5ReadString(hFile, sRun + "/start_time", &m_sStartTime);
  (void) nH5ReadString(hFile, sRun + "/end_time", &m_sEndTime);
  (void) nH5ReadString(hFile, sRun + "/comment", &m_sComment);
  (void) nH5ReadDouble(hFile, sRun + "/NXinstrument/scalars/microamphours", &m_dMicroAmpHours);

  H5Fclose(hFile);

  if (0 != nApplyResponse())
    {
      vClear();
      return (8);
    }

  m_sFilename = rsFilename;
  m_lFileTime = lModTime;
  return (0);
}

int
CPCSVolume::nBuildHeader(Cimage_header* poHeader, const Cstring& rsNameForTemplate,
                         const int nChan)
{
  // Write the keywords of the original PCS d*TREK headers, plus the
  // ones the stock d*TREK classes require for an image.

  Cstring sDet = "Detector1_";
  double  a10dRot[10];
  double  a4dTemp[4];
  int     a3nDim[3];
  char    acBuf[128];

  if (0 == m_vfData.size())
    return (1);

  // Image description

  (void) poHeader->nReplaceValue(Cimage_header::ms_sSize1, m_nDim0);
  (void) poHeader->nReplaceValue(Cimage_header::ms_sSize2, m_nDim1);
  (void) poHeader->nReplaceValue(Cstring(D_K_DataType), Cstring(D_K_FloatIEEE));
  (void) poHeader->nReplaceValue(Cstring(D_K_Compression), Cstring("PCSNXS"));
  (void) poHeader->nReplaceValue(Cimage_header::ms_sByteOrder,
                                 Cstring( (eCPU_little_endian == nGetByteOrderCPU())
                                          ? "little_endian" : "big_endian"));
  (void) poHeader->nReplaceValue(Cstring("DIM"), 2);
  (void) poHeader->nReplaceValue(Cstring("FILENAME"), sFileGetBasename(rsNameForTemplate));
  (void) poHeader->nReplaceValue(Cstring("SATURATED_VALUE"), Cstring("1e+30"));
  (void) poHeader->nReplaceValue(Cstring(D_K_OriginalImageFormat), Cstring("PCS_NeXus"));

  // Detector

  a3nDim[0] = m_nDim0;
  a3nDim[1] = m_nDim1;
  a3nDim[2] = m_nNumChan;
  (void) poHeader->nReplaceValue(Cstring("DETECTOR_NUMBER"), 1);
  (void) poHeader->nReplaceValue(Cstring(D_K_DetectorNames), sDet);
  (void) poHeader->nReplaceValue(Cstring(D_K_DetectorDimensions), 3, a3nDim);
  (void) poHeader->nReplaceValue(sDet + D_K_DetectorDimensions, 2, a3nDim);
  (void) poHeader->nReplaceValue(sDet + "DETECTOR_SIZE", 3, m_a3dDetSize, 4);
  (void) poHeader->nReplaceValue(sDet + "DETECTOR_TSIZE", m_nNumChan);
  (void) poHeader->nReplaceValue(sDet + "DETECTOR_DESCRIPTION", m_sDetType);
  if ("" != m_sDetModel)
    (void) poHeader->nReplaceValue(sDet + "DETECTOR_TYPE", m_sDetModel);
  snprintf(acBuf, sizeof(acBuf), "%g %g %g %g %g %g",
           m_a3dFastVec[0], m_a3dFastVec[1], m_a3dFastVec[2],
           m_a3dSlowVec[0], m_a3dSlowVec[1], m_a3dSlowVec[2]);
  (void) poHeader->nReplaceValue(sDet + "DETECTOR_VECTORS", Cstring(acBuf));

  (void) poHeader->nReplaceValue(sDet + "GONIO_DESCRIPTION", Cstring("UNKNOWN"));
  (void) poHeader->nReplaceValue(sDet + "GONIO_TYPE", Cstring("Eulerian"));
  (void) poHeader->nReplaceValue(sDet + "GONIO_NUM_VALUES", (int) m_vsDetGonioNames.size());
  (void) poHeader->nReplaceValue(sDet + "GONIO_NAMES", sJoin(m_vsDetGonioNames));
  {
    Cstring sUnits = "";
    for (size_t i = 0; i < m_vsDetGonioNames.size(); i++)
      {
        if (0 < i) sUnits += ' ';
        sUnits += ( m_vsDetGonioNames[i].contains("slide")
                   || m_vsDetGonioNames[i].contains("dist") ) ? "mm" : "deg";
      }
    (void) poHeader->nReplaceValue(sDet + "GONIO_UNITS", sUnits);
  }
  (void) poHeader->nReplaceValue(sDet + "GONIO_VALUES", sJoin(m_vdDetGonioValues, "%.4f"));
  (void) poHeader->nReplaceValue(sDet + "GONIO_VECTORS", sJoin(m_vdDetGonioVectors, "%g"));
  (void) poHeader->nReplaceValue(sDet + "NONUNF_TYPE", Cstring("None"));
  (void) poHeader->nReplaceValue(Cstring("DTREK_NONUNF_OKVALUE"), Cstring("0"));  // 0 counts are valid

  // Spatial: cylindrical detector, beam position and pixel sizes

  (void) poHeader->nReplaceValue(sDet + "SPATIAL_BEAM_POSITION", 2, m_a2dBeamCenter, 4);
  a4dTemp[0] = m_a2dBeamCenter[0];
  a4dTemp[1] = m_a2dBeamCenter[1];
  a4dTemp[2] = m_a2dPixelSize[0];
  a4dTemp[3] = m_a2dPixelSize[1];
  (void) poHeader->nReplaceValue(sDet + "SPATIAL_DISTORTION_INFO", 4, a4dTemp, 5);
  (void) poHeader->nReplaceValue(sDet + "SPATIAL_DISTORTION_TYPE", Cstring("Cylinder_spatial"));
  (void) poHeader->nReplaceValue(sDet + "SPATIAL_DISTORTION_VECTORS", Cstring("1 0 0 1"));

  // TOF channels (keyword names as in the original PCS d*TREK)

  (void) poHeader->nReplaceValue(sDet + CTOFLaue::ms_sRotationArray,
                                 sJoin(m_vdTOFEdges, "%.4f"));
  (void) poHeader->nReplaceValue(sDet + CTOFLaue::ms_sRotationBin, m_dTimeBin, 4);
  (void) poHeader->nReplaceValue(sDet + CTOFLaue::ms_sRotationDelay, m_dTimeDelay, 4);
  (void) poHeader->nReplaceValue(sDet + CTOFLaue::ms_sRotationL1, m_dL1, 4);

  // Crystal goniometer and (usually placeholder) crystal information

  (void) poHeader->nReplaceValue(Cstring("CRYSTAL_DESCRIPTION"), Cstring("unknown"));
  (void) poHeader->nReplaceValue(Cstring("CRYSTAL_GONIO_TYPE"), m_sCrysGonioType);
  (void) poHeader->nReplaceValue(Cstring("CRYSTAL_GONIO_NUM_VALUES"), (int) m_vsCrysGonioNames.size());
  (void) poHeader->nReplaceValue(Cstring("CRYSTAL_GONIO_NAMES"), sJoin(m_vsCrysGonioNames));
  {
    Cstring sUnits = "";
    for (size_t i = 0; i < m_vsCrysGonioNames.size(); i++)
      sUnits += (0 < i) ? " deg" : "deg";
    (void) poHeader->nReplaceValue(Cstring("CRYSTAL_GONIO_UNITS"), sUnits);
  }
  (void) poHeader->nReplaceValue(Cstring("CRYSTAL_GONIO_VALUES"), sJoin(m_vdCrysGonioValues, "%.4f"));
  (void) poHeader->nReplaceValue(Cstring("CRYSTAL_GONIO_VECTORS"), sJoin(m_vdCrysGonioVectors, "%g"));
  if (6 == m_vdUnitCell.size())
    (void) poHeader->nReplaceValue(Cstring("CRYSTAL_UNIT_CELL"), sJoin(m_vdUnitCell, "%.4f"));
  if (3 == m_vdOrientAngles.size())
    (void) poHeader->nReplaceValue(Cstring("CRYSTAL_ORIENT_ANGLES"), sJoin(m_vdOrientAngles, "%.4f"));
  (void) poHeader->nReplaceValue(Cstring("CRYSTAL_SPACEGROUP"), m_nSpacegroup);
  (void) poHeader->nReplaceValue(Cstring("CRYSTAL_MOSAICITY"), m_dMosaicity, 4);

  // Source: neutrons, so no polarization; the wavelength is 1 so that
  // d*TREK's "normalized" reciprocal space is absolute (1/Angstrom).

  (void) poHeader->nReplaceValue(Cstring("SOURCE_WAVELENGTH"), Cstring("1 1"));
  (void) poHeader->nReplaceValue(Cstring("SOURCE_VECTORS"), Cstring("0 0 1 0 1 0 1 0 0"));
  (void) poHeader->nReplaceValue(Cstring("SOURCE_VALUES"), Cstring("0 0"));
  (void) poHeader->nReplaceValue(Cstring("SOURCE_POLARZ"), Cstring("0.5 1 0 0"));
  (void) poHeader->nReplaceValue(Cstring("SOURCE_SPECTRAL_DISPERSION"), Cstring("0 0"));
  (void) poHeader->nReplaceValue(Cstring("SOURCE_CROSSFIRE"), Cstring("0 0 0 0"));
  (void) poHeader->nReplaceValue(Cstring("SOURCE_SIZE"), Cstring("0 0 0 0"));

  // Scan over TOF channels; the channel index is the pseudo-rotation

  for (int i = 0; i < 10; i++)
    a10dRot[i] = 0.0;
  a10dRot[0] = 0.0;
  a10dRot[1] = (double) m_nNumChan;
  a10dRot[2] = 1.0;
  a10dRot[3] = 1.0;
  (void) poHeader->nReplaceValue(Cstring("SCAN_") + D_K_Rotation, 10, a10dRot, 4);
  if (0 <= nChan)
    {
      a10dRot[0] = (double) nChan;
      a10dRot[1] = (double) (nChan + 1);
    }
  (void) poHeader->nReplaceValue(Cstring(D_K_Rotation), 10, a10dRot, 4);
  snprintf(acBuf, sizeof(acBuf), "0 1 %d", m_nNumChan);
  (void) poHeader->nReplaceValue(Cstring(D_K_ScanSeqInfo), Cstring(acBuf));
  (void) poHeader->nReplaceValue(Cstring(D_K_ScanTemplate), rsNameForTemplate);

  // Run information

  (void) poHeader->nReplaceValue(Cstring("PCS_RUN"), m_sRunName);
  if ("" != m_sStartTime)
    (void) poHeader->nReplaceValue(Cstring("PCS_START_TIME"), m_sStartTime);
  if ("" != m_sEndTime)
    (void) poHeader->nReplaceValue(Cstring("PCS_END_TIME"), m_sEndTime);
  if ("" != m_sComment)
    (void) poHeader->nReplaceValue(Cstring("PCS_COMMENT"), m_sComment);
  if (0.0 < m_dMicroAmpHours)
    {
      snprintf(acBuf, sizeof(acBuf), "%.0f", m_dMicroAmpHours);
      (void) poHeader->nReplaceValue(Cstring("PCS_MICROAMPHOURS"), Cstring(acBuf));
    }
  if (nPCS_CHANNEL_SUM == nChan)
    (void) poHeader->nReplaceValue(Cstring("PCS_IMAGE"), Cstring("sum"));
  else if (nPCS_CHANNEL_MAX == nChan)
    (void) poHeader->nReplaceValue(Cstring("PCS_IMAGE"), Cstring("max"));
  else
    (void) poHeader->nReplaceValue(Cstring("PCS_IMAGE"), nChan);
  return (0);
}

int
CPCSVolume::nGetImage(const int nChan, float* pfData)
{
  if (nPCS_CHANNEL_SUM == nChan)
    return (nGetSum(0, m_nNumChan - 1, pfData));
  else if (nPCS_CHANNEL_MAX == nChan)
    return (nGetMax(0, m_nNumChan - 1, pfData));
  else if ( (0 > nChan) || (m_nNumChan <= nChan) )
    {
      cout << "ERROR in CPCSVolume: TOF channel " << nChan << " is out of range 0 to "
           << m_nNumChan - 1 << ".\n" << flush;
      return (1);
    }
  memcpy(pfData, &m_vfData[(size_t) nChan * m_nDim0 * m_nDim1],
         (size_t) m_nDim0 * m_nDim1 * sizeof(float));
  return (0);
}

int
CPCSVolume::nGetSum(const int nChan0, const int nChan1, float* pfData)
{
  size_t nPix = (size_t) m_nDim0 * m_nDim1;
  int    nC0  = std::max(0, std::min(nChan0, nChan1));
  int    nC1  = std::min(m_nNumChan - 1, std::max(nChan0, nChan1));

  memset(pfData, 0, nPix * sizeof(float));
  for (int nChan = nC0; nChan <= nC1; nChan++)
    {
      const float *pfChan = &m_vfData[(size_t) nChan * nPix];
      for (size_t i = 0; i < nPix; i++)
        pfData[i] += pfChan[i];
    }
  return (0);
}

int
CPCSVolume::nGetMax(const int nChan0, const int nChan1, float* pfData)
{
  size_t nPix = (size_t) m_nDim0 * m_nDim1;
  int    nC0  = std::max(0, std::min(nChan0, nChan1));
  int    nC1  = std::min(m_nNumChan - 1, std::max(nChan0, nChan1));

  memcpy(pfData, &m_vfData[(size_t) nC0 * nPix], nPix * sizeof(float));
  for (int nChan = nC0 + 1; nChan <= nC1; nChan++)
    {
      const float *pfChan = &m_vfData[(size_t) nChan * nPix];
      for (size_t i = 0; i < nPix; i++)
        if (pfChan[i] > pfData[i])
          pfData[i] = pfChan[i];
    }
  return (0);
}

int
nReadPCSNexusHeader(const Cstring& rsFilename, const Cstring& rsNameForTemplate,
                    Cimage_header* poHeader)
{
  int         nStat;
  CPCSVolume *poVol = CPCSVolume::poGetInstance();

  nStat = poVol->nLoad(rsFilename);
  if (0 != nStat)
    return (nStat);
  return (poVol->nBuildHeader(poHeader, rsNameForTemplate, nGetPCSRequestedChannel()));
}

int
nReadPCSNexusData(const Cstring& rsFilename, float* pfData, Cimage_header* poHeader)
{
  int         nStat;
  int         nDim0 = 0, nDim1 = 0;
  CPCSVolume *poVol = CPCSVolume::poGetInstance();

  nStat = poVol->nLoad(rsFilename);
  if (0 != nStat)
    return (nStat);
  (void) poHeader->nGetValue(Cimage_header::ms_sSize1, &nDim0);
  (void) poHeader->nGetValue(Cimage_header::ms_sSize2, &nDim1);
  if ( (nDim0 != poVol->nGetDim0()) || (nDim1 != poVol->nGetDim1()) )
    {
      cout << "ERROR in nReadPCSNexusData: header and data dimensions differ.\n" << flush;
      return (1);
    }
  return (poVol->nGetImage(nGetPCSRequestedChannel(), pfData));
}
