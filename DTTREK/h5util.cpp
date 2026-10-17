//
// h5util.cpp
//    Helpers around the HDF5 C API.  See h5util.h.
//

//+Include files

#include <string.h>
#include "h5util.h"

namespace h5u {

CErrorSilencer::CErrorSilencer()
{
  H5Eget_auto2(H5E_DEFAULT, &m_pfFunc, &m_pvData);
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
}

CErrorSilencer::~CErrorSilencer()
{
  H5Eset_auto2(H5E_DEFAULT, m_pfFunc, m_pvData);
}

bool
bPathExists(hid_t hLoc, const std::string& rsPath)
{
  // H5Lexists fails if an intermediate group is missing, so test each
  // component of the path in turn.

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
nListGroup(hid_t hLoc, const std::string& rsPath, std::vector<std::string>& rvsNames)
{
  H5G_info_t stInfo;
  rvsNames.clear();
  if ( ("/" != rsPath) && !bPathExists(hLoc, rsPath) )
    return (1);
  hid_t hGroup = H5Gopen2(hLoc, rsPath.c_str(), H5P_DEFAULT);
  if (0 > hGroup)
    return (2);
  if (0 > H5Gget_info(hGroup, &stInfo))
    {
      H5Gclose(hGroup);
      return (3);
    }
  for (hsize_t i = 0; i < stInfo.nlinks; i++)
    {
      ssize_t nLen = H5Lget_name_by_idx(hGroup, ".", H5_INDEX_NAME, H5_ITER_INC, i,
                                        NULL, 0, H5P_DEFAULT);
      if (0 >= nLen)
        continue;
      std::vector<char> vcName((size_t) nLen + 1, '\0');
      (void) H5Lget_name_by_idx(hGroup, ".", H5_INDEX_NAME, H5_ITER_INC, i,
                                &vcName[0], (size_t) nLen + 1, H5P_DEFAULT);
      rvsNames.push_back(std::string(&vcName[0]));
    }
  H5Gclose(hGroup);
  return (0);
}

int
nGetDims(hid_t hLoc, const std::string& rsPath, std::vector<hsize_t>& rvnDims)
{
  rvnDims.clear();
  if (!bPathExists(hLoc, rsPath))
    return (1);
  hid_t hData = H5Dopen2(hLoc, rsPath.c_str(), H5P_DEFAULT);
  if (0 > hData)
    return (2);
  hid_t hSpace = H5Dget_space(hData);
  int   nRank  = H5Sget_simple_extent_ndims(hSpace);
  if (0 < nRank)
    {
      rvnDims.resize((size_t) nRank);
      H5Sget_simple_extent_dims(hSpace, &rvnDims[0], NULL);
    }
  H5Sclose(hSpace);
  H5Dclose(hData);
  return ( (0 <= nRank) ? 0 : 3 );
}

static int
nReadWhole(hid_t hLoc, const std::string& rsPath, hid_t hMemType, size_t nElemSize,
           std::vector<char>& rvcBuf, hsize_t* pnPoints)
{
  if (!bPathExists(hLoc, rsPath))
    return (1);
  hid_t hData = H5Dopen2(hLoc, rsPath.c_str(), H5P_DEFAULT);
  if (0 > hData)
    return (2);
  hid_t    hSpace  = H5Dget_space(hData);
  hssize_t nPoints = H5Sget_simple_extent_npoints(hSpace);
  H5Sclose(hSpace);
  if (0 >= nPoints)
    {
      H5Dclose(hData);
      return (3);
    }
  rvcBuf.assign((size_t) nPoints * nElemSize, '\0');
  herr_t nErr = H5Dread(hData, hMemType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rvcBuf[0]);
  H5Dclose(hData);
  if (0 > nErr)
    return (4);
  *pnPoints = (hsize_t) nPoints;
  return (0);
}

int
nReadDoubles(hid_t hLoc, const std::string& rsPath, std::vector<double>& rvdValues)
{
  std::vector<char> vcBuf;
  hsize_t nPoints = 0;
  rvdValues.clear();
  int nStat = nReadWhole(hLoc, rsPath, H5T_NATIVE_DOUBLE, sizeof(double), vcBuf, &nPoints);
  if (0 != nStat)
    return (nStat);
  rvdValues.resize((size_t) nPoints);
  memcpy(&rvdValues[0], &vcBuf[0], (size_t) nPoints * sizeof(double));
  return (0);
}

int
nReadDouble(hid_t hLoc, const std::string& rsPath, double* pdValue)
{
  std::vector<double> vdTemp;
  if ( (0 != nReadDoubles(hLoc, rsPath, vdTemp)) || vdTemp.empty() )
    return (1);
  *pdValue = vdTemp[0];
  return (0);
}

int
nReadFloats(hid_t hLoc, const std::string& rsPath, std::vector<float>& rvfValues)
{
  rvfValues.clear();
  if (!bPathExists(hLoc, rsPath))
    return (1);
  hid_t hData = H5Dopen2(hLoc, rsPath.c_str(), H5P_DEFAULT);
  if (0 > hData)
    return (2);
  hid_t    hSpace  = H5Dget_space(hData);
  hssize_t nPoints = H5Sget_simple_extent_npoints(hSpace);
  H5Sclose(hSpace);
  if (0 >= nPoints)
    {
      H5Dclose(hData);
      return (3);
    }
  rvfValues.resize((size_t) nPoints);
  herr_t nErr = H5Dread(hData, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rvfValues[0]);
  H5Dclose(hData);
  if (0 > nErr)
    {
      rvfValues.clear();
      return (4);
    }
  return (0);
}

int
nReadSlab1D(hid_t hData, const hid_t hMemType, const hsize_t nOffset,
            const hsize_t nCount, void* pvBuffer)
{
  if (0 == nCount)
    return (0);
  hid_t hFileSpace = H5Dget_space(hData);
  hsize_t a1nOffset[1] = {nOffset};
  hsize_t a1nCount[1]  = {nCount};
  H5Sselect_hyperslab(hFileSpace, H5S_SELECT_SET, a1nOffset, NULL, a1nCount, NULL);
  hid_t hMemSpace = H5Screate_simple(1, a1nCount, NULL);
  herr_t nErr = H5Dread(hData, hMemType, hMemSpace, hFileSpace, H5P_DEFAULT, pvBuffer);
  H5Sclose(hMemSpace);
  H5Sclose(hFileSpace);
  return ( (0 > nErr) ? 1 : 0 );
}

// Read a string value from a dataset or attribute (hObj, given its type and
// space), handling fixed and variable-length strings and int8 arrays.

static int
nReadStringValue(hid_t hObj, bool bAttr, std::vector<std::string>& rvsValues)
{
  rvsValues.clear();
  hid_t    hType   = bAttr ? H5Aget_type(hObj) : H5Dget_type(hObj);
  hid_t    hSpace  = bAttr ? H5Aget_space(hObj) : H5Dget_space(hObj);
  hssize_t nPoints = H5Sget_simple_extent_npoints(hSpace);
  herr_t   nErr    = -1;
  if (0 >= nPoints)
    nPoints = 1;   // scalar

  if (H5T_STRING == H5Tget_class(hType))
    {
      if (H5Tis_variable_str(hType))
        {
          std::vector<char*> vpcBuf((size_t) nPoints, (char*) NULL);
          hid_t hMemType = H5Tcopy(H5T_C_S1);
          H5Tset_size(hMemType, H5T_VARIABLE);
          H5Tset_cset(hMemType, H5Tget_cset(hType));
          nErr = bAttr ? H5Aread(hObj, hMemType, &vpcBuf[0])
                       : H5Dread(hObj, hMemType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vpcBuf[0]);
          if (0 <= nErr)
            {
              for (size_t i = 0; i < vpcBuf.size(); i++)
                rvsValues.push_back(vpcBuf[i] ? std::string(vpcBuf[i]) : std::string(""));
              H5Dvlen_reclaim(hMemType, hSpace, H5P_DEFAULT, &vpcBuf[0]);
            }
          H5Tclose(hMemType);
        }
      else
        {
          size_t nSize = H5Tget_size(hType);
          std::vector<char> vcBuf((size_t) nPoints * nSize + 1, '\0');
          hid_t hMemType = H5Tcopy(hType);
          nErr = bAttr ? H5Aread(hObj, hMemType, &vcBuf[0])
                       : H5Dread(hObj, hMemType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vcBuf[0]);
          H5Tclose(hMemType);
          if (0 <= nErr)
            for (hssize_t i = 0; i < nPoints; i++)
              {
                std::string sTemp(&vcBuf[(size_t) i * nSize], nSize);
                size_t nEnd = sTemp.find('\0');
                if (std::string::npos != nEnd)
                  sTemp = sTemp.substr(0, nEnd);
                rvsValues.push_back(sTemp);
              }
        }
    }
  else if (H5T_INTEGER == H5Tget_class(hType))
    {
      std::vector<char> vcBuf((size_t) nPoints + 1, '\0');
      nErr = bAttr ? H5Aread(hObj, H5T_NATIVE_SCHAR, &vcBuf[0])
                   : H5Dread(hObj, H5T_NATIVE_SCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vcBuf[0]);
      if (0 <= nErr)
        rvsValues.push_back(std::string(&vcBuf[0]));
    }
  H5Sclose(hSpace);
  H5Tclose(hType);
  if ( (0 > nErr) || rvsValues.empty() )
    return (1);

  // Trim trailing blanks
  for (size_t i = 0; i < rvsValues.size(); i++)
    {
      size_t nEnd = rvsValues[i].find_last_not_of(" \t\r\n");
      rvsValues[i] = (std::string::npos == nEnd) ? "" : rvsValues[i].substr(0, nEnd + 1);
    }
  return (0);
}

int
nReadString(hid_t hLoc, const std::string& rsPath, std::string* psValue)
{
  std::vector<std::string> vsValues;
  *psValue = "";
  if (!bPathExists(hLoc, rsPath))
    return (1);
  hid_t hData = H5Dopen2(hLoc, rsPath.c_str(), H5P_DEFAULT);
  if (0 > hData)
    return (2);
  int nStat = nReadStringValue(hData, false, vsValues);
  H5Dclose(hData);
  if (0 != nStat)
    return (3);
  *psValue = vsValues[0];
  return (0);
}

static hid_t
hOpenAttr(hid_t hLoc, const std::string& rsPath, const std::string& rsAttr)
{
  if ( ("/" != rsPath) && ("." != rsPath) && !bPathExists(hLoc, rsPath) )
    return (-1);
  if (0 >= H5Aexists_by_name(hLoc, rsPath.c_str(), rsAttr.c_str(), H5P_DEFAULT))
    return (-1);
  return (H5Aopen_by_name(hLoc, rsPath.c_str(), rsAttr.c_str(), H5P_DEFAULT, H5P_DEFAULT));
}

int
nReadAttrStrings(hid_t hLoc, const std::string& rsPath, const std::string& rsAttr,
                 std::vector<std::string>& rvsValues)
{
  rvsValues.clear();
  hid_t hAttr = hOpenAttr(hLoc, rsPath, rsAttr);
  if (0 > hAttr)
    return (1);
  int nStat = nReadStringValue(hAttr, true, rvsValues);
  H5Aclose(hAttr);
  return (nStat);
}

int
nReadAttrString(hid_t hLoc, const std::string& rsPath, const std::string& rsAttr,
                std::string* psValue)
{
  std::vector<std::string> vsValues;
  *psValue = "";
  if (0 != nReadAttrStrings(hLoc, rsPath, rsAttr, vsValues))
    return (1);
  *psValue = vsValues[0];
  return (0);
}

int
nReadAttrDoubles(hid_t hLoc, const std::string& rsPath, const std::string& rsAttr,
                 std::vector<double>& rvdValues)
{
  rvdValues.clear();
  hid_t hAttr = hOpenAttr(hLoc, rsPath, rsAttr);
  if (0 > hAttr)
    return (1);
  hid_t    hSpace  = H5Aget_space(hAttr);
  hssize_t nPoints = H5Sget_simple_extent_npoints(hSpace);
  H5Sclose(hSpace);
  hid_t hType = H5Aget_type(hAttr);
  H5T_class_t eClass = H5Tget_class(hType);
  H5Tclose(hType);
  if ( (0 >= nPoints) || ( (H5T_FLOAT != eClass) && (H5T_INTEGER != eClass) ) )
    {
      H5Aclose(hAttr);
      return (2);
    }
  rvdValues.resize((size_t) nPoints);
  herr_t nErr = H5Aread(hAttr, H5T_NATIVE_DOUBLE, &rvdValues[0]);
  H5Aclose(hAttr);
  if (0 > nErr)
    {
      rvdValues.clear();
      return (3);
    }
  return (0);
}

hid_t
hCreateGroup(hid_t hLoc, const std::string& rsName, const std::string& rsNXClass)
{
  hid_t hGroup = H5Gcreate2(hLoc, rsName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if ( (0 <= hGroup) && ("" != rsNXClass) )
    (void) nWriteAttrString(hGroup, "NX_class", rsNXClass);
  return (hGroup);
}

int
nWriteAttrString(hid_t hObj, const std::string& rsAttr, const std::string& rsValue)
{
  hid_t hType = H5Tcopy(H5T_C_S1);
  H5Tset_size(hType, rsValue.size() ? rsValue.size() : 1);
  H5Tset_strpad(hType, H5T_STR_NULLTERM);
  hid_t hSpace = H5Screate(H5S_SCALAR);
  hid_t hAttr  = H5Acreate2(hObj, rsAttr.c_str(), hType, hSpace, H5P_DEFAULT, H5P_DEFAULT);
  herr_t nErr  = (0 > hAttr) ? -1 : H5Awrite(hAttr, hType, rsValue.size() ? rsValue.c_str() : "");
  if (0 <= hAttr)
    H5Aclose(hAttr);
  H5Sclose(hSpace);
  H5Tclose(hType);
  return ( (0 > nErr) ? 1 : 0 );
}

int
nWriteAttrStrings(hid_t hObj, const std::string& rsAttr, const std::vector<std::string>& rvsValues)
{
  size_t nMax = 1;
  for (size_t i = 0; i < rvsValues.size(); i++)
    nMax = (rvsValues[i].size() > nMax) ? rvsValues[i].size() : nMax;
  std::vector<char> vcBuf(rvsValues.size() * nMax, '\0');
  for (size_t i = 0; i < rvsValues.size(); i++)
    memcpy(&vcBuf[i * nMax], rvsValues[i].c_str(), rvsValues[i].size());
  hid_t hType = H5Tcopy(H5T_C_S1);
  H5Tset_size(hType, nMax);
  H5Tset_strpad(hType, H5T_STR_NULLPAD);
  hsize_t a1nDim[1] = {rvsValues.size()};
  hid_t hSpace = H5Screate_simple(1, a1nDim, NULL);
  hid_t hAttr  = H5Acreate2(hObj, rsAttr.c_str(), hType, hSpace, H5P_DEFAULT, H5P_DEFAULT);
  herr_t nErr  = (0 > hAttr) ? -1 : H5Awrite(hAttr, hType, &vcBuf[0]);
  if (0 <= hAttr)
    H5Aclose(hAttr);
  H5Sclose(hSpace);
  H5Tclose(hType);
  return ( (0 > nErr) ? 1 : 0 );
}

int
nWriteAttrDoubles(hid_t hObj, const std::string& rsAttr, const std::vector<double>& rvdValues)
{
  hsize_t a1nDim[1] = {rvdValues.size()};
  hid_t hSpace = (1 == rvdValues.size()) ? H5Screate(H5S_SCALAR)
                                          : H5Screate_simple(1, a1nDim, NULL);
  hid_t hAttr  = H5Acreate2(hObj, rsAttr.c_str(), H5T_NATIVE_DOUBLE, hSpace,
                            H5P_DEFAULT, H5P_DEFAULT);
  herr_t nErr  = (0 > hAttr) ? -1 : H5Awrite(hAttr, H5T_NATIVE_DOUBLE, &rvdValues[0]);
  if (0 <= hAttr)
    H5Aclose(hAttr);
  H5Sclose(hSpace);
  return ( (0 > nErr) ? 1 : 0 );
}

int
nWriteDoubles(hid_t hLoc, const std::string& rsName, const std::vector<double>& rvdValues,
              const std::string& rsUnits)
{
  hsize_t a1nDim[1] = {rvdValues.size()};
  hid_t hSpace = H5Screate_simple(1, a1nDim, NULL);
  hid_t hData  = H5Dcreate2(hLoc, rsName.c_str(), H5T_NATIVE_DOUBLE, hSpace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  herr_t nErr  = (0 > hData) ? -1
                             : H5Dwrite(hData, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                                        H5P_DEFAULT, &rvdValues[0]);
  if ( (0 <= hData) && ("" != rsUnits) )
    (void) nWriteAttrString(hData, "units", rsUnits);
  if (0 <= hData)
    H5Dclose(hData);
  H5Sclose(hSpace);
  return ( (0 > nErr) ? 1 : 0 );
}

int
nWriteDouble(hid_t hLoc, const std::string& rsName, const double dValue,
             const std::string& rsUnits)
{
  hid_t hSpace = H5Screate(H5S_SCALAR);
  hid_t hData  = H5Dcreate2(hLoc, rsName.c_str(), H5T_NATIVE_DOUBLE, hSpace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  herr_t nErr  = (0 > hData) ? -1
                             : H5Dwrite(hData, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                                        H5P_DEFAULT, &dValue);
  if ( (0 <= hData) && ("" != rsUnits) )
    (void) nWriteAttrString(hData, "units", rsUnits);
  if (0 <= hData)
    H5Dclose(hData);
  H5Sclose(hSpace);
  return ( (0 > nErr) ? 1 : 0 );
}

int
nWriteString(hid_t hLoc, const std::string& rsName, const std::string& rsValue)
{
  hid_t hType = H5Tcopy(H5T_C_S1);
  H5Tset_size(hType, rsValue.size() ? rsValue.size() : 1);
  H5Tset_strpad(hType, H5T_STR_NULLTERM);
  hid_t hSpace = H5Screate(H5S_SCALAR);
  hid_t hData  = H5Dcreate2(hLoc, rsName.c_str(), hType, hSpace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  herr_t nErr  = (0 > hData) ? -1
                             : H5Dwrite(hData, hType, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                        rsValue.size() ? rsValue.c_str() : "");
  if (0 <= hData)
    H5Dclose(hData);
  H5Sclose(hSpace);
  H5Tclose(hType);
  return ( (0 > nErr) ? 1 : 0 );
}

}  // namespace h5u
