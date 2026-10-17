#ifndef DT_H5UTIL_H
#define DT_H5UTIL_H
//
// h5util.h
//    Small helpers around the HDF5 C API for the NeXus readers and writers
//    of the TOF-Laue code (nxevents, the NXlauetof histogram reader).
//

//+Include files

#include <string>
#include <vector>
#include "hdf5.h"

namespace h5u {

// Silence the HDF5 error stack while probing for optional objects.
class CErrorSilencer
{
  H5E_auto2_t  m_pfFunc;
  void        *m_pvData;
public:
  CErrorSilencer();
  ~CErrorSilencer();
};

// True if every component of the path exists (links only; not attributes).
bool bPathExists(hid_t hLoc, const std::string& rsPath);

// Names of the links in a group, in name order.
int  nListGroup(hid_t hLoc, const std::string& rsPath, std::vector<std::string>& rvsNames);

// Dimensions of a dataset (empty on error).
int  nGetDims(hid_t hLoc, const std::string& rsPath, std::vector<hsize_t>& rvnDims);

// Whole datasets, converted to the requested type
int  nReadDoubles(hid_t hLoc, const std::string& rsPath, std::vector<double>& rvdValues);
int  nReadDouble(hid_t hLoc, const std::string& rsPath, double* pdValue);
int  nReadFloats(hid_t hLoc, const std::string& rsPath, std::vector<float>& rvfValues);

// Part of a 1D dataset: nCount values from nOffset, converted to hMemType.
int  nReadSlab1D(hid_t hData, const hid_t hMemType, const hsize_t nOffset,
                 const hsize_t nCount, void* pvBuffer);

// A string dataset (fixed or variable length, or int8 array); the first
// element if it is an array.
int  nReadString(hid_t hLoc, const std::string& rsPath, std::string* psValue);

// Attributes of an object (dataset or group) at rsPath
int  nReadAttrString(hid_t hLoc, const std::string& rsPath, const std::string& rsAttr,
                     std::string* psValue);
int  nReadAttrStrings(hid_t hLoc, const std::string& rsPath, const std::string& rsAttr,
                      std::vector<std::string>& rvsValues);
int  nReadAttrDoubles(hid_t hLoc, const std::string& rsPath, const std::string& rsAttr,
                      std::vector<double>& rvdValues);

// Writing
hid_t hCreateGroup(hid_t hLoc, const std::string& rsName, const std::string& rsNXClass);
int   nWriteAttrString(hid_t hObj, const std::string& rsAttr, const std::string& rsValue);
int   nWriteAttrStrings(hid_t hObj, const std::string& rsAttr,
                        const std::vector<std::string>& rvsValues);
int   nWriteAttrDoubles(hid_t hObj, const std::string& rsAttr,
                        const std::vector<double>& rvdValues);
int   nWriteDoubles(hid_t hLoc, const std::string& rsName, const std::vector<double>& rvdValues,
                    const std::string& rsUnits = "");
int   nWriteDouble(hid_t hLoc, const std::string& rsName, const double dValue,
                   const std::string& rsUnits = "");
int   nWriteString(hid_t hLoc, const std::string& rsName, const std::string& rsValue);

}  // namespace h5u

#endif   // DT_H5UTIL_H
