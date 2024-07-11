//
// Copyright (c) 1996, 1995 Molecular Structure Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No. 
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// CUIimageviewprops.h    Initial author: J.W. Pflugrath           16-Dec-1995
//    This file is the header file for user interface class CUIimageviewprops.
/*
 *
 * Copyright (C) 2014 Rigaku Americas Corporation
 *                    9009 New Trails Drive
 *                    The Woodlands, TX, USA  77381
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *    * Redistributions of source code must retain the above copyright
 *      notice(s), this list of conditions and the following disclaimer.
 *    * Redistributions in binary form must reproduce the above copyright
 *      notice(s), this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *    * Neither the name of the Rigaku Americas Corporation nor the 
 *      names of its contributors may be used to endorse or promote products
 *      derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL RIGAKU AMERICAS CORPORATION BE LIABLE 
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA OR PROFITS; OR BUSINESS INTERUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 * DAMAGE.
 *
*/ 

//
// README: Portions of this file are merged at file generation
// time. Edits can be made *only* in between specified code blocks, look
// for keywords <Begin user code> and <End user code>.
//
//////////////////////////////////////////////////////////////
//
// Header file for CUIimageviewprops
//
//    This file is generated by Builder Xcessory 3.5.1.
//
//    This class is a user interface "component", as described
//    in "Object-Oriented Programming with C++ and OSF/Motif",
//    by Douglas Young, Prentice Hall, 1992. ISBN 0-13-630252-1
//
//////////////////////////////////////////////////////////////
#ifndef CUIimageviewprops_H
#define CUIimageviewprops_H
#include "UIComponent.h"
#include "CUIom2.h"
#include "CUIom8.h"
#include "CUIom9.h"

//
// Globally included information (change thru Output File Names Dialog).
//


//
// Class Specific Includes (change thru the class in Resource Editor).
//



// Begin user code block <head> 

#include <string.h>
#include "CXhelp.h"

typedef struct _tagImageProps 
{
  float fRotStart;
  float fRotEnd;
  float fExpTime;
  double fBorder;
  float fScaleMin;
  float fScaleMax;
  float fImageMin;
  float fImageMax;
  float fAvg;
  float fSD;
  float fAspectRatio;        // Pixel size aspect ratio (fast / slow)
  int   nDim[2];             // Image dimensions
  int   nOrig[2];            // Origin in image pixels of displayed image
  int   nExt[2];             // Extents in image pixels of displayed image
  int   nIntegrateBox[2];    // Integration box size in pixels
  int   nAutoRescale;        // Flag whether to auto re-scale new images or not
  int   nOrient;             // Orient flag (0 -> 7)
  int   nAspectMode;         // Flag to preserve (=1) or not pixel size aspect
  int   nRescaleMode;        // Flag to rescale or not
  int   nDisplayReflns;      // Flag to display (>0) reflections if available
  int   nDisplayPixValues;   // Flag to display (>0) pixel values on zoomed imgs
  float fSatPixValue;        // Saturated pixel value
  float fTooLowPixValue;     // Too low a pixel value. Default: -2
  Cstring sSatPixColor;      // Color to plot saturated pixels in, but must be
                             //   in the color table already.
  Cstring sTooLowPixColor;   // Color to plot too low pixels in, but must be
                             //   in the color table already. Default: yellow
} tagImageProps;

// End user code block <head> 

class CUIimageviewprops : public UIComponent
{

  public:

    CUIimageviewprops(const char *, Widget);
    CUIimageviewprops(const char *);
    virtual ~CUIimageviewprops();
    virtual void create(Widget);
    const char *const  className();
    
    // Begin user code block <public> 

    XtPointer  m_pObj;
    void (*m_prvOKCallback) (XtPointer, const tagImageProps&);
    void vSetValues(const tagImageProps&);
    void vGetValues(tagImageProps*);

    // End user code block <public> 
  protected:
    // Classes created by this class
    CUIom2* _omRescale;
    CUIom2* _omUseAspect;
    CUIom8* _omOrient;
    CUIom9* _omColor;
    
    // Widgets created by this class
    Widget _CUIimageviewprops;
    Widget _form;
    Widget _lbExpTimeO;
    Widget _lbRotEndO;
    Widget _lbRotStartO;
    Widget _lbImgSizeO;
    Widget _lbExpTime;
    Widget _lbRotEnd;
    Widget _lbRotStart;
    Widget _frame2;
    Widget _form2;
    Widget _tfReserve;
    Widget _tfAspect;
    Widget _lbReserve;
    Widget _lbRescale;
    Widget _lbAspect;
    Widget _lbUseAspect;
    Widget _lbOrient;
    Widget _frame1;
    Widget _form1;
    Widget _tfBorder;
    Widget _lbBorder;
    Widget _tfSatValue;
    Widget _lbSatColor;
    Widget _lbSatValue;
    Widget _tfBoxSize;
    Widget _lbBoxSize;
    Widget _lbBlackPix;
    Widget _lbWhitePix;
    Widget _tfBlackPix;
    Widget _tfWhitePix;
    Widget _arrowButton;
    Widget _lbMinMaxO;
    Widget _lbAvgSDO;
    Widget _lbMinMax;
    Widget _lbAvgSD;
    Widget _lbImgSize;
    
    // These virtual functions are called from the private callbacks 
    // or event handlers intended to be overridden in derived classes
    // to define actions
    
    virtual void vOKApplyCancelCB(Widget, XtPointer, XtPointer);
    virtual void vTextFieldCB(Widget, XtPointer, XtPointer);
    virtual void vArrowCB(Widget, XtPointer, XtPointer);
    
    // Begin user code block <protected> 

    tagImageProps m_tImageProps;

    void vUpdateWidgets(void);
    void vReadWidgets(void);

    // End user code block <protected> 
  private: 
    
    //
    // Default application and class resources.
    //
    static String         _defaultCUIimageviewpropsResources[];
    static UIAppDefault   _appDefaults[];
    static Boolean        _initAppDefaults;
    //
    // Callback client data.
    //
    UICallbackStruct  *_clientDataStructs;
    
    //
    // Callbacks to interface with Motif.
    //
    static void vOKApplyCancelCBCallback(Widget, XtPointer, XtPointer);
    static void vTextFieldCBCallback(Widget, XtPointer, XtPointer);
    static void vArrowCBCallback(Widget, XtPointer, XtPointer);
    
    // Begin user code block <private> 

    CXhelp *m_poXhelp;

    void vSetLabel(Widget, char *);

    // End user code block <private> 
};

// Begin user code block <tail> 
// End user code block <tail> 
#endif
