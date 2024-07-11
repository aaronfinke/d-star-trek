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
///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
// Some parts of this file are derived from the example code from the book:
//
//           Object-Oriented Programming with C++ and OSF/Motif
//         by
//           Douglas Young
//           Prentice Hall, 1992
//           ISBN 0-13-630252-1	
//
//         Copyright 1991 by Prentice Hall
//         All Rights Reserved
//
//  Permission to use, copy, modify, and distribute this software for 
//  any purpose except publication and without fee is hereby granted, provided 
//  that the above copyright notice appear in all copies of the software.
///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////
// UIComponent.cc: Base class for all C++/Motif UI components
///////////////////////////////////////////////////////////////

#if (XlibSpecificationRelease >= 5)
#include "UIComponent.h"
#else
#define XLIB_ILLEGAL_ACCESS
#include "UIComponent.h"
#undef XLIB_ILLEGAL_ACCESS
#endif
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

UIComponent::UIComponent ( const char *name )
{
    _w = NULL;
    assert ( name != NULL );  // Make sure programmers provide name
    _name = XtNewString( name );
}

void UIComponent::widgetDestroyedCallback ( Widget, 
					   XtPointer clientData, 
					   XtPointer )
{
    UIComponent * obj = (UIComponent *) clientData;	
    
    obj->widgetDestroyed();
}

void UIComponent::widgetDestroyed()
{
    _w = NULL;
}

const char* const UIComponent::className()
{
    return "UIComponent";
}

void UIComponent::installDestroyHandler()
{
    assert ( _w != NULL );
    XtAddCallback ( _w, 
		   XmNdestroyCallback,
		   &UIComponent::widgetDestroyedCallback, 
		   (XtPointer) this );
}

void UIComponent::manage()
{
    assert ( _w != NULL );
    assert ( XtHasCallbacks ( _w, XmNdestroyCallback ) ==
	    XtCallbackHasSome );
    XtManageChild ( _w );
}

void UIComponent::unmanage()
{
    assert ( _w != NULL );
    XtUnmanageChild ( _w );
}


UIComponent::~UIComponent()
{
    // Make sure the widget hasn't already been destroyed
    
    if ( _w ) 
    {
	// Remove destroy callback so Xt can't call the callback
	// with a pointer to an object that has already been freed
	
	XtRemoveCallback ( _w, 
			  XmNdestroyCallback,
			  &UIComponent::widgetDestroyedCallback,
			  (XtPointer) this );	

	XtDestroyWidget ( _w );
    }

    XtFree(_name);
}

void UIComponent::getResources ( const XtResourceList resources, 
				const int numResources )
{
    // Check for errors
    
    assert ( _w != NULL );
    assert ( resources != NULL );
    
    // Retrieve the requested resources relative to the 
    // parent of this object's base widget
    
    XtGetSubresources ( XtParent( _w ), 
		       (XtPointer) this, 
		       _name,
		       className(),
		       resources, 
		       numResources,
		       NULL, 
		       0 );
}

void UIComponent::setDefaultResources ( const Widget w, 
                                         const String *resourceSpec )
{
   int         i;	
   Display    *dpy = XtDisplay ( w );	  // Retrieve the display pointer
   XrmDatabase rdb = NULL;             // A resource data base

   // Create an empty resource database

   rdb = XrmGetStringDatabase ( "" );

   // Add the Component resources, prepending the name of the component

   i = 0;
   while ( resourceSpec[i] != NULL )
   {
       char buf[1000];

       snprintf(buf, sizeof(buf), "*%s%s", _name, resourceSpec[i++]);
       XrmPutLineResource( &rdb, buf );
   }

   // Merge them into the Xt database, with lowest precendence

   if ( rdb )
   {
#if (XlibSpecificationRelease>=5)
        XrmDatabase db = XtDatabase(dpy);
	XrmCombineDatabase(rdb, &db, FALSE);
#else
        XrmMergeDatabases ( dpy->db, &rdb );
        dpy->db = rdb;
#endif
    }
}

// This method gets all the resources from the app-defaults file 
// (resource databse) and fills in the table (defs) if the app default 
// value exists.

void
UIComponent::initAppDefaults(const Widget parent,
			     const char * cname,
			     UIAppDefault * defs)
{
    XrmQuark		cQuark = XrmStringToQuark(cname);	// class quark
    XrmQuark		rsc[6];
    XrmRepresentation	rep;
    XrmValue		val;
    XrmDatabase		rdb;
    int			rscIdx;

    // Get the database

#if (XlibSpecificationRelease >= 5)
    if ((rdb = XrmGetDatabase(XtDisplay(parent))) == NULL)
    {
	return;			// Can't get the database
    }
#else
    Display *dpy = XtDisplay(parent);
    if ((rdb = dpy->db) == NULL)
    {
	return;
    }
#endif

    // Look for each resource in the table

    while (defs->wName)
    {
	rscIdx = 0;
	rsc[rscIdx++] = cQuark;
	if (defs->wName[0] == '\0')
	{
	    rsc[rscIdx++] = cQuark;
	}
	else
	{
	    rsc[rscIdx++] = XrmStringToQuark(defs->wName);
	}

	if (defs->cInstName && defs->cInstName[0] != '\0')
	{
	    rsc[rscIdx++] = XrmStringToQuark(defs->cInstName);
	}

	rsc[rscIdx++] = XrmStringToQuark(defs->wRsc);
	rsc[rscIdx++] = NULLQUARK;

	if (XrmQGetResource(rdb, rsc, rsc, &rep, &val))
	{
	    defs->value = XtNewString((char*)val.addr);
	}
	defs++;
    }
}

// This method applies the app defaults for the class to a specific
// instance. All the widgets in the path are loosly coupled (use *).
// To override a specific instance, use a tightly coupled app defaults
// resource line (use .).
void
UIComponent::setAppDefaults(const Widget w,
				UIAppDefault * defs,
				const char* inst_name,
				Boolean override_inst)
{
   Display*	dpy = XtDisplay ( w );	// Retrieve the display pointer
   XrmDatabase	rdb = NULL;		// A resource data base
   char		lineage[1000];
   char		buf[1000];
   Widget       parent;
   char 	*wName;

   // Create an empty resource database

   rdb = XrmGetStringDatabase ( "" );

   // Get the lineage

   lineage[0] = '\0';
   parent = w;

   while (parent)
   {
       WidgetClass wclass = XtClass(parent);

       if (wclass == applicationShellWidgetClass) break;

       strcpy(buf, lineage);
       snprintf(lineage, sizeof(lineage), "*%s%s", XtName(parent), buf);

       parent = XtParent(parent);
   }

   // Add the Component resources, prepending the name of the component

   while (defs->wName != NULL)
   {
       if (defs->value == NULL ||
	   override_inst && strcmp(inst_name, defs->wName) ||
	   ! override_inst && defs->cInstName != NULL)
       {
	   defs++;
	   continue;
       }

       // Build up string after lineage
       if (defs->cInstName != NULL)
       {
	   // Don't include class instance name if it is also the instance
	   // being affected. 

	   if (*defs->cInstName != '\0')
	   {
	       snprintf(buf, sizeof(buf), "%s*%s*%s.%s: %s",
		       lineage, defs->wName, defs->cInstName, defs->wRsc,
		       defs->value);
	   }
	   else
	   {
	       snprintf(buf, sizeof(buf), "%s*%s.%s: %s",
		       lineage, defs->wName, defs->wRsc, defs->value);
	   }
       }
       else if (*defs->wName != '\0')
       {
	   snprintf(buf, sizeof(buf), "%s*%s*%s.%s: %s",
		   lineage, inst_name, defs->wName, defs->wRsc, defs->value);
       }
       else
       {
	   snprintf(buf, sizeof(buf), "%s*%s.%s: %s",
		   lineage, inst_name, defs->wRsc, defs->value);
       }

       XrmPutLineResource( &rdb, buf );
       defs++;
   }

   // Merge them into the Xt database, with lowest precendence
   if ( rdb )
   {
#if (XlibSpecificationRelease >= 5)
        XrmDatabase db = XtDatabase(dpy);
	XrmCombineDatabase(rdb, &db, FALSE);
#else
        XrmMergeDatabases ( dpy->db, &rdb );
        dpy->db = rdb;
#endif
    }
}

