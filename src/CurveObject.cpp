// CurveObject.cpp: implementation of the CurveObject class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
#pragma warning(disable : 4786)	// disable truncation warning (Only used by debugger)
#endif
#include "CurveObject.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CurveObject::CurveObject()
	//
	// default constructor for cxxExchComp
	//
{
	x.clear();
	y.clear();
	int nxy, npoints, npoints_plot, prev_npoints;

	std::string id, color, symbol;
	int y_axis; 
	double line_w, symbol_size;
}


CurveObject::~CurveObject()
{
}

