// ChartObject.cpp: implementation of the ChartObject class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
#pragma warning(disable : 4786)	// disable truncation warning (Only used by debugger)
#endif

#include <iostream>
#include "ChartObject.h"
#define PHR_ISFINITE(x) _finite(x)
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ChartObject::ChartObject()
	//
	// default constructor for cxxExchComp
	//
{
	ncurves_changed[0] = ncurves_changed[1] = ncurves_changed[2] = 0;

	Color_vector.push_back("Red");
	Color_vector.push_back("Green");
	Color_vector.push_back("Blue");
	Color_vector.push_back("Orange");
	Color_vector.push_back("Magenta");
	Color_vector.push_back("Yellow");
	Color_vector.push_back("Black");

	Symbol_vector.push_back("Square");
	Symbol_vector.push_back("Diamond");
	Symbol_vector.push_back("Triangle");
	Symbol_vector.push_back("Circle");
	Symbol_vector.push_back("XCross");
	Symbol_vector.push_back("Plus");
	Symbol_vector.push_back("Star");
	Symbol_vector.push_back("TriangleDown");
	Symbol_vector.push_back("HDash");
	Symbol_vector.push_back("VDash");
	Symbol_vector.push_back("None");

	nCSV_headers = 0;

	user_graph_headings.clear();
	chart_title.clear();
	axis_titles.clear();
	
	int i;
	for (i = 0; i < 5; i++)
	{
		axis_scale_x[i] = NA;
		axis_scale_y[i] = NA;
		axis_scale_y2[i] = NA;
	}

	chart_type = 0;
	RowOffset = 0;
	ColumnOffset = 0;
	connect_simulations = false;
	graph_initial_solutions = true;
	show_chart = true;
}

ChartObject::~ChartObject()
{
	// all data cleans itself up
}

bool ChartObject::Set_axis_scale(std::vector<std::string> strings, 
								 std::vector<int> types,
								 std::ostringstream &estream)
{
	double *scale_ptr = NULL;
	std::string type;

	if (strings.size() == 0)
	{
		estream << "No axis defined for scales" << std::endl;
		return false;
	}

	std::string str = strings[0];

	// determine axis
	switch (str[0])
	{
	case 'X':
	case 'x':
		type = "x";
		scale_ptr = this->axis_scale_x;
		break;
	case 'Y':
	case 'y':
		type = "y";
		scale_ptr = this->axis_scale_y;
		break;
	case 'S':
	case 's':
		type = "sy";
		scale_ptr = this->axis_scale_y2;
		break;
	default:
		estream << "Found " << str;
		estream << ", but expect axis type \'x\', \'y\', or \'sy\'.";
		estream << std::endl;
		return false;
		break;
	}
	size_t j;
	for (j = 1; j < strings.size() && j < 5; j++)
	{
		std::string s = strings[j];
		if (s[0] == 'a' || s[0] == 'A')
		{
			scale_ptr[j - 1] = NAN;
		}
		else if (types[j] == DIGIT)
		{
			scale_ptr[j - 1] = atof(s.c_str());
		}
		else
		{
			estream << "Found " << s;
			estream << ", but expect number or 'a(uto)'.";
			estream << std::endl;
			return false;
		}
	}
	if (strings.size() == 5)
	{
		std::string s = strings[4];
		if (s[0] == 't' || s[0] == 'T')
		{
			if ((!PHR_ISFINITE(scale_ptr[0]) && scale_ptr[0] <=0) ||
				(!PHR_ISFINITE(scale_ptr[1]) && scale_ptr[1] <=0))
			{
				estream << "MIN and MAX must be > 0 for log " << type << "-scale.";
				estream << std::endl;
			}

		}
	}
	if (PHR_ISFINITE(scale_ptr[0]) && PHR_ISFINITE(scale_ptr[1]))
	{
		if (scale_ptr[0] > scale_ptr[1])
		{
			estream << "Maximum must be larger than minimum of axis_scale " << type << "-scale." << std::endl;
			estream << "Switching values for MIN and MAX. " << std::endl;
		}
	}
	return true;
}