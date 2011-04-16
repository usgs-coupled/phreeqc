#if !defined(CHARTOBJECT_H_INCLUDED)
#define CHARTOBJECT_H_INCLUDED
#include <vector>
#include <string>
#include <sstream>
#include "CurveObject.h"
#include "global_structures.h"

#include <float.h>

//#define NA (float) -9.9999	            /* NA = not available */
class ChartObject
{

  public:
	ChartObject();
	~ChartObject();

	std::vector<std::string> &Get_user_graph_headings()
	{
		return this->user_graph_headings;
	}
	std::string &Get_chart_title()
	{
		return this->chart_title;
	}
	std::vector<std::string> &Get_axis_titles()
	{
		return this->axis_titles;
	}
	bool Set_axis_scale(std::vector<std::string>, std::vector<int> types, std::ostringstream &);

  protected:

	bool ncurves_changed[3];
	std::vector<std::string> Symbol_vector;
	std::vector<std::string> Color_vector; 
	int nCSV_headers;
	std::vector<std::string> user_graph_headings;
	std::string chart_title;
	std::vector<std::string> axis_titles;
	double axis_scale_x[5];
	double axis_scale_y[5];
	double axis_scale_y2[5];

	int chart_type;
	int RowOffset, ColumnOffset;
	bool connect_simulations;
	bool graph_initial_solutions;
	bool show_chart;
	std::vector<CurveObject> Curves;
	struct rate user_graph;

  public:

};

#endif // !defined(CHARTOBJECT_H_INCLUDED)
