#if !defined(CURVEOBJECT_H_INCLUDED)
#define CURVEOBJECT_H_INCLUDED
#include <vector>
#include <string>

class CurveObject
{

public:
	CurveObject();
	~CurveObject();

protected:
	//float *x, *y;
	std::vector<float> x, y;
	int nxy, npoints, npoints_plot, prev_npoints;

	std::string id, color, symbol;
	int y_axis; 
	float line_w, symbol_size;

public:

};

#endif // !defined(CURVEOBJECT_H_INCLUDED)
