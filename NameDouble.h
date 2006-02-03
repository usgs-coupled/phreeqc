#if !defined(NAMEDOUBLE_H_INCLUDED)
#define NAMEDOUBLE_H_INCLUDED

#define EXTERNAL extern
#include "global.h"
#include <cassert> // assert
#include <map>     // std::map
#include <string>  // std::string
#include <list>    // std::list
#include <vector>  // std::vector

#include "char_star.h"
#include "Parser.h"
class cxxNameDouble : public std::map <char *, double, CHARSTAR_LESS> 
{

public:
	enum ND_TYPE {
		ND_ELT_MOLES     = 1,
		ND_SPECIES_LA    = 2,
		ND_SPECIES_GAMMA = 3
	};

	cxxNameDouble();
	cxxNameDouble(struct elt_list *);
	cxxNameDouble(struct conc *);
        cxxNameDouble(struct master_activity *ma, int count, cxxNameDouble::ND_TYPE);
	~cxxNameDouble();

	struct elt_list *exch_comp();

	void dump_xml(std::ostream& os, unsigned int indent = 0)const;

	void dump_raw(std::ostream& s_oss, unsigned int indent)const;

	void read_raw(CParser& parser);

	enum ND_TYPE type;

protected:
	//std::map <char *, double, CHARSTAR_LESS> totals; 


public:
	//static std::map<int, cxxNameDouble>& map;

};

#endif // !defined(NAMEDOUBLE_H_INCLUDED)
