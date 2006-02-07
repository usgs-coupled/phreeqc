#if !defined(EXCHCOMP_H_INCLUDED)
#define EXCHCOMP_H_INCLUDED

#include "NameDouble.h"
#define EXTERNAL extern
#include "global.h"
#include <cassert> // assert
#include <map>     // std::map
#include <string>  // std::string
#include <list>    // std::list
#include <vector>  // std::vector

#include "char_star.h"

class cxxExchComp 
{

public:
	enum EXCHANGE_CONSTRUCTOR {
		MAKE_LIST        = 1
	};

	cxxExchComp();
	cxxExchComp(struct exch_comp *);
	~cxxExchComp();


	struct master *get_master();

	static struct exch_comp *cxxExchComp2exch_comp(std::list<cxxExchComp>& el);

	void dump_xml(std::ostream& os, unsigned int indent = 0)const;

	void dump_raw(std::ostream& s_oss, unsigned int indent)const;

	void read_raw(CParser& parser);

protected:
 	char * formula;
	char * exchange_name;                     // 
 	double moles;
	cxxNameDouble formula_totals;
	cxxNameDouble totals; 
	double la;
	double charge_balance;
	char   *phase_name;
	double phase_proportion;
	char   *rate_name;
	double formula_z;                         // charge on formula

public:

};

#endif // !defined(EXCHCOMP_H_INCLUDED)
