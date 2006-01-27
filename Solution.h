#if !defined(SOLUTION_H_INCLUDED)
#define SOLUTION_H_INCLUDED

#include "NumKeyword.h"
//#include "Parser.h"
#include "Conc.h"
#include "Isotope.h"
#include "Pe_Data.h"
#include <cassert> // assert
#include <map>     // std::map
#include <string>  // std::string
#include <list>    // std::list
#include <vector>  // std::vector


class cxxSolution : public cxxNumKeyword
{

public:
	cxxSolution();
	cxxSolution(struct solution *);
	//cxxSolution(const cxxSolution&);
	~cxxSolution();

	//static cxxSolution& read(CParser& parser);

	void add(cxxConc conc)       { this->totals.push_back(conc); }
	void add(cxxIsotope isotope) { this->isotopes.push_back(isotope); }

	double get_ph()const {return this->ph;}
	void set_ph(double pH) {this->ph = pH;}

	double get_solution_pe()const {return this->solution_pe;}
	void set_solution_pe(double solution_pe) {this->solution_pe = solution_pe;}

	double get_tc()const {return this->tc;}
	void set_tc(double tc) {this->tc = tc;}

	double get_density()const {return this->density;}
	void set_density(double density) {this->density = density;}

	char * get_units()const {return this->units;}
	void set_units(char * units) {this->units = units;}

	//char * get_redox()const {return this->pe[this->default_pe].get_name();}

	long double get_mass_water()const {return this->mass_water;};
	void set_mass_water(long double mass_water) {this->mass_water = mass_water;};

	void dump_xml(std::ostream& os, unsigned int indent = 0)const;

protected:
	friend class cxxConc; // for this->pe access
	int new_def;
	double tc;
	double ph;
	double solution_pe;
	double mu;
	double ah2o;
	double density;
	double total_h;
	double total_o;
	double cb;
	double mass_water;
	double total_alkalinity;
	double total_co2;
	char * units;
	std::vector<cxxConc> totals; /// std::set<cxxConc> m_totals; ////std::list<cxxConc> m_totals;
	//std::vector<cxxPe_Data> pe;
	std::map <char *, struct reaction *> pe;
	int default_pe;
	std::list<cxxIsotope> isotopes;
	std::map <char *, double> species_gamma;
	std::map <char *, double> master_activity;

public:
	//static std::map<int, cxxSolution>& map;

};

#endif // !defined(SOLUTION_H_INCLUDED)
