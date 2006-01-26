#if !defined(SOLUTION_H_INCLUDED)
#define SOLUTION_H_INCLUDED

#include "NumKeyword.h"
#include "Parser.h"
#include "Conc.h"
#include "Isotope.h"
#include "Pe_Data.h"
#include <cassert> // assert
#include <map>     // std::map
#include <string>  // std::string
#include <list>    // std::list
#include <vector>  // std::vector


class CSolution : public CNumKeyword
{
public:
	CSolution();
	~CSolution();

	static CSolution& read(CParser& parser);

	void add(CConc conc)       { m_totals.push_back(conc); }
	void add(CIsotope isotope) { m_isotopes.push_back(isotope); }

	double get_ph()const {return m_ph;}
	void set_ph(double pH) {m_ph = pH;}

	double get_solution_pe()const {return m_solution_pe;}
	void set_solution_pe(double solution_pe) {m_solution_pe = solution_pe;}

	double get_tc()const {return m_tc;}
	void set_tc(double tc) {m_tc = tc;}

	double get_density()const {return m_density;}
	void set_density(double density) {m_density = density;}

	std::string get_units()const {return m_units;}
	void set_units(std::string units) {m_units = units;}

	std::string get_redox()const {return m_pe[m_default_pe].get_name();}

	long double get_mass_water()const {return m_mass_water;};
	void set_mass_water(long double mass_water) {m_mass_water = mass_water;};

	void dump_xml(std::ostream& os, unsigned int indent = 0)const;

protected:
	friend class CConc; // for m_pe access
	double m_ph;
	double m_tc;
	double m_solution_pe;
	double m_mu;
	double m_ah2o;
	double m_density;
	long double m_mass_water;
	std::string m_units;


	std::vector<CConc> m_totals; /// std::set<CConc> m_totals; ////std::list<CConc> m_totals;
	std::vector<CPe_Data> m_pe;
	int m_default_pe;

	std::list<CIsotope> m_isotopes;

public:
	static std::map<int, CSolution>& s_map;

};

#endif // !defined(SOLUTION_H_INCLUDED)
