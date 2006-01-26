#if !defined(CONC_H_INCLUDED)
#define CONC_H_INCLUDED

#include "Parser.h"
#include "utilities.h"

#include <string>

// forward declarations
class CSolution;  // reqd for read and dump_xml

class CConc
{
public:
	CConc(void);
	~CConc(void);

	enum STATUS_TYPE {
		ERROR  = 0,
		OK     = 1
	};

public:

	STATUS_TYPE read(CParser& parser, CSolution& sol);

	void dump_xml(const CSolution& solution, std::ostream& os, unsigned int indent = 0)const;

	double get_input_conc()const {return m_input_conc;}
	void set_input_conc(double input_conc) {m_input_conc = input_conc;}

	std::string get_equation_name()const {return m_equation_name;}
	void set_equation_name(std::string equation_name) {m_equation_name = equation_name;}

	std::string get_description()const {return m_description;}
	void set_description(std::string description) {m_description = description;}

	std::string get_units()const {return m_units;}
	void set_units(std::string units) {m_units = units;}

	int get_n_pe()const {return m_n_pe;}
	void set_n_pe(int n_pe) {m_n_pe = n_pe;}

	bool operator<(const CConc& conc)const    { return (m_description < conc.m_description); }

private:
	std::string m_description;
	// int m_skip;
	double m_moles;
	double m_input_conc;
	std::string m_units;
	std::string m_equation_name;
	// struct phase *m_phase;
	double m_phase_si;
	int m_n_pe;
	std::string m_as;
	double m_gfw;
};

#endif // CONC_H_INCLUDED
