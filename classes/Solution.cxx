// Solution.cxx: implementation of the CSolutionxx class.
//
//////////////////////////////////////////////////////////////////////

#include "Solution.h"

#include <cassert>     // assert
#include <algorithm>   // std::sort 

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

static std::map<int, CSolution> ss_map;
std::map<int, CSolution>& CSolution::s_map = ss_map;

CSolution::CSolution()
: CNumKeyword()
, m_units("mMol/kgw")
, m_pe(CPe_Data::alloc())
{
	m_ph          = 7.0;
	m_tc          = 25.0;
	m_density     = 1.0;

	m_solution_pe = 4.0;
	m_mu          = 1e-7;
	m_ah2o        = 1.0;
	m_mass_water  = 1.0;
	m_default_pe  = -1;
}

CSolution::~CSolution()
{
}

CSolution& CSolution::read(CParser& parser)
{
	static std::vector<std::string> vopts;
	if (vopts.empty()) {
		vopts.reserve(11);
		vopts.push_back("temp");         // 0
		vopts.push_back("temperature");  // 1
		vopts.push_back("dens");         // 2
		vopts.push_back("density");      // 3
		vopts.push_back("units");        // 4
		vopts.push_back("redox");        // 5
		vopts.push_back("ph");           // 6
		vopts.push_back("pe");           // 7
		vopts.push_back("unit");         // 8
		vopts.push_back("isotope");      // 9
		vopts.push_back("water");        // 10
	}
	// const int count_opt_list = vopts.size();

	CSolution numkey;

	// Read solution number and description
	numkey.read_number_description(parser);

	// Malloc space for solution data
	//// g_solution_map[numkey.n_user()] = numkey;
	s_map[numkey.n_user()] = numkey;

	std::istream::pos_type ptr;
	std::istream::pos_type next_char;
	std::string token;
	CParser::TOKEN_TYPE j;
	
	CSolution& sol = s_map[numkey.n_user()];
	int default_pe = 0;

	for (;;)
	{
		int opt = parser.get_option(vopts, next_char);
		if (opt == CParser::OPTION_DEFAULT)
		{
			ptr = next_char;
			if (parser.copy_token(token, ptr) == CParser::TT_DIGIT) {
				opt = 9;
			}
		}

		switch (opt)
		{
		case CParser::OPTION_EOF:
			break;
		case CParser::OPTION_KEYWORD:
			break;
		case CParser::OPTION_ERROR:
			opt = CParser::OPTION_EOF;
			parser.error_msg("Unknown input in SOLUTION keyword.", CParser::OT_CONTINUE);
			parser.error_msg(parser.line().c_str(), CParser::OT_CONTINUE);
			break;

		case 0: // temp
		case 1: // temperature			
			if (!(parser.get_iss() >> sol.m_tc))
			{
				sol.m_tc = 25;
			}
			break;

		case 2: // dens
		case 3: // density
			parser.get_iss() >> sol.m_density;
			break;

		case 4: // units
		case 8: // unit
			if (parser.copy_token(token, next_char) == CParser::TT_EMPTY) break;
			if (parser.check_units(token, false, false, sol.m_units, true) == CParser::OK) {
				sol.m_units = token;				
			} else {
				parser.incr_input_error();
			}
			break;

		case 5: // redox
			if (parser.copy_token(token, next_char) == CParser::TT_EMPTY) break;
			if (parser.parse_couple(token) == CParser::OK) {
				default_pe = CPe_Data::store(sol.m_pe, token);
			} else {
				parser.incr_input_error();
			}
			break;

		case 6: // ph
			{
				CConc conc;
				if (conc.read(parser, sol) == CConc::ERROR) {
					parser.incr_input_error();
					break;
				}
				sol.m_ph = conc.get_input_conc();
				if (conc.get_equation_name().empty()) {
					break;
				}
				conc.set_description("H(1)");
				sol.add(conc);
			}
			break;

		case 7: // pe
			{
				CConc conc;
				if (conc.read(parser, sol) == CConc::ERROR) {
					parser.incr_input_error();
					break;
				}
				sol.m_solution_pe = conc.get_input_conc();
				if (conc.get_equation_name().empty()) {
					break;
				}
				conc.set_description("E");
				sol.add(conc);
			}
			break;

		case 9: // isotope
			{
				CIsotope isotope;
				if (isotope.read(parser) == CIsotope::OK) {
					sol.add(isotope);
				}
			}
			break;

		case 10: // water
			j = parser.copy_token(token, next_char);
			if (j == CParser::TT_EMPTY) {
				sol.m_mass_water = 1.0;
			} else if (j != CParser::TT_DIGIT) {
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for mass of water in solution.", CParser::OT_CONTINUE);
			} else {
				std::istringstream(token) >> sol.m_mass_water;
			}
			break;

		case CParser::OPTION_DEFAULT:
			{
				//  Read concentration
				CConc conc;
				if (conc.read(parser, sol) == CConc::ERROR) {
					parser.incr_input_error();
				} else {
					sol.add(conc);
				}
			}
			break;
		}
		if (opt == CParser::OPTION_EOF || opt == CParser::OPTION_KEYWORD) break;
	}
#ifdef SKIP
	//
	// Sort totals by description
	//
	std::sort(sol.m_totals.begin(), sol.m_totals.end());
#endif

	//
	// fix up default units and default pe
	//
	std::string token1;
	std::vector<CConc>::iterator iter = sol.m_totals.begin();
	for (; iter != sol.m_totals.end(); ++iter)
	{
		token = (*iter).get_description();
		utilities::str_tolower(token);
		if ((*iter).get_units().empty()) {
			(*iter).set_units(sol.m_units);
		} else {
			bool alk = false;
			if (token.find("alk") == 0) alk = true;
			token1 = (*iter).get_units();
			if (parser.check_units(token1, alk, true, sol.get_units(), true) == CParser::ERROR) {
				parser.incr_input_error();
			} else {
				(*iter).set_units(token1);
			}
		}
		if ((*iter).get_n_pe() < 0) {
			(*iter).set_n_pe(default_pe);
		}
	}
	sol.m_default_pe = default_pe;
	return sol;
}

void CSolution::dump_xml(std::ostream& os, unsigned int indent)const
{
	unsigned int i;

	for(i = 0; i < indent; ++i) os << utilities::INDENT;		
	os << "<solution>\n";

	CNumKeyword::dump_xml(os, indent);

	for(i = 0; i < indent + 1; ++i) os << utilities::INDENT;		
	os << "<temp>" << get_tc() << "</temp>" << "\n";

	for(i = 0; i < indent + 1; ++i) os << utilities::INDENT;		
	os << "<pH>" << get_ph() << "</pH>" << "\n";

	for(i = 0; i < indent + 1; ++i) os << utilities::INDENT;		
	os << "<pe>" << get_solution_pe() << "</pe>" << "\n";

	assert(m_pe.size() > 0);
	assert(m_default_pe >= 0);
	assert(m_pe.size() > (unsigned int) m_default_pe);
	m_pe[m_default_pe].dump_xml(os, indent + 1);

	for(i = 0; i < indent + 1; ++i) os << utilities::INDENT;		
	os << "<units>" << get_units() << "</units>" << "\n";

	for(i = 0; i < indent + 1; ++i) os << utilities::INDENT;		
	os << "<density>" << get_density() << "</density>" << "\n";

	// foreach conc
	if (!m_totals.empty())
	{
		for(i = 0; i < indent + 1; ++i) os << utilities::INDENT;		
		os << "<totals>\n";

		std::vector<CConc>::const_iterator iter = m_totals.begin();
		for(; iter != m_totals.end(); ++iter)
		{
			(*iter).dump_xml(*this, os, indent + 2);
		}

		for(i = 0; i < indent + 1; ++i) os << utilities::INDENT;		
		os << "</totals>\n";
	}

	// foreach isotope
	if (!m_isotopes.empty())
	{
		for(i = 0; i < indent + 1; ++i) os << utilities::INDENT;		
		os << "<isotopes>\n";

		std::list<CIsotope>::const_iterator iter = m_isotopes.begin();
		for(; iter != m_isotopes.end(); ++iter)
		{
			(*iter).dump_xml(os, indent + 2);
		}

		for(i = 0; i < indent + 1; ++i) os << utilities::INDENT;		
		os << "</isotopes>\n";
	}

	for(i = 0; i < indent + 1; ++i) os << utilities::INDENT;		
	os << "<water>" << get_mass_water() << "</water>" << "\n";

	for(i = 0; i < indent; ++i) os << utilities::INDENT;		
	os << "</solution>" << "\n";
}















