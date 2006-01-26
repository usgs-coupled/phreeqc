#include "Conc.h"
#include "Solution.h"
#include "utilities.h"
#include <cassert>

CConc::CConc(void)
: m_description("")
// , m_skip(0);
, m_moles(0.0)
, m_input_conc(0.0)
, m_units("")
, m_equation_name("")
// , m_phase(0)
, m_phase_si(0.0)
, m_n_pe(-1)
, m_as("")
, m_gfw(0.0)
{
}

CConc::~CConc(void)
{
}

CConc::STATUS_TYPE CConc::read(CParser& parser, CSolution& solution)
{
	// std::string& str = parser.line(); 
	std::string str = parser.line();

	// defaults set in ctor

	// Remove space between "kg" and "solution" or "water" in units
	utilities::replace("Kg", "kg", str);
	utilities::replace("KG", "kg", str);
	while (utilities::replace("kg ", "kg", str));

	std::istream::pos_type ptr = 0;

	//
	// Read master species list for mass balance equation
	//
	std::string token;
	std::string token1;
	int count_redox_states = 0;
	CParser::TOKEN_TYPE j;
	while ( ((j = parser.copy_token(token, ptr)) == CParser::TT_UPPER ) ||
		( token[0] == '[' ) ||
		( utilities::strcmp_nocase_arg1(token.c_str(), "ph") == 0 ) ||
		( utilities::strcmp_nocase_arg1(token.c_str(), "pe") == 0 ) )
	{
		++count_redox_states;
		utilities::replace("(+", "(", token);
		if (count_redox_states > 1) token1 += " ";
		token1 += token;
	}
	if (count_redox_states == 0) {
		parser.incr_input_error();
		parser.error_msg("No element or master species given for concentration input.", CParser::OT_CONTINUE);
		return CConc::ERROR;
	}
	m_description = token1;

	// Determine if reading alkalinity, allow equivalents for units
	utilities::str_tolower(token1);
	bool alk = false;
	if (token1.find("alk") == 0) {
		alk = true;
	}

	// Read concentration
	if (!(std::istringstream(token) >> m_input_conc)) {
		std::ostringstream err;
		err << "Concentration data error for " << token1 << " in solution input.";
		parser.error_msg(err, CParser::OT_CONTINUE);
		return CConc::ERROR;
	}
	if ( (j = parser.copy_token(token, ptr)) == CParser::TT_EMPTY) return CConc::OK;

	// Read optional data
	token1 = token;
	
	// Check for units info
	if (parser.check_units(token1, alk, false, solution.get_units(), false) == CParser::OK) {
		if (parser.check_units(token1, alk, false, solution.get_units(), true) == CParser::OK) {
			m_units = token1;
			if ( (j = parser.copy_token(token, ptr)) == CParser::TT_EMPTY) return CConc::OK;
		} else {
			return CConc::ERROR;
		}
	}

	// Check for "as" followed by formula to be used for gfw
	token1 = token;
	utilities::str_tolower(token1);
	if (token1.compare("as") == 0)
	{
		parser.copy_token(token, ptr);
		m_as = token;
		if ( (j = parser.copy_token(token, ptr)) == CParser::TT_EMPTY) return CConc::OK;
	}
	// Check for "gfw" followed by gram formula weight
	else if (token1.compare("gfw") == 0)
	{
		if (parser.copy_token(token, ptr) != CParser::TT_DIGIT) {
			parser.error_msg("Expecting gram formula weight.", CParser::OT_CONTINUE);
			return CConc::ERROR;
		} else {
			parser.get_iss() >> m_gfw;
			if ( (j = parser.copy_token(token, ptr)) == CParser::TT_EMPTY) return CConc::OK;
		}
	}

	// Check for redox couple for pe
	if  ( utilities::strcmp_nocase_arg1(token.c_str(), "pe") == 0 ) {
		m_n_pe = CPe_Data::store(solution.m_pe, token);
		if ( (j = parser.copy_token(token, ptr)) == CParser::TT_EMPTY) return CConc::OK;
	} else if (token.find("/") != std::string::npos) {
		if (parser.parse_couple(token) == CParser::OK) {
			m_n_pe = CPe_Data::store(solution.m_pe, token);
			if ( (j = parser.copy_token(token, ptr)) == CParser::TT_EMPTY) return CConc::OK;
		} else {
			return CConc::ERROR;
		}
	}

	// Must have phase
	m_equation_name = token;
	if ( (j = parser.copy_token(token, ptr)) == CParser::TT_EMPTY) return CConc::OK;

	// Check for saturation index
	if (!(std::istringstream(token) >> m_phase_si))
	{
		parser.error_msg("Expected saturation index.", CParser::OT_CONTINUE);
		return CConc::ERROR;
	}
	return CConc::OK;
}


void CConc::dump_xml(const CSolution& solution, std::ostream& os, unsigned int indent)const
{
	unsigned int i;
	for(i = 0; i < indent; ++i) os << utilities::INDENT;
	os << "<conc>\n";

	for(i = 0; i < indent + 1; ++i) os << utilities::INDENT;
	os << "<element_list>" << m_description << "</element_list>\n";

	for(i = 0; i < indent + 1; ++i) os << utilities::INDENT;
	os << "<concentration>" << m_input_conc << "</concentration>\n";

	if (!m_units.empty()) {
		for(i = 0; i < indent + 1; ++i) os << utilities::INDENT;
		os << "<units>" << m_units << "</units>\n";
	}

	if ( !m_as.empty() ) {
		for(i = 0; i < indent + 1; ++i) os << utilities::INDENT;
		os << "<as>" << m_as << "</as>\n";
	}
	else if (m_gfw > 0.0) {
		for(i = 0; i < indent + 1; ++i) os << utilities::INDENT;
		os << "<gfw>" << m_gfw << "</gfw>\n";
	}
	////if (m_n_pe > 0) {
		solution.m_pe[m_n_pe].dump_xml(os, indent + 1);
	////}

	if (!m_equation_name.empty()) {
		if (utilities::strcmp_nocase_arg1(m_equation_name.c_str(), "charge") == 0)
		{
			for(i = 0; i < indent + 1; ++i) os << utilities::INDENT;
			os << "<charge/>\n";
		}
		else
		{
			for(i = 0; i < indent + 1; ++i) os << utilities::INDENT;
			os << "<phase_name>" << m_equation_name << "</phase_name>\n";

			for(i = 0; i < indent + 1; ++i) os << utilities::INDENT;
			os << "<saturation_index>" << m_phase_si << "</saturation_index>\n";
		}
	}

	for(i = 0; i < indent; ++i) os << utilities::INDENT;
	os << "</conc>\n";
}
