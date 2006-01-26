#include "Isotope.h"
#include "utilities.h"
#include <cassert>

CIsotope::CIsotope(void)
: m_isotope_number(0.0)
, m_ratio_uncertainty_defined(false)
{
}

CIsotope::~CIsotope(void)
{
}

std::string CIsotope::get_name()const
{
	std::ostringstream oss;
	oss << m_isotope_number << m_elt_name;
	return oss.str();
}

CIsotope::STATUS CIsotope::read(CParser& parser)
{
	if ( !(parser.get_iss() >> m_isotope_number) ) {
		assert(parser.get_iss().fail());
		parser.incr_input_error();
		parser.error_msg("Expected isotope name to"
			" begin with an isotopic number.", CParser::OT_CONTINUE);
		return ERROR;
	}
	assert(parser.get_iss().good() || parser.get_iss().eof());

	// read and save element name
	std::istringstream::int_type c = parser.get_iss().peek();
	if ( c == std::char_traits<char>::eof() || !(::isupper(c)) ) {
		parser.error_msg("Expecting element name.", CParser::OT_CONTINUE);
		parser.error_msg(parser.line().c_str(), CParser::OT_CONTINUE);
		parser.incr_input_error();
		return ERROR;
	}
	assert(parser.get_iss().good() || parser.get_iss().eof());
	if ( !(parser.get_iss() >> m_elt_name) ) {
		// should never get here
		return ERROR;
	}
	assert(parser.get_iss().good() || parser.get_iss().eof());
	assert(!m_elt_name.empty() && ::isupper(m_elt_name[0]));

	// read and store isotope ratio
	if ( !(parser.get_iss() >> m_ratio) ) {
		assert(parser.get_iss().fail());
		parser.incr_input_error();
		parser.error_msg("Expected numeric value for isotope ratio.", CParser::OT_CONTINUE);
		return ERROR;
	}
	assert(parser.get_iss().good() || parser.get_iss().eof());

	// read and store isotope ratio
	m_ratio_uncertainty_defined = false;
	if ( !(parser.get_iss() >> m_ratio_uncertainty)) {
		if ( !parser.get_iss().eof() ) {
			parser.incr_input_error();
			parser.error_msg("Expected numeric value for uncertainty in isotope ratio.", CParser::OT_CONTINUE);
			return ERROR;
		}
	} else {
		m_ratio_uncertainty_defined = true;
	}
	assert(parser.get_iss().good() || parser.get_iss().eof());
	return OK;
}

void CIsotope::dump_xml(std::ostream& os, unsigned int indent)const
{
	unsigned int i;

	for(i = 0; i < indent; ++i) os << utilities::INDENT;
	os << "<isotope name=\"" << get_name() << "\" value=\"" << m_ratio << "\"";
	if ( m_ratio_uncertainty_defined /* utilities::isnan(m_ratio_uncertainty) */ ) {
        os << "/>\n";
	}
	else {
        os << ">\n";

		for(i = 0; i < indent + 1; ++i) os << utilities::INDENT;
		os << "<uncertainity_limit>" << m_ratio_uncertainty << "</uncertainity_limit>\n";

		for(i = 0; i < indent; ++i) os << utilities::INDENT;
        os << "</isotope>\n";
	}
}

bool CIsotope::operator<(const CIsotope& isotope)const
{
	int i = utilities::strcmp_nocase(m_elt_name.c_str(), isotope.m_elt_name.c_str());
	if (i != 0) return (i < 0);
	return ( m_isotope_number < isotope.m_isotope_number );
}

