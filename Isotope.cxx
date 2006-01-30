#include "Isotope.h"
#include "Utils.h"
#define EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "phrqproto.h"
#include <cassert>
#include <strstream>                       // std::ostrstream

cxxIsotope::cxxIsotope(void)
: isotope_number(0.0)
{
}

cxxIsotope::cxxIsotope(struct isotope *isotope_ptr)
{
	isotope_number          = isotope_ptr->isotope_number;
	elt_name                = isotope_ptr->elt_name;
	total                   = isotope_ptr->total;
	ratio                   = isotope_ptr->ratio;
	ratio_uncertainty       = isotope_ptr->ratio_uncertainty;
	master                  = isotope_ptr->master;
	primary                 = isotope_ptr->primary;
	// Don't think these need to be stored
	//x_ratio_uncertainty
	//coef
}

cxxIsotope::~cxxIsotope(void)
{
}

struct isotope *cxxIsotope::list2isotope(std::list <cxxIsotope> &isolist)
	// takes a std::list of isotope structures
	// returns array of isotope structures
{
	struct isotope *iso;
	if (isolist.size() <= 0) {
		return NULL;
	} else {
		iso = (struct isotope *) PHRQ_malloc((size_t) ((isolist.size()) * sizeof(struct isotope)));
		if (iso == NULL) malloc_error();
		int i = 0;
		for (std::list <cxxIsotope>::iterator it = isolist.begin(); it != isolist.end(); ++it) {
			iso[i].isotope_number          = it->isotope_number;
			iso[i].elt_name                = it->elt_name;
			iso[i].total                   = it->total;
			iso[i].ratio                   = it->ratio;
			iso[i].ratio_uncertainty       = it->ratio_uncertainty;
			iso[i].master                  = it->master;
			iso[i].primary                 = it->primary;
			i++;
		}			
	}
	return(iso);
}

std::string cxxIsotope::get_name()const
{
	//std::ostringstream oss;
	std::ostrstream oss;
	oss << this->isotope_number << this->elt_name;
	return oss.str();
}
void cxxIsotope::dump_xml(std::ostream& os, unsigned int indent)const
{
	unsigned int i;

	for(i = 0; i < indent; ++i) os << Utilities::INDENT;
	os << "<isotope name=\"" << get_name() << "\" value=\"" << this->ratio << "\"";
	#ifdef SKIP
	if ( this->ratio_uncertainty_defined /* Utilities::isnan(this->ratio_uncertainty) */ ) {
        os << "/>\n";
	}
	else {
        os << ">\n";

		for(i = 0; i < indent + 1; ++i) os << Utilities::INDENT;
		os << "<uncertainity_limit>" << this->ratio_uncertainty << "</uncertainity_limit>\n";

		for(i = 0; i < indent; ++i) os << Utilities::INDENT;
        os << "</isotope>\n";
	}
	#endif
}
bool cxxIsotope::operator<(const cxxIsotope& isotope)const
{
        //int i = Utilities::strcmp_nocase(this->elt_name.c_str(), isotope.elt_name.c_str());
	int i = Utilities::strcmp_nocase(this->elt_name, isotope.elt_name);
	if (i != 0) return (i < 0);
	return ( this->isotope_number < isotope.isotope_number );
}
#ifdef SKIP
cxxIsotope::STATUS cxxIsotope::read(CParser& parser)
{
	if ( !(parser.get_iss() >> this->isotope_number) ) {
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
	if ( !(parser.get_iss() >> this->elt_name) ) {
		// should never get here
		return ERROR;
	}
	assert(parser.get_iss().good() || parser.get_iss().eof());
	assert(!this->elt_name.empty() && ::isupper(this->elt_name[0]));

	// read and store isotope ratio
	if ( !(parser.get_iss() >> this->ratio) ) {
		assert(parser.get_iss().fail());
		parser.incr_input_error();
		parser.error_msg("Expected numeric value for isotope ratio.", CParser::OT_CONTINUE);
		return ERROR;
	}
	assert(parser.get_iss().good() || parser.get_iss().eof());

	// read and store isotope ratio
	this->ratio_uncertainty_defined = false;
	if ( !(parser.get_iss() >> this->ratio_uncertainty)) {
		if ( !parser.get_iss().eof() ) {
			parser.incr_input_error();
			parser.error_msg("Expected numeric value for uncertainty in isotope ratio.", CParser::OT_CONTINUE);
			return ERROR;
		}
	} else {
		this->ratio_uncertainty_defined = true;
	}
	assert(parser.get_iss().good() || parser.get_iss().eof());
	return OK;
}
#endif
