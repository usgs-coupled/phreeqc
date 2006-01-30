#include "Conc.h"
#include "ISolution.h"
#include "Utils.h"
#include <cassert>
#define EXTERNAL extern
#include "global.h"
#include "phrqproto.h"
#include "phqalloc.h"

cxxConc::cxxConc(void)
: description(NULL)
, moles(0.0)
, input_conc(0.0)
, units(NULL)
, equation_name(NULL)
, phase_si(0.0)
, n_pe(-1)
, as(NULL)
, gfw(0.0)
	//, skip(0);
	//, phase(NULL)
{
}
cxxConc::cxxConc(struct conc *conc_ptr)
{
	description         = conc_ptr->description;
	moles               = conc_ptr->moles; 
	input_conc          = conc_ptr->input_conc;
	units               = conc_ptr->units;
	equation_name       = conc_ptr->equation_name;
	phase_si            = conc_ptr->phase_si;
	n_pe                = conc_ptr->n_pe;
	as                  = conc_ptr->as;
	gfw                 = conc_ptr->gfw;
	//skip                = conc_ptr->skip;
	//phase               = conc_ptr->phase;
}

cxxConc::~cxxConc(void)
{
}

struct conc *cxxConc::concarray(const std::map <char *, double> &totals)
	// for Solutions, not ISolutions
	// takes a map of (elt name, moles)
	// returns list of conc structures
{
	struct conc *c;
	c = (struct conc *) PHRQ_malloc((size_t) ((totals.size() + 1) * sizeof(struct conc)));
	if (c == NULL) malloc_error();
	int i = 0;
	for (std::map <char *, double>::const_iterator it = totals.begin(); it != totals.end(); ++it) {
		c[i].description         = (char *)it->first;
		c[i].moles               = it->second;
		c[i].input_conc          = it->second;
		c[i].units               = NULL;
		c[i].equation_name       = NULL;
		c[i].phase_si            = 0.0;
		c[i].n_pe                = 0;
		c[i].as                  = NULL;
		c[i].gfw                 = 0.0;
		c[i].skip                = 0;
		c[i].phase               = NULL;
		i++;
	}			
	c[i].description = NULL;
	return(c);
}

struct conc *cxxConc::concarray(const std::vector <cxxConc> &totals)
	// for ISolutions
	// takes a std::vector cxxConc structures
	// returns list of conc structures
{
	struct conc *c;
	c = (struct conc *) PHRQ_malloc((size_t) ((totals.size() + 1) * sizeof(struct conc)));
	if (c == NULL) malloc_error();
	int i = 0;
	for (std::vector<cxxConc>::const_iterator it = totals.begin(); it != totals.end(); ++it) {
		c[i].description         = it->description;
		c[i].moles               = it->moles;
		c[i].input_conc          = it->input_conc;
		c[i].units               = it->units;
		c[i].equation_name       = it->equation_name;
		c[i].phase_si            = it->phase_si;
		c[i].n_pe                = it->n_pe;
		c[i].as                  = it->as;
		c[i].gfw                 = it->gfw;
		c[i].skip                = 0;
		c[i].phase               = NULL;
		i++;
	}			
	c[i].description = NULL;
	return(c);
}

#ifdef SKIP
cxxConc::STATUS_TYPE cxxConc::read(CParser& parser, cxxISolution& solution)
{
	// std::string& str = parser.line(); 
	std::string str = parser.line();

	// defaults set in ctor

	// Remove space between "kg" and "solution" or "water" in units
	Utilities::replace("Kg", "kg", str);
	Utilities::replace("KG", "kg", str);
	while (Utilities::replace("kg ", "kg", str));

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
		( Utilities::strcmp_nocase_arg1(token.c_str(), "ph") == 0 ) ||
		( Utilities::strcmp_nocase_arg1(token.c_str(), "pe") == 0 ) )
	{
		++count_redox_states;
		Utilities::replace("(+", "(", token);
		if (count_redox_states > 1) token1 += " ";
		token1 += token;
	}
	if (count_redox_states == 0) {
		parser.incr_input_error();
		parser.error_msg("No element or master species given for concentration input.", CParser::OT_CONTINUE);
		return cxxConc::ERROR;
	}
	description = token1;

	// Determine if reading alkalinity, allow equivalents for units
	Utilities::str_tolower(token1);
	bool alk = false;
	if (token1.find("alk") == 0) {
		alk = true;
	}

	// Read concentration
	if (!(std::istringstream(token) >> this->input_conc)) {
		std::ostringstream err;
		err << "Concentration data error for " << token1 << " in solution input.";
		parser.error_msg(err, CParser::OT_CONTINUE);
		return cxxConc::ERROR;
	}
	if ( (j = parser.copy_token(token, ptr)) == CParser::TT_EMPTY) return cxxConc::OK;

	// Read optional data
	token1 = token;
	
	// Check for units info
	if (parser.check_units(token1, alk, false, solution.get_units(), false) == CParser::OK) {
		if (parser.check_units(token1, alk, false, solution.get_units(), true) == CParser::OK) {
			this->units = token1;
			if ( (j = parser.copy_token(token, ptr)) == CParser::TT_EMPTY) return cxxConc::OK;
		} else {
			return cxxConc::ERROR;
		}
	}

	// Check for "as" followed by formula to be used for gfw
	token1 = token;
	Utilities::str_tolower(token1);
	if (token1.compare("as") == 0)
	{
		parser.copy_token(token, ptr);
		this->as = token;
		if ( (j = parser.copy_token(token, ptr)) == CParser::TT_EMPTY) return cxxConc::OK;
	}
	// Check for "gfw" followed by gram formula weight
	else if (token1.compare("gfw") == 0)
	{
		if (parser.copy_token(token, ptr) != CParser::TT_DIGIT) {
			parser.error_msg("Expecting gram formula weight.", CParser::OT_CONTINUE);
			return cxxConc::ERROR;
		} else {
			parser.get_iss() >> this->gfw;
			if ( (j = parser.copy_token(token, ptr)) == CParser::TT_EMPTY) return cxxConc::OK;
		}
	}

	// Check for redox couple for pe
	if  ( Utilities::strcmp_nocase_arg1(token.c_str(), "pe") == 0 ) {
		this->n_pe = cxxPe_Data::store(solution.pe, token);
		if ( (j = parser.copy_token(token, ptr)) == CParser::TT_EMPTY) return cxxConc::OK;
	} else if (token.find("/") != std::string::npos) {
		if (parser.parse_couple(token) == CParser::OK) {
			this->n_pe = cxxPe_Data::store(solution.pe, token);
			if ( (j = parser.copy_token(token, ptr)) == CParser::TT_EMPTY) return cxxConc::OK;
		} else {
			return cxxConc::ERROR;
		}
	}

	// Must have phase
	this->equation_name = token;
	if ( (j = parser.copy_token(token, ptr)) == CParser::TT_EMPTY) return cxxConc::OK;

	// Check for saturation index
	if (!(std::istringstream(token) >> this->phase_si))
	{
		parser.error_msg("Expected saturation index.", CParser::OT_CONTINUE);
		return cxxConc::ERROR;
	}
	return cxxConc::OK;
}
#endif

#ifdef SKIP
void cxxConc::dump_xml(const cxxISolution& solution, std::ostream& os, unsigned int indent)const
{
	unsigned int i;
	for(i = 0; i < indent; ++i) os << Utilities::INDENT;
	os << "<conc>\n";

	for(i = 0; i < indent + 1; ++i) os << Utilities::INDENT;
	os << "<element_list>" << this->description << "</element_list>\n";

	for(i = 0; i < indent + 1; ++i) os << Utilities::INDENT;
	os << "<concentration>" << this->input_conc << "</concentration>\n";

	//if (!this->units.empty()) {
	//for(i = 0; i < indent + 1; ++i) os << Utilities::INDENT;
	//os << "<units>" << this->units << "</units>\n";
	//}

	//if ( !this->as.empty() ) {
	//for(i = 0; i < indent + 1; ++i) os << Utilities::INDENT;
	//os << "<as>" << this->as << "</as>\n";
	//}
	//else if (this->gfw > 0.0) {
	//	for(i = 0; i < indent + 1; ++i) os << Utilities::INDENT;
	//	os << "<gfw>" << this->gfw << "</gfw>\n";
	//}
	////if (this->n_pe > 0) {
		solution.pe[this->n_pe].dump_xml(os, indent + 1);
	////}

	//if (!this->equation_name.empty()) {
	//	if (Utilities::strcmp_nocase_arg1(this->equation_name.c_str(), "charge") == 0)
	//	{
	//		for(i = 0; i < indent + 1; ++i) os << Utilities::INDENT;
	//		os << "<charge/>\n";
	//	}
	//	else
	//	{
	//		for(i = 0; i < indent + 1; ++i) os << Utilities::INDENT;
	//		os << "<phase_name>" << this->equation_name << "</phase_name>\n";

	//			for(i = 0; i < indent + 1; ++i) os << Utilities::INDENT;
	//os << "<saturation_index>" << this->phase_si << "</saturation_index>\n";
	//}
	//}

	for(i = 0; i < indent; ++i) os << Utilities::INDENT;
	os << "</conc>\n";
}
#endif
