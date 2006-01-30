// Solution.cxx: implementation of the cxxSolution class.
//
//////////////////////////////////////////////////////////////////////
#include "Solution.h"
#include "Utilities.h"   // define before global.h
#define EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "phrqproto.h"
#include <cassert>     // assert
#include <algorithm>   // std::sort 

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//static std::map<int, cxxSolution> ss_map;
//std::map<int, cxxSolution>& cxxSolution::s_map = ss_map;

cxxSolution::cxxSolution()
	//
	// default constructor for cxxSolution 
	//
: cxxNumKeyword()
{
	tc          = 25.0;
	ph          = 7.0;
	pe          = 4.0;
	mu          = 1e-7;
	ah2o        = 1.0;
	total_h     = 111.1;
	total_o     = 55.55;
	cb          = 0.0;
	mass_water  = 1.0;
	total_alkalinity = 0.0;
}

cxxSolution::cxxSolution(struct solution *solution_ptr)
	//
	// constructor for cxxSolution from struct solution
	//
: cxxNumKeyword()
{
	description = solution_ptr->description;
	n_user      = solution_ptr->n_user;
	n_user_end  = solution_ptr->n_user_end;
	tc          = solution_ptr->tc;
	ph          = solution_ptr->ph;
	pe          = solution_ptr->solution_pe;
	mu          = solution_ptr->mu;
	ah2o        = solution_ptr->ah2o;
	total_h     = solution_ptr->total_h;
	total_o     = solution_ptr->total_o;
	cb          = solution_ptr->cb;
	mass_water  = solution_ptr->mass_water;
	total_alkalinity     = solution_ptr->total_alkalinity;
	// Totals, just save description and moles
	for (int i = 0; solution_ptr->totals[i].description != NULL; i++) {
		totals[solution_ptr->totals[i].description] = solution_ptr->totals[i].moles;
	}
	// Isotopes
	for (int i = 0; i < solution_ptr->count_isotopes; i++) {
		cxxIsotope iso = cxxIsotope(&solution_ptr->isotopes[i]);
		isotopes.push_back(iso);
	}
	// Master_activity
	for (int i = 0; i < solution_ptr->count_master_activity; i++) {
		master_activity[solution_ptr->master_activity[i].description] = solution_ptr->master_activity[i].la;
	}
	// Species_gammas
	for (int i = 0; i < solution_ptr->count_species_gamma; i++) {
		species_gamma[solution_ptr->species_gamma[i].description] = solution_ptr->species_gamma[i].la;
	}
}

cxxSolution::~cxxSolution()
{
}

struct solution *cxxSolution::cxxSolution2solution()
	//
	// Builds a solution structure from instance of cxxSolution 
	//
{
	struct solution *soln_ptr = solution_alloc();
	
	soln_ptr->description        = string_duplicate(this->description.c_str());
	soln_ptr->n_user             = this->n_user;
	soln_ptr->n_user_end         = this->n_user_end;
	soln_ptr->new_def            = 0;
	soln_ptr->tc                 = this->tc;
	soln_ptr->ph                 = this->ph;
	soln_ptr->solution_pe        = this->pe;
	soln_ptr->mu                 = this->mu;
	soln_ptr->ah2o               = this->ah2o;
	soln_ptr->total_h            = this->total_h;
	soln_ptr->total_o            = this->total_o;
	soln_ptr->cb                 = this->cb;
	soln_ptr->mass_water         = this->mass_water;
	soln_ptr->total_alkalinity   = this->total_alkalinity;
	soln_ptr->units              = moles_per_kilogram_string;
	soln_ptr->default_pe         = 0;
	// pe_data

	// totals
	soln_ptr->totals = (struct conc *) free_check_null(soln_ptr->totals);
	soln_ptr->totals = cxxConc::concarray((const std::map<char *, double>) this->totals);
	// master_activity
	{
		soln_ptr->master_activity = (struct master_activity *) PHRQ_realloc(soln_ptr->master_activity, (size_t) ((master_activity.size() + 1) * sizeof(struct master_activity)));
		int i = 0;
		if (soln_ptr->master_activity == NULL) malloc_error();
		for (std::map <char *, double>::iterator it = master_activity.begin(); it != master_activity.end(); ++it) {
			soln_ptr->master_activity[i].description = (char *)it->first;
			soln_ptr->master_activity[i].la = it->second;
			i++;
		}
		soln_ptr->master_activity[i].description = NULL;
		soln_ptr->count_master_activity = this->master_activity.size() + 1;
	}
	// species_gamma
	if (species_gamma.size() >= 0) {
		soln_ptr->species_gamma = (struct master_activity *) PHRQ_malloc((size_t) ((species_gamma.size() + 1) * sizeof(struct master_activity)));
		int i = 0;
		if (soln_ptr->species_gamma == NULL) malloc_error();
		for (std::map <char *, double>::iterator it = species_gamma.begin(); it != species_gamma.end(); ++it) {
			soln_ptr->species_gamma[i].description = (char *)it->first;
			soln_ptr->species_gamma[i].la = it->second;
			i++;
		}
		soln_ptr->species_gamma[i].description = NULL;
		soln_ptr->count_species_gamma = this->species_gamma.size() + 1;
	}
	// isotopes
	soln_ptr->isotopes = (struct isotope *) free_check_null(soln_ptr->isotopes);
	soln_ptr->isotopes = cxxIsotope::list2isotope(this->isotopes);

	return(soln_ptr);
}



#ifdef SKIP
cxxSolution& cxxSolution::read(CParser& parser)
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

	cxxSolution numkey;

	// Read solution number and description
	numkey.read_number_description(parser);

	// Malloc space for solution data
	//// g_solution_map[numkey.n_user()] = numkey;
	s_map[numkey.n_user()] = numkey;

	std::istream::pos_type ptr;
	std::istream::pos_type next_char;
	std::string token;
	CParser::TOKEN_TYPE j;
	
	cxxSolution& sol = s_map[numkey.n_user()];
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
			if (!(parser.get_iss() >> sol.tc))
			{
				sol.tc = 25;
			}
			break;

		case 2: // dens
		case 3: // density
			parser.get_iss() >> sol.density;
			break;

		case 4: // units
		case 8: // unit
			if (parser.copy_token(token, next_char) == CParser::TT_EMPTY) break;
			if (parser.check_units(token, false, false, sol.units, true) == CParser::OK) {
				sol.units = token;				
			} else {
				parser.incr_input_error();
			}
			break;

		case 5: // redox
			if (parser.copy_token(token, next_char) == CParser::TT_EMPTY) break;
			if (parser.parse_couple(token) == CParser::OK) {
				default_pe = cxxPe_Data::store(sol.pe, token);
			} else {
				parser.incr_input_error();
			}
			break;

		case 6: // ph
			{
				cxxConc conc;
				if (conc.read(parser, sol) == cxxConc::ERROR) {
					parser.incr_input_error();
					break;
				}
				sol.ph = conc.get_input_conc();
				if (conc.get_equation_name().empty()) {
					break;
				}
				conc.set_description("H(1)");
				sol.add(conc);
			}
			break;

		case 7: // pe
			{
				cxxConc conc;
				if (conc.read(parser, sol) == cxxConc::ERROR) {
					parser.incr_input_error();
					break;
				}
				sol.solution_pe = conc.get_input_conc();
				if (conc.get_equation_name().empty()) {
					break;
				}
				conc.set_description("E");
				sol.add(conc);
			}
			break;

		case 9: // isotope
			{
				cxxIsotope isotope;
				if (isotope.read(parser) == cxxIsotope::OK) {
					sol.add(isotope);
				}
			}
			break;

		case 10: // water
			j = parser.copy_token(token, next_char);
			if (j == CParser::TT_EMPTY) {
				sol.mass_water = 1.0;
			} else if (j != CParser::TT_DIGIT) {
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for mass of water in solution.", CParser::OT_CONTINUE);
			} else {
				std::istringstream(token) >> sol.mass_water;
			}
			break;

		case CParser::OPTION_DEFAULT:
			{
				//  Read concentration
				cxxConc conc;
				if (conc.read(parser, sol) == cxxConc::ERROR) {
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
	std::sort(sol.totals.begin(), sol.totals.end());
#endif

	//
	// fix up default units and default pe
	//
	std::string token1;
	std::vector<cxxConc>::iterator iter = sol.totals.begin();
	for (; iter != sol.totals.end(); ++iter)
	{
		token = (*iter).get_description();
		Utilities::str_tolower(token);
		if ((*iter).get_units().empty()) {
			(*iter).set_units(sol.units);
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
	sol.default_pe = default_pe;
	return sol;
}
#endif
void cxxSolution::dump_xml(std::ostream& os, unsigned int indent)const
{
	unsigned int i;

	for(i = 0; i < indent; ++i) os << Utilities::INDENT;		
	os << "<solution>\n";

	cxxNumKeyword::dump_xml(os, indent);

	for(i = 0; i < indent + 1; ++i) os << Utilities::INDENT;		
	os << "<temp>" << this->get_tc() << "</temp>" << "\n";

	for(i = 0; i < indent + 1; ++i) os << Utilities::INDENT;		
	os << "<pH>" << this->get_ph() << "</pH>" << "\n";

	for(i = 0; i < indent + 1; ++i) os << Utilities::INDENT;		
	os << "<pe>" << this->get_pe() << "</pe>" << "\n";

	//assert(this->pe.size() > 0);
	//assert(this->default_pe >= 0);
	//assert(this->pe.size() > (unsigned int) this->default_pe);
	//this->pe[this->default_pe].dump_xml(os, indent + 1);

	//for(i = 0; i < indent + 1; ++i) os << Utilities::INDENT;		
	//os << "<units>" << this->get_units() << "</units>" << "\n";

	//for(i = 0; i < indent + 1; ++i) os << Utilities::INDENT;		
	//os << "<density>" << this->get_density() << "</density>" << "\n";

	// foreach conc
	/*
	if (!this->totals.empty())
	{
		for(i = 0; i < indent + 1; ++i) os << Utilities::INDENT;		
		os << "<totals>\n";

		std::vector<cxxConc>::const_iterator iter = this->totals.begin();
		for(; iter != this->totals.end(); ++iter)
		{
			(*iter).dump_xml(*this, os, indent + 2);
		}

		for(i = 0; i < indent + 1; ++i) os << Utilities::INDENT;		
		os << "</totals>\n";
	}
	*/

	// foreach isotope
	if (!this->isotopes.empty())
	{
		for(i = 0; i < indent + 1; ++i) os << Utilities::INDENT;		
		os << "<isotopes>\n";

		std::list<cxxIsotope>::const_iterator iter = this->isotopes.begin();
		for(; iter != this->isotopes.end(); ++iter)
		{
			(*iter).dump_xml(os, indent + 2);
		}

		for(i = 0; i < indent + 1; ++i) os << Utilities::INDENT;		
		os << "</isotopes>\n";
	}

	for(i = 0; i < indent + 1; ++i) os << Utilities::INDENT;		
	os << "<water>" << this->get_mass_water() << "</water>" << "\n";

	for(i = 0; i < indent; ++i) os << Utilities::INDENT;		
	os << "</solution>" << "\n";
}
#include "ISolution.h"
void test_classes(void)
{
	int i;
	for (i=0; i < count_solution; i++) {
		if (solution[i]->new_def == true) {
			cxxISolution sol(solution[i]);
			solution[i] = (struct solution *) solution_free(solution[i]);
			solution[i] = sol.cxxISolution2solution();
		} else {
			cxxSolution sol(solution[i]);
			solution[i] = (struct solution *) solution_free(solution[i]);
			solution[i] = sol.cxxSolution2solution();
		}
	}
} 
