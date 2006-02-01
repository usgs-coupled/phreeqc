// Solution.cxx: implementation of the cxxSolution class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
#pragma warning(disable : 4786)   // disable truncation warning (Only used by debugger)
#endif

#include "Solution.h"
#include "Utils.h"   // define before global.h
#define EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "phrqproto.h"
#include <cassert>     // assert
#include <algorithm>   // std::sort 



//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

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
	int i;

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
	for (i = 0; solution_ptr->totals[i].description != NULL; i++) {
		if (solution_ptr->totals[i].description == NULL) continue;
		totals[solution_ptr->totals[i].description] = solution_ptr->totals[i].moles;
	}

	// Isotopes
	for (i = 0; i < solution_ptr->count_isotopes; i++) {
		cxxIsotope iso = cxxIsotope(&solution_ptr->isotopes[i]);
		isotopes.push_back(iso);
	}

	// Master_activity
	for (i = 0; i < solution_ptr->count_master_activity ; i++) {
		if (solution_ptr->master_activity[i].description == NULL) continue;
		master_activity[solution_ptr->master_activity[i].description] = solution_ptr->master_activity[i].la;
	}

	// Species_gammas
	for (i = 0; i < solution_ptr->count_species_gamma; i++) {
		if (solution_ptr->species_gamma[i].description == NULL) continue;
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
	int i;

	struct solution *soln_ptr = solution_alloc();
	
	soln_ptr->description        = this->get_description();
	soln_ptr->n_user             = this->n_user;
	soln_ptr->n_user_end         = this->n_user_end;
	soln_ptr->new_def            = FALSE;
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
	soln_ptr->density            = 1.0;
	soln_ptr->units              = moles_per_kilogram_string;
	soln_ptr->default_pe         = 0;
	// pe_data

	// totals
	soln_ptr->totals = (struct conc *) free_check_null(soln_ptr->totals);
	//soln_ptr->totals = cxxConc::concarray((const std::map<char *, double>) this->totals);
	soln_ptr->totals = cxxConc::concarray(this->totals);

	// master_activity
	soln_ptr->master_activity = (struct master_activity *) PHRQ_realloc(soln_ptr->master_activity, (size_t) ((master_activity.size() + 1) * sizeof(struct master_activity)));
	if (soln_ptr->master_activity == NULL) malloc_error();
	i = 0;
	for (std::map <char *, double>::iterator it = master_activity.begin(); it != master_activity.end(); it++) {
		soln_ptr->master_activity[i].description = (char *)it->first;
		soln_ptr->master_activity[i].la = it->second;
		i++;
	}
	soln_ptr->master_activity[i].description = NULL;
	soln_ptr->count_master_activity = this->master_activity.size() + 1;

	// species_gamma
	if (species_gamma.size() >= 0) {
		soln_ptr->species_gamma = (struct master_activity *) PHRQ_malloc((size_t) ((species_gamma.size()) * sizeof(struct master_activity)));
		int i = 0;
		if (soln_ptr->species_gamma == NULL) malloc_error();
		for (std::map <char *, double>::iterator it = species_gamma.begin(); it != species_gamma.end(); ++it) {
			soln_ptr->species_gamma[i].description = (char *)it->first;
			soln_ptr->species_gamma[i].la = it->second;
			i++;
		}
		soln_ptr->count_species_gamma = this->species_gamma.size();
	} else {
		soln_ptr->species_gamma = NULL;
		soln_ptr->count_species_gamma = 0;
	}
	// isotopes
	soln_ptr->isotopes = (struct isotope *) free_check_null(soln_ptr->isotopes);
	soln_ptr->isotopes = cxxIsotope::list2isotope(this->isotopes);
	soln_ptr->count_isotopes = this->isotopes.size();

	return(soln_ptr);
}
void cxxSolution::dump_xml(std::ostream& s_oss, unsigned int indent)const
{
	const char    ERR_MESSAGE[] = "Packing solution message: %s, element not found\n";
	int i;
	s_oss.precision(DBL_DIG - 1);
	std::string indent0(""), indent1("");
	for(i = 0; i < indent; ++i) indent0.append(Utilities::INDENT);
	for(i = 0; i < indent + 1; ++i) indent1.append(Utilities::INDENT);

	// Solution element and attributes
	s_oss << indent0;
	s_oss << "<solution " << std::endl;

	//s_oss << indent1;
	//s_oss << "soln_new_def=\"" << this->new_def << "\"" << std::endl;

	s_oss << indent1;
	s_oss << "soln_n_user=\"" << this->n_user << "\" " << std::endl;

	s_oss << indent1;
	s_oss << "soln_description=\"" << this->description << "\"" << std::endl;

	s_oss << indent1;
	s_oss << "soln_tc=\"" << this->tc << "\"" << std::endl;

	s_oss << indent1;
	s_oss << "soln_ph=\"" << this->ph << "\"" << std::endl;

	s_oss << indent1;
	s_oss << "soln_solution_pe=\"" << this->pe << "\"" << std::endl;

	s_oss << indent1;
	s_oss << "soln_mu=\"" << this->mu << "\"" << std::endl;

	s_oss << indent1;
	s_oss << "soln_ah2o=\""	 << this->ah2o << "\"" << std::endl;

	s_oss << indent1;
	s_oss << "soln_total_h=\"" << this->total_h << "\"" << std::endl;

	s_oss << indent1;
	s_oss << "soln_total_o=\"" << this->total_o << "\"" << std::endl;

	s_oss << indent1;
	s_oss << "soln_cb=\"" << this->cb << "\"" << std::endl;

	s_oss << indent1;
	s_oss << "soln_mass_water=\"" << this->mass_water << "\"" << std::endl;

	s_oss << indent1;
	s_oss << "soln_total_alkalinity=\"" << this->total_alkalinity << "\"" << std::endl;

	s_oss << indent1;
	s_oss << "\">" << std::endl;

	// soln_total conc structures
	for (std::map <char *, double, CHARSTAR_LESS>::const_iterator it = totals.begin(); it != totals.end(); ++it) {
		s_oss << indent1;
		s_oss << "<soln_total";
		s_oss << " conc_desc=\"" << it->first << "\"";
		s_oss << " conc_moles=\"" << it->second << "\"" ;
		s_oss << "\">" << std::endl;
	}

	// master_activity map
	for (std::map <char *, double>::const_iterator it = master_activity.begin(); it != master_activity.end(); ++it) {
		s_oss << indent1;
		s_oss << "<soln_m_a";
		s_oss << " m_a_desc=\"" << it->first << "\"" ;
		s_oss << " m_a_la=\"" << it->second << "\"" ;
		s_oss << "\">" << std::endl;
	}

	// species_gamma map
	for (std::map <char *, double>::const_iterator it = species_gamma.begin(); it != species_gamma.end(); ++it) {
		s_oss << indent1;
		s_oss << "<soln_s_g";
		s_oss << " m_a_desc=\"" << it->first << "\"" ;
		s_oss << " m_a_la=\"" << it->second << "\"" ;
		s_oss << "\">" << std::endl;
	}

	for (std::list<cxxIsotope>::const_iterator it = this->isotopes.begin(); it != isotopes.end(); ++it) {
		it->dump_xml(s_oss, indent + 1);
	}

	// End of solution
	s_oss << indent0;
	s_oss << "</solution>" << std::endl;

	return;
}

void cxxSolution::dump_raw(std::ostream& s_oss, unsigned int indent)const
{
	const char    ERR_MESSAGE[] = "Packing solution message: %s, element not found\n";
	int i;
	s_oss.precision(DBL_DIG - 1);
	std::string indent0(""), indent1("");
	for(i = 0; i < indent; ++i) indent0.append(Utilities::INDENT);
	for(i = 0; i < indent + 1; ++i) indent1.append(Utilities::INDENT);

	// Solution element and attributes
	s_oss << indent0;
	s_oss << "SOLUTION_RAW       " << this->n_user  << " " << this->description << std::endl;

	s_oss << indent1;
	s_oss << "-temp              " << this->tc << std::endl;

	s_oss << indent1;
	s_oss << "-pH                " << this->ph << std::endl;

	s_oss << indent1;
	s_oss << "-pe                " << this->pe << std::endl;

	// new identifier
	s_oss << indent1;
	s_oss << "-mu                " << this->mu << std::endl;

	// new identifier
	s_oss << indent1;
	s_oss << "-ah2o              " << this->ah2o << std::endl;

	// new identifier
	s_oss << indent1;
	s_oss << "-total_h           " << this->total_h << std::endl;

	// new identifier
	s_oss << indent1;
	s_oss << "-total_o           " << this->total_o << std::endl;

	// new identifier
	s_oss << indent1;
	s_oss << "-cb                " << this->cb << std::endl;

	// new identifier
	s_oss << indent1;
	s_oss << "-mass_water        " << this->mass_water << std::endl;

	// new identifier
	s_oss << indent1;
	s_oss << "-total_alkalinity  " << this->total_alkalinity << std::endl;

	// soln_total conc structures
	for (std::map <char *, double, CHARSTAR_LESS>::const_iterator it = totals.begin(); it != totals.end(); ++it) {
		s_oss << indent1;
		s_oss << "-tot   "  << it->first << "   " <<  it->second << std::endl;
	}

	// master_activity map
	for (std::map <char *, double>::const_iterator it = master_activity.begin(); it != master_activity.end(); ++it) {
		s_oss << indent1;
		s_oss << "-act  "  << it->first << "   " << it->second << std::endl;
	}

	// species_gamma map
	for (std::map <char *, double>::const_iterator it = species_gamma.begin(); it != species_gamma.end(); ++it) {
		s_oss << indent1;
		s_oss << "-gam  "  << it->first << "   " << it->second << std::endl;
	}

	for (std::list<cxxIsotope>::const_iterator it = this->isotopes.begin(); it != isotopes.end(); ++it) {
		it->dump_raw(s_oss, indent + 1);
	}

	// End of solution
	//s_oss << indent0;
	//s_oss << "SOLUTION_RAW_END" << std::endl;

	return;
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

	std::istream::pos_type ptr;
	std::istream::pos_type next_char;
	std::string token;
	CParser::TOKEN_TYPE j;
	
	//cxxSolution& sol = s_map[numkey.n_user()];
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


#include "ISolution.h"
void test_classes(void)
{
	int i;
	for (i=0; i < count_solution; i++) {
		if (solution[i]->new_def == TRUE) {
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
