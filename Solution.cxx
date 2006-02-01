// Solution.cxx: implementation of the cxxSolution class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
#pragma warning(disable : 4786)   // disable truncation warning (Only used by debugger)
#endif

#include "Utils.h"   // define first
#include "Solution.h"
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
		//cxxIsotope iso = cxxIsotope(&solution_ptr->isotopes[i]);
		cxxIsotope iso(&solution_ptr->isotopes[i]);
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
	//const char    ERR_MESSAGE[] = "Packing solution message: %s, element not found\n";
	unsigned int i;
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
	//const char    ERR_MESSAGE[] = "Packing solution message: %s, element not found\n";
	unsigned int i;
	s_oss.precision(DBL_DIG - 1);
	std::string indent0(""), indent1(""), indent2("");
	for(i = 0; i < indent; ++i) indent0.append(Utilities::INDENT);
	for(i = 0; i < indent + 1; ++i) indent1.append(Utilities::INDENT);
	for(i = 0; i < indent + 2; ++i) indent2.append(Utilities::INDENT);

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
	s_oss << indent1;
	s_oss << "-totals" << std::endl;
	for (std::map <char *, double, CHARSTAR_LESS>::const_iterator it = totals.begin(); it != totals.end(); ++it) {
		s_oss << indent2;
		s_oss << it->first << "   " <<  it->second << std::endl;
	}

	// master_activity map
	s_oss << indent1;
	s_oss << "-activities" << std::endl;
	for (std::map <char *, double>::const_iterator it = master_activity.begin(); it != master_activity.end(); ++it) {
		s_oss << indent2;
		s_oss << it->first << "   " << it->second << std::endl;
	}

	// species_gamma map
	s_oss << indent1;
	s_oss << "-gammas" << std::endl;
	for (std::map <char *, double>::const_iterator it = species_gamma.begin(); it != species_gamma.end(); ++it) {
		s_oss << indent2;
		s_oss << it->first << "   " << it->second << std::endl;
	}

	// Isotopes
		s_oss << indent1;
	s_oss << "-Isotopes" << std::endl;
	for (std::list<cxxIsotope>::const_iterator it = this->isotopes.begin(); it != isotopes.end(); ++it) {
		it->dump_raw(s_oss, indent + 2);
	}

	return;
}

void cxxSolution::read_raw(CParser& parser)
{
	static std::vector<std::string> vopts;
	if (vopts.empty()) {
		vopts.reserve(21);
		vopts.push_back("totals");                    // 0 
		vopts.push_back("activities");		      // 1 
		vopts.push_back("gammas");    		      // 2 
		vopts.push_back("isotopes");  		      // 3 
		vopts.push_back("temp");             	      // 4 
		vopts.push_back("tc");               	      // 5 
		vopts.push_back("temperature");      	      // 6 
		vopts.push_back("ph");               	      // 7 
		vopts.push_back("pe");               	      // 8 
		vopts.push_back("mu");               	      // 9 
		vopts.push_back("ionic_strength");   	      // 10
		vopts.push_back("ah2o");             	      // 11
		vopts.push_back("activity_water");   	      // 12
		vopts.push_back("total_h");          	      // 13
		vopts.push_back("total_o");          	      // 14
		vopts.push_back("mass_water");       	      // 15
		vopts.push_back("mass_h2o");                  // 16
		vopts.push_back("total_alkalinity"); 	      // 17
		vopts.push_back("total_alk");        	      // 18
		vopts.push_back("cb");           	      // 19
		vopts.push_back("charge_balance");            // 20
	}

	cxxSolution numkey;

	// Read solution number and description
	numkey.read_number_description(parser);

	//skip keyword
	parser.get_line();

	std::istream::pos_type ptr;
	std::istream::pos_type next_char;
	std::string token;
	char * cstring;
	int opt_save;

	opt_save = CParser::OPT_ERROR;
	bool tc_defined(false); 
	bool ph_defined(false); 
	bool pe_defined(false); 
	bool mu_defined(false); 
	bool ah2o_defined(false); 
	bool total_h_defined(false); 
	bool total_o_defined(false); 
	bool cb_defined(false); 
	bool mass_water_defined(false); 
	bool total_alkalinity_defined(false);

	for (;;)
	{
		int opt = parser.get_option(vopts, next_char);
		if (opt == CParser::OPT_DEFAULT)
		{
			opt = opt_save;
		}

		switch (opt)
		{
		case CParser::OPT_EOF:
			break;
		case CParser::OPT_KEYWORD:
			break;
		case CParser::OPT_DEFAULT:
		case CParser::OPT_ERROR:
			opt = CParser::OPT_EOF;
			parser.error_msg("Unknown input in SOLUTION_RAW keyword.", CParser::OT_CONTINUE);
			parser.error_msg(parser.line().c_str(), CParser::OT_CONTINUE);
			break;

		case 0: // totals
			if (parser.copy_token(token, ptr) != CParser::TT_EMPTY) {
				if ( parser.addPair(this->totals, next_char) != CParser::PARSER_OK) {
					parser.incr_input_error();
					parser.error_msg("Expected element name and moles for totals.", CParser::OT_CONTINUE);
				}			
			}	
			opt_save = 0;
			break;

		case 1: // activities
			if ( parser.addPair(this->master_activity, next_char) != CParser::PARSER_OK) {
				parser.incr_input_error();
				parser.error_msg("Expected species name and log activity for activities.", CParser::OT_CONTINUE);
			}				
			opt_save = 1;
			break;

		case 2: // gammas
			if ( parser.addPair(this->species_gamma, next_char) != CParser::PARSER_OK) {
				parser.incr_input_error();
				parser.error_msg("Expected species name and activity coefficient for gammas.", CParser::OT_CONTINUE);
			}				
			opt_save = 2;
			break;

		case 3: // isotopes
			{
				cxxIsotope iso;
				iso.read_raw(parser);
				if ( iso.read_raw(parser) != CParser::PARSER_OK) {
					parser.incr_input_error();
					parser.error_msg("Expected data for isotopes.", CParser::OT_CONTINUE);
				} else {
					this->isotopes.push_back(iso);
				}
			}
			opt_save = 3;
			break;

		case 4: // temp
		case 5: // tc
		case 6: // temperature			
			if (!(parser.get_iss() >> this->tc))
			{
				this->tc = 25.0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for temperature.", CParser::OT_CONTINUE);
			}
			tc_defined = true;
			break;

		case 7: // ph
			if (!(parser.get_iss() >> this->ph))
			{
				this->ph = 7.0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for pH.", CParser::OT_CONTINUE);
			}
			ph_defined = true;
			break;

		case 8: // pe
			if (!(parser.get_iss() >> this->pe))
			{
				this->pe = 4.0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for pe.", CParser::OT_CONTINUE);
			}
			pe_defined = true;
			break;

		case 9: // mu
		case 10: // ionic_strength
			if (!(parser.get_iss() >> this->mu))
			{
				this->mu = 1e-7;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for ionic strength.", CParser::OT_CONTINUE);
			}
			mu_defined = true;
			break;

		case 11: // ah2o
		case 12: // activity_water
			if (!(parser.get_iss() >> this->ah2o))
			{
				this->ah2o = 1.0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for activity of water.", CParser::OT_CONTINUE);
			}
			ah2o_defined = true;
			break;

		case 13: // total_h
			if (!(parser.get_iss() >> this->total_h))
			{
				this->total_h = 111.1;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for total hydrogen.", CParser::OT_CONTINUE);
			}
			total_h_defined = true;
			break;

		case 14: // total_o
			if (!(parser.get_iss() >> this->total_o))
			{
				this->total_o = 55.55;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for total oxygen.", CParser::OT_CONTINUE);
			}
			total_o_defined = true;
			break;

		case 15: // mass_water
		case 16: // mass_h2o
			if (!(parser.get_iss() >> this->mass_water))
			{
				this->mass_water = 1.0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for mass of water.", CParser::OT_CONTINUE);
			}
			mass_water_defined = true;
			break;

		case 17: // total_alkalinity
		case 18: // total_alk
			if (!(parser.get_iss() >> this->total_alkalinity))
			{
				this->total_alkalinity = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for total_alkalinity.", CParser::OT_CONTINUE);
			}
			total_alkalinity_defined = true;
			break;

		case 19: // cb
		case 20: // charge_balance
			if (!(parser.get_iss() >> this->cb))
			{
				this->cb = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for charge balance.", CParser::OT_CONTINUE);
			}
			cb_defined = true;
			break;

		}
		if (opt == CParser::OPT_EOF || opt == CParser::OPT_KEYWORD) break;
	}
#ifdef SKIP
	//
	// Sort totals by description
	//
	std::sort(sol.totals.begin(), sol.totals.end());
#endif

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
#include <iostream>     // std::cout std::cerr
#include <strstream>
#include <fstream>
void test_classes(void)
{
	int i;
	for (i=0; i < count_solution; i++) {
		if (solution[i]->new_def == TRUE) {
			cxxISolution sol(solution[i]);
			solution[i] = (struct solution *) solution_free(solution[i]);
			solution[i] = sol.cxxISolution2solution();
		} else {
			std::strstream oss;
			cxxSolution sol(solution[i]);
			solution[i] = (struct solution *) solution_free(solution[i]);
			solution[i] = sol.cxxSolution2solution();
			sol.dump_raw(oss, 0);
			cxxSolution sol1;

			//std::ostringstream oss_out;
			//std::ostringstream oss_err;
			//std::istringstream iss_in();

			std::fstream myfile("t");


			CParser cparser(myfile, std::cout, std::cerr);

			sol1.read_raw(cparser);
			
		}
	}
} 
