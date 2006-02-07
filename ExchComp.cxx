// Solution.cxx: implementation of the cxxExchComp class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
#pragma warning(disable : 4786)   // disable truncation warning (Only used by debugger)
#endif

#include "Utils.h"   // define first
#include "ExchComp.h"
#define EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "phrqproto.h"
#include <cassert>     // assert
#include <algorithm>   // std::sort 

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cxxExchComp::cxxExchComp()
        //
        // default constructor for cxxExchComp 
        //
{
        moles                   = 0.0;
        la                      = 0.0;
        charge_balance          = 0.0;
        phase_name              = NULL;
        phase_proportion        = 0.0;
	rate_name               = NULL;
}

cxxExchComp::cxxExchComp(struct exch_comp *exch_comp_ptr)
        //
        // constructor for cxxExchComp from struct exch_comp
        //
: 
formula_totals(exch_comp_ptr->formula_totals),
totals(exch_comp_ptr->totals)
{
	formula                  = exch_comp_ptr->formula;
	formula_z                = exch_comp_ptr->formula_z;
	moles                    = exch_comp_ptr->moles;
        // totals in constructor
	//formula_totals in constructor
	la                       = exch_comp_ptr->la;
	charge_balance           = exch_comp_ptr->charge_balance;
	phase_name               = exch_comp_ptr->phase_name;
	phase_proportion         = exch_comp_ptr->phase_proportion;
	rate_name                = exch_comp_ptr->rate_name;
}

cxxExchComp::~cxxExchComp()
{
}

struct master *cxxExchComp::get_master()
{	
	struct master *master_ptr = NULL;
	for (std::map <char *, double, CHARSTAR_LESS>::iterator it = totals.begin(); it != totals.end(); it++) {

		/* Find master species */
		char *eltName = it->first;
		struct element *elt_ptr = element_store(eltName);
		if (elt_ptr->master == NULL) {
                        std::ostringstream error_oss;
			error_oss << "Master species not in data base for " << elt_ptr->name << std::endl;
			Utilities::error_msg(error_oss.str(), CONTINUE);
		}
		if (elt_ptr->master->type != EX) continue;
		master_ptr = elt_ptr->master;
		break;
	}
	if (master_ptr == NULL) {
		std::ostringstream error_oss;
		error_oss << "Exchange formula does not contain an exchange master species, " << this->formula << std::endl;
		Utilities::error_msg(error_oss.str(), CONTINUE);
	}
	return(master_ptr);
}

struct exch_comp *cxxExchComp::cxxExchComp2exch_comp(std::list<cxxExchComp>& el)
        //
        // Builds exch_comp structure from of cxxExchComp 
        //
{
	struct exch_comp *exch_comp_ptr = (struct exch_comp *) PHRQ_malloc((size_t) (el.size() * sizeof(struct exch_comp)));
	if (exch_comp_ptr == NULL) malloc_error();

	int i = 0;
	for (std::list<cxxExchComp>::iterator it = el.begin(); it != el.end(); ++it) {
		exch_comp_ptr[i].formula		=  it->formula;
		exch_comp_ptr[i].formula_z		=  it->formula_z;
		exch_comp_ptr[i].formula_totals         =  it->formula_totals.elt_list();
		exch_comp_ptr[i].moles			=  it->moles;
		exch_comp_ptr[i].totals                 =  it->formula_totals.elt_list();
		exch_comp_ptr[i].la			=  it->la;
		exch_comp_ptr[i].charge_balance		=  it->charge_balance;
		exch_comp_ptr[i].phase_name		=  it->phase_name;
		exch_comp_ptr[i].phase_proportion	=  it->phase_proportion;
		exch_comp_ptr[i].rate_name            	=  it->rate_name;
		i++;
	}
        return(exch_comp_ptr);
}

#ifdef SKIP
struct exch_comp *cxxExchComp::cxxExchComp2exch_comp()
        //
        // Builds exch_comp structure from of cxxExchComp 
        //
{
	struct exch_comp *exch_comp_ptr = (struct exch_comp *) PHRQ_malloc((size_t) (this->totals.size() * sizeof(struct exch_comp)));
	if (exch_comp_ptr == NULL) malloc_error();

	int i = 0;
	for (std::list::iterator it = this->totals.begin(); it != totals.end(); ++it) {
		exch_comp_ptr->formula		        =  formula;
		exch_comp_ptr->formula_z		=  formula_z;
		exch_comp_ptr->formula_totals           =  formula_totals.elt_list();
		exch_comp_ptr->moles			=  moles;
		exch_comp_ptr->totals                   =  formula_totals.elt_list();
		exch_comp_ptr->la			=  la;
		exch_comp_ptr->charge_balance		=  charge_balance;
		exch_comp_ptr->phase_name		=  phase_name;
		exch_comp_ptr->phase_proportion		=  phase_proportion;
		exch_comp_ptr->rate_name            	=  rate_name;
		i++;
	}
        return(exch_comp_ptr);
}
#endif

void cxxExchComp::dump_xml(std::ostream& s_oss, unsigned int indent)const
{
        //const char    ERR_MESSAGE[] = "Packing exch_comp message: %s, element not found\n";
        unsigned int i;
	s_oss.precision(DBL_DIG - 1);
        std::string indent0(""), indent1(""), indent2("");
        for(i = 0; i < indent; ++i) indent0.append(Utilities::INDENT);
        for(i = 0; i < indent + 1; ++i) indent1.append(Utilities::INDENT);
        for(i = 0; i < indent + 2; ++i) indent2.append(Utilities::INDENT);

        // Exch_Comp element and attributes

        s_oss << indent0 << "formula=\"" << this->formula << "\"" << std::endl;
        s_oss << indent0 << "exchange_name=\"" << this->exchange_name << "\"" << std::endl;
        s_oss << indent0 << "moles=\"" << this->moles  << "\"" << std::endl;
        s_oss << indent0 << "la=\"" << this->la     << "\"" << std::endl;
        s_oss << indent0 << "charge_balance=\"" << this->charge_balance << "\"" << std::endl;
	if (this->phase_name != NULL) {
		s_oss << indent0 << "phase_name=\"" << this->phase_name << "\"" << std::endl;
	}
	if (this->rate_name != NULL) {
		s_oss << indent0 << "rate_name=\"" << this->rate_name << "\"" << std::endl;
	}
        s_oss << indent0 << "phase_proportion=\"" << this->phase_proportion  << "\"" << std::endl;

        // totals
        s_oss << indent0;
        s_oss << "<totals " << std::endl;
	this->totals.dump_xml(s_oss, indent + 1);

        // formula_totals
        s_oss << indent0;
        s_oss << "<formula_totals " << std::endl;
	this->formula_totals.dump_xml(s_oss, indent + 1);
}

void cxxExchComp::dump_raw(std::ostream& s_oss, unsigned int indent)const
{
        //const char    ERR_MESSAGE[] = "Packing exch_comp message: %s, element not found\n";
        unsigned int i;
	s_oss.precision(DBL_DIG - 1);
        std::string indent0(""), indent1(""), indent2("");
        for(i = 0; i < indent; ++i) indent0.append(Utilities::INDENT);
        for(i = 0; i < indent + 1; ++i) indent1.append(Utilities::INDENT);
        for(i = 0; i < indent + 2; ++i) indent2.append(Utilities::INDENT);

        // Exch_Comp element and attributes

        s_oss << indent0 << "-formula               " << this->formula << std::endl;
        s_oss << indent0 << "-exchange_name         " << this->exchange_name << std::endl;
        s_oss << indent0 << "-moles                 " << this->moles  << std::endl;
        s_oss << indent0 << "-la                    " << this->la     << std::endl;
        s_oss << indent0 << "-charge_balance        " << this->charge_balance << std::endl;
	if (this->phase_name != NULL) {
		s_oss << indent0 << "-phase_name            " << this->phase_name << std::endl;
	}
	if (this->rate_name != NULL) {
		s_oss << indent0 << "-rate_name             " << this->rate_name << std::endl;
	}
        s_oss << indent0 << "-phase_proportion              " << this->phase_proportion  << std::endl;

        // totals
        s_oss << indent1;
        s_oss << "-totals" << std::endl;
	this->totals.dump_raw(s_oss, indent + 2);
	/*
        for (std::map <char *, double, CHARSTAR_LESS>::const_iterator it = totals.begin(); it != totals.end(); ++it) {
                s_oss << indent2;
                s_oss << it->first << "   " <<  it->second << std::endl;
        }
	*/
        // formula_totals
        s_oss << indent1;
        s_oss << "-formula_totals" << std::endl;
	this->formula_totals.dump_raw(s_oss, indent + 2);
}

#ifdef SKIP
void cxxExchComp::read_raw(CParser& parser)
{
        static std::vector<std::string> vopts;
        if (vopts.empty()) {
                vopts.reserve(10);
                vopts.push_back("formula");                   // 0 
                vopts.push_back("exchange_name");             // 1 
                vopts.push_back("moles");                     // 2 
                vopts.push_back("la");                        // 3 
                vopts.push_back("charge_balance");            // 4 
                vopts.push_back("phase_name");                // 5 
                vopts.push_back("rate_name");                 // 6 
                vopts.push_back("phase_proportion");          // 7 
                vopts.push_back("totals");                    // 8
        }

        std::istream::pos_type ptr;
        std::istream::pos_type next_char;
        std::string token;
        int opt_save;

        opt_save = CParser::OPT_ERROR;
        bool formula_defined(false); 
        bool exchange_name_defined(false); 
        bool moles_defined(false); 
        bool la_defined(false); 
        bool charge_balance_defined(false); 

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
                        parser.error_msg("Unknown input in EXCH_COMP_RAW keyword.", CParser::OT_CONTINUE);
                        parser.error_msg(parser.line().c_str(), CParser::OT_CONTINUE);
                        break;

                case 0: // totals
                        if ( parser.addPair(this->totals, next_char) != CParser::PARSER_OK) {
                                parser.incr_input_error();
                                parser.error_msg("Expected element name and moles for totals.", CParser::OT_CONTINUE);
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
				if ( iso.read_raw(parser) != CParser::PARSER_OK) {
                                        parser.incr_input_error();
                                        parser.error_msg("Expected data for isotopes.", CParser::OT_CONTINUE);
                                } else {
					if (iso.get_isotope_name() != NULL) {
						this->isotopes.push_back(iso);
					}
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
	// all members must be defined
        if (tc_defined == false) {
		parser.incr_input_error();
		parser.error_msg("Temp not defined for EXCH_COMP_RAW input.", CParser::OT_CONTINUE);
	}
	if (ph_defined == false) {
		parser.incr_input_error();
		parser.error_msg("pH not defined for EXCH_COMP_RAW input.", CParser::OT_CONTINUE);
	}
	if (pe_defined == false) {
		parser.incr_input_error();
		parser.error_msg("pe not defined for EXCH_COMP_RAW input.", CParser::OT_CONTINUE);
	}
	if (mu_defined == false) {
		parser.incr_input_error();
		parser.error_msg("Ionic strength not defined for EXCH_COMP_RAW input.", CParser::OT_CONTINUE);
	}
	if (ah2o_defined == false) {
		parser.incr_input_error();
		parser.error_msg("Activity of water not defined for EXCH_COMP_RAW input.", CParser::OT_CONTINUE);
	}
	if (total_h_defined == false) {
		parser.incr_input_error();
		parser.error_msg("Total hydrogen not defined for EXCH_COMP_RAW input.", CParser::OT_CONTINUE);
	}
	if (total_o_defined == false) {
		parser.incr_input_error();
		parser.error_msg("Total oxygen not defined for EXCH_COMP_RAW input.", CParser::OT_CONTINUE);
	}
	if (cb_defined == false) {
		parser.incr_input_error();
		parser.error_msg("Charge balance not defined for EXCH_COMP_RAW input.", CParser::OT_CONTINUE);
	}
	if (mass_water_defined == false) {
		parser.incr_input_error();
		parser.error_msg("Temp not defined for EXCH_COMP_RAW input.", CParser::OT_CONTINUE);
	}
	if (total_alkalinity_defined == false) {
		parser.incr_input_error();
		parser.error_msg("Total alkalinity not defined for EXCH_COMP_RAW input.", CParser::OT_CONTINUE);
	}
        return;
}
#endif
#ifdef SKIP
cxxExchComp& cxxExchComp::read(CParser& parser)
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

        cxxExchComp numkey;

        // Read exch_comp number and description
        numkey.read_number_description(parser);

        std::istream::pos_type ptr;
        std::istream::pos_type next_char;
        std::string token;
        CParser::TOKEN_TYPE j;
        
        //cxxExchComp& sol = s_map[numkey.n_user()];
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
                        parser.error_msg("Unknown input in EXCH_COMP keyword.", CParser::OT_CONTINUE);
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
                                sol.exch_comp_pe = conc.get_input_conc();
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
                                parser.error_msg("Expected numeric value for mass of water in exch_comp.", CParser::OT_CONTINUE);
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

