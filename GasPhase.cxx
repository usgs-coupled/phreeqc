// GasPhase.cxx: implementation of the cxxGasPhase class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
#pragma warning(disable : 4786)   // disable truncation warning (Only used by debugger)
#endif

#include "Utils.h"   // define first
#include "GasPhase.h"
#define EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "phrqproto.h"
#include <cassert>     // assert
#include <algorithm>   // std::sort 

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cxxGasPhase::cxxGasPhase()
        //
        // default constructor for cxxGasPhase 
        //
: cxxNumKeyword()
{
        total_p                     = 0;
	volume                      = 0;
        gasPhaseComps.type          = cxxNameDouble::ND_NAME_COEF;
}

cxxGasPhase::cxxGasPhase(struct gas_phase *gas_phase_ptr)
        //
        // constructor for cxxGasPhase from struct gas_phase
        //
: 
cxxNumKeyword()
{
        int i;

        this->set_description(gas_phase_ptr->description);
        n_user                       = gas_phase_ptr->n_user;     
        n_user_end                   = gas_phase_ptr->n_user_end;  
	if (gas_phase_ptr->type == PRESSURE) {
		type                 = cxxGasPhase::GP_PRESSURE;
	} else {
		type                 = cxxGasPhase::GP_VOLUME;
	}	
        total_p                      = gas_phase_ptr->total_p; 
        volume                       = gas_phase_ptr->volume; 

        // gas_phase components
        for (i = 0; i < gas_phase_ptr->count_comps; i++) {
		gasPhaseComps[gas_phase_ptr->comps[i].name] = gas_phase_ptr->comps[i].moles;
        }
}

cxxGasPhase::~cxxGasPhase()
{
}

struct gas_comp * cxxGasPhase::cxxGasPhaseComp2gas_comp()
{
	struct gas_comp *gas_comp_ptr(NULL);
	if (this->gasPhaseComps.size() > 0) {
		int i = 0;
		int n;
		gas_comp_ptr = (struct gas_comp *) PHRQ_malloc((size_t) (this->gasPhaseComps.size() * sizeof(struct gas_comp)));
		if (gas_comp_ptr == NULL) malloc_error();
		for (cxxNameDouble::iterator it = this->gasPhaseComps.begin(); it != this->gasPhaseComps.end(); it++) {
			gas_comp_ptr[i].name = it->first;
			gas_comp_ptr[i].phase = phase_bsearch(it->first, &n, TRUE);
			gas_comp_ptr[i].p_read = 0;
			gas_comp_ptr[i].moles = it->second;
			gas_comp_ptr[i].initial_moles = 0;
			i++;
		}
	}
	return(gas_comp_ptr);
}

struct gas_phase *cxxGasPhase::cxxGasPhase2gas_phase()
        //
        // Builds a gas_phase structure from instance of cxxGasPhase 
        //
{
        struct gas_phase *gas_phase_ptr = gas_phase_alloc();
        
        gas_phase_ptr->description                 = this->get_description();  
        gas_phase_ptr->n_user                      = this->n_user;            
        gas_phase_ptr->n_user_end                  = this->n_user_end;        
        gas_phase_ptr->new_def                     = FALSE;
        gas_phase_ptr->solution_equilibria         = FALSE;
        gas_phase_ptr->n_solution                  = -2;
	if (this->type == cxxGasPhase::GP_PRESSURE) {
		gas_phase_ptr->type                = PRESSURE;
	} else {
		gas_phase_ptr->type                = VOLUME;
	}
        gas_phase_ptr->total_p                     = this->total_p;
        gas_phase_ptr->volume                      = this->volume;
        gas_phase_ptr->temperature                 = 273.15;

        // comps
        gas_phase_ptr->count_comps                 = this->gasPhaseComps.size();
        gas_phase_ptr->comps                       = (struct gas_comp *) free_check_null(gas_phase_ptr->comps);
	gas_phase_ptr->comps                       = this->cxxGasPhaseComp2gas_comp();

        return(gas_phase_ptr);
}

#ifdef SKIP
void cxxGasPhase::dump_xml(std::ostream& s_oss, unsigned int indent)const
{
        //const char    ERR_MESSAGE[] = "Packing gas_phase message: %s, element not found\n";
        unsigned int i;
        s_oss.precision(DBL_DIG - 1);
        std::string indent0(""), indent1(""), indent2("");
        for(i = 0; i < indent; ++i) indent0.append(Utilities::INDENT);
        for(i = 0; i < indent + 1; ++i) indent1.append(Utilities::INDENT);
        for(i = 0; i < indent + 2; ++i) indent2.append(Utilities::INDENT);

        // GasPhase element and attributes
        s_oss << indent0;
        s_oss << "<gas_phase " << std::endl;

        s_oss << indent1;
        s_oss << "pitzer_gas_phase_gammas=\"" << this->pitzer_gas_phase_gammas << "\"" << std::endl;

        // components
        s_oss << indent1;
        s_oss << "<component " << std::endl;
        for (std::list<cxxGasPhaseComp>::const_iterator it = gas_phaseComps.begin(); it != gas_phaseComps.end(); ++it) {
                it->dump_xml(s_oss, indent + 2);
        }

        return;
}
#endif

void cxxGasPhase::dump_raw(std::ostream& s_oss, unsigned int indent)const
{
        //const char    ERR_MESSAGE[] = "Packing gas_phase message: %s, element not found\n";
        unsigned int i;
        s_oss.precision(DBL_DIG - 1);
        std::string indent0(""), indent1(""), indent2("");
        for(i = 0; i < indent; ++i) indent0.append(Utilities::INDENT);
        for(i = 0; i < indent + 1; ++i) indent1.append(Utilities::INDENT);
        for(i = 0; i < indent + 2; ++i) indent2.append(Utilities::INDENT);

        // GasPhase element and attributes
        s_oss << indent0;
        s_oss << "GAS_PHASE_RAW       " << this->n_user  << " " << this->description << std::endl;

        s_oss << indent1;
        s_oss << "-type               " << this->type  << std::endl;

        s_oss << indent1;
        s_oss << "-total_p            " << this->total_p  << std::endl;

        s_oss << indent1;
        s_oss << "-volume             " << this->volume  << std::endl;

        // gasPhaseComps 
	s_oss << indent1;
	s_oss << "-component" << std::endl;
	this->gasPhaseComps.dump_raw(s_oss, indent + 2);
}

void cxxGasPhase::read_raw(CParser& parser)
{

	int i;
        static std::vector<std::string> vopts;
        if (vopts.empty()) {
                vopts.reserve(15);
                vopts.push_back("type");                  //0
                vopts.push_back("total_p");               //1
                vopts.push_back("volume");                //2
                vopts.push_back("component");             //3
        }                                                      
                                                               
        std::istream::pos_type ptr;                            
        std::istream::pos_type next_char;                      
        std::string token;                                     
        int opt_save;
        bool useLastLine(false);

        // Read gas_phase number and description
        this->read_number_description(parser);

        opt_save = CParser::OPT_ERROR;
        bool type_defined(false); 
        bool total_p_defined(false); 
        bool volume_defined(false); 

        for (;;)
        {
                int opt;
                if (useLastLine == false) {
                        opt = parser.get_option(vopts, next_char);
                } else {
                        opt = parser.getOptionFromLastLine(vopts, next_char);
                }
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
                        parser.error_msg("Unknown input in GAS_PHASE_COMP_RAW keyword.", CParser::OT_CONTINUE);
                        parser.error_msg(parser.line().c_str(), CParser::OT_CONTINUE);
                        useLastLine = false;
                        break;

                case 0: // type
			if (!(parser.get_iss() >> i))
                        //if (!(parser.get_iss() >> (cxxGasPhase::GP_TYPE) this->type))
                        {
                                this->type = cxxGasPhase::GP_PRESSURE;
                                parser.incr_input_error();
                                parser.error_msg("Expected enum for type.", CParser::OT_CONTINUE);
                        } else {
				this->type = (cxxGasPhase::GP_TYPE) i;
			}
                        type_defined = true;
                        useLastLine = false;
                        break;

                case 1: // total_p
                        if (!(parser.get_iss() >> this->total_p))
                        {
                                this->total_p = 0;
                                parser.incr_input_error();
                                parser.error_msg("Expected numeric value for total_p.", CParser::OT_CONTINUE);
                        }
                        total_p_defined = true;
                        useLastLine = false;
                        break;

                case 2: // volume
                        if (!(parser.get_iss() >> this->volume))
                        {
                                this->volume = 0;
                                parser.incr_input_error();
                                parser.error_msg("Expected numeric value for volume.", CParser::OT_CONTINUE);
                        }
                        volume_defined = true;
                        useLastLine = false;
                        break;

                case 3: // component
                        if ( this->gasPhaseComps.read_raw(parser, next_char) != CParser::PARSER_OK) {
                                parser.incr_input_error();
                                parser.error_msg("Expected gas component name and moles for gasPhaseComps.", CParser::OT_CONTINUE);
                        }                               
                        opt_save = 3;
                        useLastLine = false;
                        break;
                }
                if (opt == CParser::OPT_EOF || opt == CParser::OPT_KEYWORD) break;
        }
        // members that must be defined
        if (type_defined == false) {
                parser.incr_input_error();
                parser.error_msg("Type not defined for GAS_PHASE_RAW input.", CParser::OT_CONTINUE);
        }
        if (total_p_defined == false) {
                parser.incr_input_error();
                parser.error_msg("Total_p not defined for GAS_PHASE_RAW input.", CParser::OT_CONTINUE);
        }
        if (volume_defined == false) {
                parser.incr_input_error();
                parser.error_msg("Volume not defined for GAS_PHASE_RAW input.", CParser::OT_CONTINUE);
        }
}
#ifdef SKIP
cxxGasPhase& cxxGasPhase::read(CParser& parser)
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

        cxxGasPhase numkey;

        // Read gas_phase number and description
        numkey.read_number_description(parser);

        std::istream::pos_type ptr;
        std::istream::pos_type next_char;
        std::string token;
        CParser::TOKEN_TYPE j;
        
        //cxxGasPhase& sol = s_map[numkey.n_user()];
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
                        parser.error_msg("Unknown input in GAS_PHASE keyword.", CParser::OT_CONTINUE);
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
                                sol.gas_phase_pe = conc.get_input_conc();
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
                                parser.error_msg("Expected numeric value for mass of water in gas_phase.", CParser::OT_CONTINUE);
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
