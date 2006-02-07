// NameDouble.cxx: implementation of the cxxNameDouble class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
#pragma warning(disable : 4786)   // disable truncation warning (Only used by debugger)
#endif

#include "Utils.h"   // define first
#include "Conc.h"
#include "NameDouble.h"
#define EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "phrqproto.h"
#include <cassert>     // assert
#include <algorithm>   // std::sort 
#include <map>   // std::sort 

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cxxNameDouble::cxxNameDouble() 
        //
        // default constructor for cxxNameDouble 
        //
{
}

cxxNameDouble::cxxNameDouble(struct elt_list *elt_list_ptr)
        //
        // constructor for cxxNameDouble from list of elt_list
        //
{
	int i;
	for (i = 0; elt_list_ptr[i].elt != NULL; i++) {
		(*this)[elt_list_ptr[i].elt->name] = elt_list_ptr[i].coef;
	}
	this->type = ND_ELT_MOLES;
}
cxxNameDouble::cxxNameDouble(struct conc *tots)
        //
        // constructor for cxxNameDouble from list of elt_list
        //
{
	int i;
        for (i = 0; tots[i].description != NULL; i++) {
                (*this)[tots[i].description] = tots[i].moles;
        }
	this->type = ND_ELT_MOLES;
}
cxxNameDouble::cxxNameDouble(struct master_activity *ma, int count, cxxNameDouble::ND_TYPE)
        //
        // constructor for cxxNameDouble from list of elt_list
        //
{
	int i;
        for (i = 0; i < count; i++) {
                if (ma[i].description == NULL) continue;
                (*this)[ma[i].description] = ma[i].la;
        }
	this->type = ND_SPECIES_LA;
}
cxxNameDouble::~cxxNameDouble()
{
}

struct elt_list *cxxNameDouble::elt_list()
        //
        // Builds a exch_comp structure from instance of cxxNameDouble 
        //
{
	assert (this->type == cxxNameDouble::ND_ELT_MOLES);
	struct elt_list *elt_list_ptr = (struct elt_list *) PHRQ_malloc((size_t)((this->size() + 1) *sizeof(struct elt_list)));
	int i = 0;
	for (iterator it = this->begin(); it != this->end(); ++it) {
		elt_list_ptr[i].elt = element_store(it->first);
		elt_list_ptr[i].coef = it->second;
	}
	elt_list_ptr[i].elt = NULL;
	elt_list_ptr[i].coef = 0;
	return(elt_list_ptr);
}

struct master_activity *cxxNameDouble::master_activity()const
        //
        // Builds a list of master_activity structures from instance of cxxNameDouble 
        //
{
	int i = 0;
	assert (this->type == cxxNameDouble::ND_SPECIES_LA || this->type == cxxNameDouble::ND_SPECIES_GAMMA);
	struct master_activity *master_activity_ptr = NULL;
	switch ((*this).type) {
	case cxxNameDouble::ND_SPECIES_LA:
		master_activity_ptr= (struct master_activity *) PHRQ_malloc((size_t) (((*this).size() + 1) * sizeof(struct master_activity)));
		if (master_activity_ptr == NULL) malloc_error();
		for (const_iterator it = (*this).begin(); it != (*this).end(); it++) {
			master_activity_ptr[i].description = (char *)it->first;
			master_activity_ptr[i].la = it->second;
			i++;
		}
		master_activity_ptr[i].description = NULL;
		break;
	case cxxNameDouble::ND_SPECIES_GAMMA:
		if ((*this).size() > 0) {
			master_activity_ptr = (struct master_activity *) PHRQ_malloc((size_t) (((*this).size()) * sizeof(struct master_activity)));
			if (master_activity_ptr == NULL) malloc_error();
			for (const_iterator it = (*this).begin(); it != (*this).end(); it++) {
				master_activity_ptr[i].description = (char *)it->first;
				master_activity_ptr[i].la = it->second;
				i++;
			}
		}
		break;
	case cxxNameDouble::ND_ELT_MOLES:
		break;
	}
	return(master_activity_ptr);
}

struct conc *cxxNameDouble::conc()const
	// for Solutions, not ISolutions
	// takes a map of (elt name, moles)
	// returns list of conc structures
{
	struct conc *c;
	assert (this->type == cxxNameDouble::ND_ELT_MOLES);
	c = (struct conc *) PHRQ_malloc((size_t) (((*this).size() + 1) * sizeof(struct conc)));
	if (c == NULL) malloc_error();
	int i = 0;
	for (const_iterator it = (*this).begin(); it != (*this).end(); ++it) {
		c[i].description         = (char *)it->first;
		c[i].moles               = it->second;
		c[i].input_conc          = it->second;
		c[i].units               = NULL;
		c[i].equation_name       = NULL;
		c[i].phase_si            = 0.0;
		c[i].n_pe                = 0;
		c[i].as                  = NULL;
		c[i].gfw                 = 0.0;
		//c[i].skip                = 0;
		c[i].phase               = NULL;
		i++;
	}			
	c[i].description = NULL;
	return(c);
}

void cxxNameDouble::dump_xml(std::ostream& s_oss, unsigned int indent)const
{
        //const char    ERR_MESSAGE[] = "Packing exch_comp message: %s, element not found\n";
        unsigned int i;
        s_oss.precision(DBL_DIG - 1);
        std::string indent0(""), indent1("");
        for(i = 0; i < indent; ++i) indent0.append(Utilities::INDENT);
        for(i = 0; i < indent + 1; ++i) indent1.append(Utilities::INDENT);
	std::string xmlElement, xmlAtt1, xmlAtt2;

	switch ((*this).type) {
	case cxxNameDouble::ND_SPECIES_LA:
		xmlElement = "<soln_m_a ";
		xmlAtt1 = " m_a_desc=\"";
		xmlAtt1 = " m_a_la=\"";
		break;
	case cxxNameDouble::ND_SPECIES_GAMMA:
		xmlElement = "<soln_s_g ";
		xmlAtt1 = " m_a_desc=\"";
		xmlAtt1 = " m_a_la=\"";
		break;
	case cxxNameDouble::ND_ELT_MOLES:
		xmlElement = "<soln_total ";
		xmlAtt1 = " conc_desc=\"";
		xmlAtt1 = " conc_moles=\"";
		break;
	}
	
	for (const_iterator it = (*this).begin(); it != (*this).end(); it++) {
		s_oss << indent0;
		s_oss << xmlElement << xmlAtt1 << it->first << xmlAtt2 << it->second << "/>" << std::endl;
	}
}

void cxxNameDouble::dump_raw(std::ostream& s_oss, unsigned int indent)const
{
        //const char    ERR_MESSAGE[] = "Packing exch_comp message: %s, element not found\n";
        unsigned int i;
        s_oss.precision(DBL_DIG - 1);
        std::string indent0("");
        for(i = 0; i < indent; ++i) indent0.append(Utilities::INDENT);
	
	for (const_iterator it = (*this).begin(); it != (*this).end(); it++) {
		s_oss << indent0;
		s_oss << it->first << "   " << it->second << std::endl;
	}
}



#ifdef SKIP
void cxxNameDouble::read_raw(CParser& parser)
{
        static std::vector<std::string> vopts;
        if (vopts.empty()) {
                vopts.reserve(21);
                vopts.push_back("totals");                    // 0 
                vopts.push_back("activities");                // 1 
                vopts.push_back("gammas");                    // 2 
                vopts.push_back("isotopes");                  // 3 
                vopts.push_back("temp");                      // 4 
                vopts.push_back("tc");                        // 5 
                vopts.push_back("temperature");               // 6 
                vopts.push_back("ph");                        // 7 
                vopts.push_back("pe");                        // 8 
                vopts.push_back("mu");                        // 9 
                vopts.push_back("ionic_strength");            // 10
                vopts.push_back("ah2o");                      // 11
                vopts.push_back("activity_water");            // 12
                vopts.push_back("total_h");                   // 13
                vopts.push_back("total_o");                   // 14
                vopts.push_back("mass_water");                // 15
                vopts.push_back("mass_h2o");                  // 16
                vopts.push_back("total_alkalinity");          // 17
                vopts.push_back("total_alk");                 // 18
                vopts.push_back("cb");                        // 19
                vopts.push_back("charge_balance");            // 20
        }

        std::istream::pos_type ptr;
        std::istream::pos_type next_char;
        std::string token;
        int opt_save;

        // Read exch_comp number and description
        this->read_number_description(parser);

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
