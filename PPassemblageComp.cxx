// PPassemblageComp.cxx: implementation of the cxxPPassemblageComp class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
#pragma warning(disable : 4786)   // disable truncation warning (Only used by debugger)
#endif

#include "Utils.h"   // define first
#include "PPassemblageComp.h"
#define EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "phrqproto.h"
#include <cassert>     // assert
#include <algorithm>   // std::sort 

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cxxPPassemblageComp::cxxPPassemblageComp()
        //
        // default constructor for cxxPPassemblageComp 
        //
{
        name                             = NULL;
        add_formula                      = NULL;
        si                               = 0;
        moles                            = 0;
        delta                            = 0;
        initial_moles                    = 0;
        dissolve_only                    = false;
}

cxxPPassemblageComp::cxxPPassemblageComp(struct pure_phase *pure_phase_ptr)
        //
        // constructor for cxxPPassemblageComp from struct pure_phase
        //

{
        name                      = pure_phase_ptr->name;
        add_formula               = pure_phase_ptr->add_formula;
        si                        = pure_phase_ptr->si;
        moles                     = pure_phase_ptr->moles;
        delta                     = pure_phase_ptr->delta;
        initial_moles             = pure_phase_ptr->initial_moles;
        dissolve_only             = ( pure_phase_ptr->dissolve_only == TRUE);
}

cxxPPassemblageComp::~cxxPPassemblageComp()
{
}

struct phase *cxxPPassemblageComp::get_phase() {
        int i;
        return phase_bsearch(this->name, &i, FALSE);
}

struct pure_phase *cxxPPassemblageComp::cxxPPassemblageComp2pure_phase(std::list<cxxPPassemblageComp>& el)
        //
        // Builds pure_phase structure from of cxxPPassemblageComp 
        //
{
        struct pure_phase *pure_phase_ptr = (struct pure_phase *) PHRQ_malloc((size_t) (el.size() * sizeof(struct pure_phase)));
        if (pure_phase_ptr == NULL) malloc_error();

        int i = 0;
        for (std::list<cxxPPassemblageComp>::iterator it = el.begin(); it != el.end(); ++it) {
                pure_phase_ptr[i].phase                 =  it->get_phase();
                pure_phase_ptr[i].name                  =  it->name;
                pure_phase_ptr[i].add_formula           =  it->add_formula;
                pure_phase_ptr[i].si                    =  it->si;
                pure_phase_ptr[i].moles                 =  it->moles;
                pure_phase_ptr[i].delta                 =  it->delta;
                pure_phase_ptr[i].initial_moles         =  it->initial_moles;
                pure_phase_ptr[i].dissolve_only         =  (int) it->dissolve_only;
                i++;
        }
        return(pure_phase_ptr);
}

void cxxPPassemblageComp::dump_xml(std::ostream& s_oss, unsigned int indent)const
{
        //const char    ERR_MESSAGE[] = "Packing pure_phase message: %s, element not found\n";
        unsigned int i;
        s_oss.precision(DBL_DIG - 1);
        std::string indent0(""), indent1(""), indent2("");
        for(i = 0; i < indent; ++i) indent0.append(Utilities::INDENT);
        for(i = 0; i < indent + 1; ++i) indent1.append(Utilities::INDENT);
        for(i = 0; i < indent + 2; ++i) indent2.append(Utilities::INDENT);

        // Pure_Phase element and attributes

        s_oss << indent0 << "name=\"" << this->name << "\"" << std::endl;
        s_oss << indent0 << "add_formula=\"" << this->add_formula  << "\"" << std::endl;
        s_oss << indent0 << "si=\"" << this->si     << "\"" << std::endl;
        s_oss << indent0 << "moles=\"" << this->moles << "\"" << std::endl;
        s_oss << indent0 << "delta=\"" << this->delta << "\"" << std::endl;
        s_oss << indent0 << "initial_moles=\"" << this->initial_moles << "\"" << std::endl;
        s_oss << indent0 << "dissolve_only=\"" << this->dissolve_only  << "\"" << std::endl;

}

void cxxPPassemblageComp::dump_raw(std::ostream& s_oss, unsigned int indent)const
{
        //const char    ERR_MESSAGE[] = "Packing pure_phase message: %s, element not found\n";
        unsigned int i;
        s_oss.precision(DBL_DIG - 1);
        std::string indent0(""), indent1(""), indent2("");
        for(i = 0; i < indent; ++i) indent0.append(Utilities::INDENT);
        for(i = 0; i < indent + 1; ++i) indent1.append(Utilities::INDENT);
        for(i = 0; i < indent + 2; ++i) indent2.append(Utilities::INDENT);

        // Pure_Phase element and attributes

        if (this->name != NULL) s_oss << indent0 << "-name                  " << this->name << std::endl;
        if (this->add_formula != NULL) s_oss << indent0 << "-add_formula           " << this->add_formula  << std::endl;
        s_oss << indent0 << "-si                    " << this->si     << std::endl;
        s_oss << indent0 << "-moles                 " << this->moles << std::endl;
        s_oss << indent0 << "-delta                 " << this->delta << std::endl;
        s_oss << indent0 << "-initial_moles         " << this->initial_moles << std::endl;
        s_oss << indent0 << "-dissolve_only         " << this->dissolve_only  << std::endl;
}

void cxxPPassemblageComp::read_raw(CParser& parser)
{
        std::string str;
        
        static std::vector<std::string> vopts;
        if (vopts.empty()) {
                vopts.reserve(10);
                vopts.push_back("name");                     // 0                 
                vopts.push_back("add_formula");              // 1
                vopts.push_back("si");                       // 2
                vopts.push_back("moles");                    // 3
                vopts.push_back("delta");                    // 4
                vopts.push_back("initial_moles");            // 5     
                vopts.push_back("dissolve_only");            // 6
        }

        std::istream::pos_type ptr;
        std::istream::pos_type next_char;
        std::string token;
        int opt_save;

        opt_save = CParser::OPT_ERROR;
        bool name_defined(false); 
        bool si_defined(false);
        bool moles_defined(false); 
        bool delta_defined(false); 
        bool initial_moles_defined(false); 
        bool dissolve_only_defined(false);

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
                        opt = CParser::OPT_KEYWORD;
                        // Allow return to Exchange for more processing
                        //parser.error_msg("Unknown input in PURE_PHASE read.", CParser::OT_CONTINUE);
                        //parser.error_msg(parser.line().c_str(), CParser::OT_CONTINUE);
                        break;

                case 0: // name
                        if (!(parser.get_iss() >> str))
                        {
                                this->name = NULL;
                                parser.incr_input_error();
                                parser.error_msg("Expected string value for name.", CParser::OT_CONTINUE);
                        } else {
                                this->name = string_hsave(str.c_str());
                        }
                        name_defined = true;
                        break;

                case 1: // add_formula
                        if (!(parser.get_iss() >> str))
                        {
                                this->add_formula = NULL;
                                parser.incr_input_error();
                                parser.error_msg("Expected string value for add_formula.", CParser::OT_CONTINUE);
                        } else {
                                this->add_formula = string_hsave(str.c_str());
                        }
                        break;

                case 2: // si
                        if (!(parser.get_iss() >> this->si))
                        {
                                this->si = 0;
                                parser.incr_input_error();
                                parser.error_msg("Expected numeric value for si.", CParser::OT_CONTINUE);
                        }
                        si_defined = true;
                        break;

                case 3: // moles
                        if (!(parser.get_iss() >> this->moles))
                        {
                                this->moles = 0;
                                parser.incr_input_error();
                                parser.error_msg("Expected numeric value for moles.", CParser::OT_CONTINUE);
                        }
                        moles_defined = true;
                        break;

                case 4: // delta
                        if (!(parser.get_iss() >> this->delta))
                        {
                                this->delta = 0;
                                parser.incr_input_error();
                                parser.error_msg("Expected numeric value for delta.", CParser::OT_CONTINUE);
                        }
                        delta_defined = true;
                        break;

                case 5: // initial_moles
                        if (!(parser.get_iss() >> this->initial_moles))
                        {
                                this->initial_moles = 0;
                                parser.incr_input_error();
                                parser.error_msg("Expected numeric value for initial_moles.", CParser::OT_CONTINUE);
                        }
                        initial_moles_defined = true;
                        break;


                case 6: // dissolve_only
                        if (!(parser.get_iss() >> this->dissolve_only))
                        {
                                this->dissolve_only = false;
                                parser.incr_input_error();
                                parser.error_msg("Expected boolean value for dissolve_only.", CParser::OT_CONTINUE);
                        }
                        dissolve_only_defined = true;
                        break;
                }
                if (opt == CParser::OPT_EOF || opt == CParser::OPT_KEYWORD) break;
        }
        // members that must be defined
        if (name_defined == false) {
                parser.incr_input_error();
                parser.error_msg("Name not defined for PPassemblageComp input.", CParser::OT_CONTINUE);
        }
        if (si_defined == false) {
                parser.incr_input_error();
                parser.error_msg("Si not defined for PPassemblageComp input.", CParser::OT_CONTINUE);
        }
        if (moles_defined == false) {
                parser.incr_input_error();
                parser.error_msg("Moles not defined for PPassemblageComp input.", CParser::OT_CONTINUE);
        }
        if (delta_defined == false) {
                parser.incr_input_error();
                parser.error_msg("Delta not defined for PPassemblageComp input.", CParser::OT_CONTINUE);
        }
        if (initial_moles_defined == false) {
                parser.incr_input_error();
                parser.error_msg("Initial_moles not defined for PPassemblageComp input.", CParser::OT_CONTINUE);
        }
        if (dissolve_only_defined == false) {
                parser.incr_input_error();
                parser.error_msg("Dissolve_only not defined for PPassemblageComp input.", CParser::OT_CONTINUE);
        }
}

