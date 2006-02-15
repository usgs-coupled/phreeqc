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
        totals.type = cxxNameDouble::ND_ELT_MOLES;
        master_activity.type = cxxNameDouble::ND_SPECIES_LA;
        species_gamma.type = cxxNameDouble::ND_SPECIES_GAMMA;
}

cxxSolution::cxxSolution(struct solution *solution_ptr)
        //
        // constructor for cxxSolution from struct solution
        //
: 
cxxNumKeyword(),
totals(solution_ptr->totals),
master_activity(solution_ptr->master_activity, solution_ptr->count_master_activity, cxxNameDouble::ND_SPECIES_LA),
species_gamma(solution_ptr->species_gamma, solution_ptr->count_species_gamma, cxxNameDouble::ND_SPECIES_GAMMA)
{
        int i;

        this->set_description(solution_ptr->description);
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

        // Totals filled in constructor, just save description and moles 

        // Isotopes
        for (i = 0; i < solution_ptr->count_isotopes; i++) {
                cxxIsotope iso(&solution_ptr->isotopes[i]);
                isotopes.push_back(iso);
        }

        // Master_activity in constructor
        // Species_gamma in constructor
}

cxxSolution::~cxxSolution()
{
}

struct solution *cxxSolution::cxxSolution2solution()
        //
        // Builds a solution structure from instance of cxxSolution 
        //
{

        struct solution *solution_ptr = solution_alloc();
        
        solution_ptr->description        = this->get_description();
        solution_ptr->n_user             = this->n_user;
        solution_ptr->n_user_end         = this->n_user_end;
        solution_ptr->new_def            = FALSE;
        solution_ptr->tc                 = this->tc;
        solution_ptr->ph                 = this->ph;
        solution_ptr->solution_pe        = this->pe;
        solution_ptr->mu                 = this->mu;
        solution_ptr->ah2o               = this->ah2o;
        solution_ptr->total_h            = this->total_h;
        solution_ptr->total_o            = this->total_o;
        solution_ptr->cb                 = this->cb;
        solution_ptr->mass_water         = this->mass_water;
        solution_ptr->total_alkalinity   = this->total_alkalinity;
        solution_ptr->density            = 1.0;
        solution_ptr->units              = moles_per_kilogram_string;
        solution_ptr->default_pe         = 0;
        // pe_data

        // totals
        solution_ptr->totals = (struct conc *) free_check_null(solution_ptr->totals);
        solution_ptr->totals = this->totals.conc();

        // master_activity
        solution_ptr->master_activity = (struct master_activity *) free_check_null(solution_ptr->master_activity);
        solution_ptr->master_activity = this->master_activity.master_activity();
        solution_ptr->count_master_activity = this->master_activity.size() + 1;

        // species_gamma
        solution_ptr->species_gamma = this->species_gamma.master_activity();
        solution_ptr->count_species_gamma = this->species_gamma.size();

        // isotopes
        solution_ptr->isotopes = (struct isotope *) free_check_null(solution_ptr->isotopes);
        solution_ptr->isotopes = cxxIsotope::list2isotope(this->isotopes);
        solution_ptr->count_isotopes = this->isotopes.size();

        return(solution_ptr);
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
        s_oss << "soln_ah2o=\""  << this->ah2o << "\"" << std::endl;

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
        this->totals.dump_xml(s_oss, indent + 1);
        /*
        {
        for (std::map <char *, double, CHARSTAR_LESS>::const_iterator it = totals.begin(); it != totals.end(); ++it) {
                s_oss << indent1;
                s_oss << "<soln_total";
                s_oss << " conc_desc=\"" << it->first << "\"";
                s_oss << " conc_moles=\"" << it->second << "\"" ;
                s_oss << "\">" << std::endl;
        }
        }
        */
        // master_activity map
        this->master_activity.dump_xml(s_oss, indent + 1);
        /*
        {
        for (std::map <char *, double>::const_iterator it = master_activity.begin(); it != master_activity.end(); ++it) {
                s_oss << indent1;
                s_oss << "<soln_m_a";
                s_oss << " m_a_desc=\"" << it->first << "\"" ;
                s_oss << " m_a_la=\"" << it->second << "\"" ;
                s_oss << "\">" << std::endl;
        }
        }
        */
        // species_gamma map
        this->species_gamma.dump_xml(s_oss, indent + 1);
        /*
        {
        for (std::map <char *, double>::const_iterator it = species_gamma.begin(); it != species_gamma.end(); ++it) {
                s_oss << indent1;
                s_oss << "<soln_s_g";
                s_oss << " m_a_desc=\"" << it->first << "\"" ;
                s_oss << " m_a_la=\"" << it->second << "\"" ;
                s_oss << "\">" << std::endl;
        }
        }
        */

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
        this->totals.dump_raw(s_oss, indent + 2);
        /*
        for (std::map <char *, double, CHARSTAR_LESS>::const_iterator it = totals.begin(); it != totals.end(); ++it) {
                s_oss << indent2;
                s_oss << it->first << "   " <<  it->second << std::endl;
        }
        */

        // master_activity map
        s_oss << indent1;
        s_oss << "-activities" << std::endl;
        this->master_activity.dump_raw(s_oss, indent + 2);
        /*
        {
                for (std::map <char *, double>::const_iterator it = master_activity.begin(); it != master_activity.end(); ++it) {
                        s_oss << indent2;
                        s_oss << it->first << "   " << it->second << std::endl;
                }
        }
        */
        // species_gamma map
        s_oss << indent1;
        s_oss << "-gammas" << std::endl;
        this->species_gamma.dump_raw(s_oss, indent + 2);
        /*
        {
                {
                        for (std::map <char *, double>::const_iterator it = species_gamma.begin(); it != species_gamma.end(); ++it) {
                                s_oss << indent2;
                                s_oss << it->first << "   " << it->second << std::endl;
                        }
                }
        }
        */
        // Isotopes
        s_oss << indent1;
        s_oss << "-Isotopes" << std::endl;
        {
                for (std::list<cxxIsotope>::const_iterator it = this->isotopes.begin(); it != isotopes.end(); ++it) {
                        it->dump_raw(s_oss, indent + 2);
                }
        }

        return;
}

void cxxSolution::read_raw(CParser& parser)
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

        // Read solution number and description
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
                        parser.error_msg("Unknown input in SOLUTION_RAW keyword.", CParser::OT_CONTINUE);
                        parser.error_msg(parser.line().c_str(), CParser::OT_CONTINUE);
                        break;

                case 0: // totals
                        if ( this->totals.read_raw(parser, next_char) != CParser::PARSER_OK) {
                                parser.incr_input_error();
                                parser.error_msg("Expected element name and moles for totals.", CParser::OT_CONTINUE);
                        }                       
                        opt_save = 0;
                        break;

                case 1: // activities
                        if ( this->master_activity.read_raw(parser, next_char) != CParser::PARSER_OK) {
                                parser.incr_input_error();
                                parser.error_msg("Expected species name and log activity for activities.", CParser::OT_CONTINUE);
                        }                               
                        opt_save = 1;
                        break;

                case 2: // gammas
                        if ( this->species_gamma.read_raw(parser, next_char) != CParser::PARSER_OK) {
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
                        opt_save = CParser::OPT_DEFAULT;
                        break;

                case 7: // ph
                        if (!(parser.get_iss() >> this->ph))
                        {
                                this->ph = 7.0;
                                parser.incr_input_error();
                                parser.error_msg("Expected numeric value for pH.", CParser::OT_CONTINUE);
                        }
                        ph_defined = true;
                        opt_save = CParser::OPT_DEFAULT;
                        break;

                case 8: // pe
                        if (!(parser.get_iss() >> this->pe))
                        {
                                this->pe = 4.0;
                                parser.incr_input_error();
                                parser.error_msg("Expected numeric value for pe.", CParser::OT_CONTINUE);
                        }
                        pe_defined = true;
                        opt_save = CParser::OPT_DEFAULT;
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
                        opt_save = CParser::OPT_DEFAULT;
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
                        opt_save = CParser::OPT_DEFAULT;
                        break;

                case 13: // total_h
                        if (!(parser.get_iss() >> this->total_h))
                        {
                                this->total_h = 111.1;
                                parser.incr_input_error();
                                parser.error_msg("Expected numeric value for total hydrogen.", CParser::OT_CONTINUE);
                        }
                        total_h_defined = true;
                        opt_save = CParser::OPT_DEFAULT;
                        break;

                case 14: // total_o
                        if (!(parser.get_iss() >> this->total_o))
                        {
                                this->total_o = 55.55;
                                parser.incr_input_error();
                                parser.error_msg("Expected numeric value for total oxygen.", CParser::OT_CONTINUE);
                        }
                        total_o_defined = true;
                        opt_save = CParser::OPT_DEFAULT;
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
                        opt_save = CParser::OPT_DEFAULT;
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
                        opt_save = CParser::OPT_DEFAULT;
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
                        opt_save = CParser::OPT_DEFAULT;
                        break;

                }
                if (opt == CParser::OPT_EOF || opt == CParser::OPT_KEYWORD) break;
        }
        // all members must be defined
        if (tc_defined == false) {
                parser.incr_input_error();
                parser.error_msg("Temp not defined for SOLUTION_RAW input.", CParser::OT_CONTINUE);
        }
        if (ph_defined == false) {
                parser.incr_input_error();
                parser.error_msg("pH not defined for SOLUTION_RAW input.", CParser::OT_CONTINUE);
        }
        if (pe_defined == false) {
                parser.incr_input_error();
                parser.error_msg("pe not defined for SOLUTION_RAW input.", CParser::OT_CONTINUE);
        }
        if (mu_defined == false) {
                parser.incr_input_error();
                parser.error_msg("Ionic strength not defined for SOLUTION_RAW input.", CParser::OT_CONTINUE);
        }
        if (ah2o_defined == false) {
                parser.incr_input_error();
                parser.error_msg("Activity of water not defined for SOLUTION_RAW input.", CParser::OT_CONTINUE);
        }
        if (total_h_defined == false) {
                parser.incr_input_error();
                parser.error_msg("Total hydrogen not defined for SOLUTION_RAW input.", CParser::OT_CONTINUE);
        }
        if (total_o_defined == false) {
                parser.incr_input_error();
                parser.error_msg("Total oxygen not defined for SOLUTION_RAW input.", CParser::OT_CONTINUE);
        }
        if (cb_defined == false) {
                parser.incr_input_error();
                parser.error_msg("Charge balance not defined for SOLUTION_RAW input.", CParser::OT_CONTINUE);
        }
        if (mass_water_defined == false) {
                parser.incr_input_error();
                parser.error_msg("Temp not defined for SOLUTION_RAW input.", CParser::OT_CONTINUE);
        }
        if (total_alkalinity_defined == false) {
                parser.incr_input_error();
                parser.error_msg("Total alkalinity not defined for SOLUTION_RAW input.", CParser::OT_CONTINUE);
        }
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
#include "Exchange.h"
#include "Surface.h"
#include "PPassemblage.h"
#include "KineticsCxx.h"
#include "SSassemblage.h"
#include "GasPhase.h"
#include "Reaction.h"
#include "Mix.h"
#include "Temperature.h"
#include <iostream>     // std::cout std::cerr
//#include <strstream>
#include <sstream>
#include <fstream>
void test_classes(void)
{
        int i;
        bool b(true);
        i = (int) b;
        for (i=0; i < count_solution; i++) {
                if (solution[i]->new_def == TRUE) {
                        cxxISolution sol(solution[i]);
                        solution[i] = (struct solution *) solution_free(solution[i]);
                        solution[i] = sol.cxxISolution2solution();
                        struct solution *soln_ptr;
                        soln_ptr = solution[i];
                        soln_ptr = solution[i];
                } else {
                        std::ostringstream oss;
                        cxxSolution sol(solution[i]);
                        solution[i] = (struct solution *) solution_free(solution[i]);
                        sol.dump_raw(oss, 0);

                        //std::fstream myfile("t");
                        //CParser cparser(myfile, std::cout, std::cerr);
                        cxxSolution sol1;
                        std::string keyInput = oss.str();
                        std::istringstream iss(keyInput);

                        CParser cparser(iss, oss, std::cerr);
                        //For testing, need to read line to get started
                        std::vector<std::string> vopts;
                        std::istream::pos_type next_char;
                        cparser.get_option(vopts, next_char);


                        sol1.read_raw(cparser);

                        solution[i] = sol1.cxxSolution2solution();
                }
        }
        for (i=0; i < count_exchange; i++) {
                if (exchange[i].new_def != TRUE) {
                        std::ostringstream oss;
                        cxxExchange ex(&(exchange[i]));
                        ex.dump_raw(oss, 0);
                        std::cerr << oss.str();

                        cxxExchange ex1;
                        std::string keyInput = oss.str();
                        std::istringstream iss(keyInput);

                        CParser cparser(iss, oss, std::cerr);
                        //For testing, need to read line to get started
                        std::vector<std::string> vopts;
                        std::istream::pos_type next_char;
                        cparser.get_option(vopts, next_char);

                        ex1.read_raw(cparser);

                        struct exchange *exchange_ptr = ex1.cxxExchange2exchange();
                        exchange_free(&exchange[i]);
                        exchange_copy(exchange_ptr, &exchange[i], exchange_ptr->n_user);
                        exchange_free(exchange_ptr);
                        free_check_null(exchange_ptr);
                }
        }
        for (i=0; i < count_surface; i++) {
                if (surface[i].new_def != TRUE) {
                        std::ostringstream oss;
                        cxxSurface ex(&(surface[i]));
                        ex.dump_raw(oss, 0);
                        std::cerr << oss.str();


                        cxxSurface ex1;
                        std::string keyInput = oss.str();
                        std::istringstream iss(keyInput);

                        CParser cparser(iss, oss, std::cerr);
                        //For testing, need to read line to get started
                        std::vector<std::string> vopts;
                        std::istream::pos_type next_char;
                        cparser.get_option(vopts, next_char);

                        ex1.read_raw(cparser);

                        struct surface *surface_ptr = ex1.cxxSurface2surface();
                        surface_free(&surface[i]);
                        surface_copy(surface_ptr, &surface[i], surface_ptr->n_user);
                        surface_free(surface_ptr);
                        free_check_null(surface_ptr);

                }

        }
        for (i=0; i < count_pp_assemblage; i++) {
                if (pp_assemblage[i].new_def != TRUE) {
                        std::ostringstream oss;
                        cxxPPassemblage ex(&(pp_assemblage[i]));
                        ex.dump_raw(oss, 0);
                        std::cerr << oss.str();


                        cxxPPassemblage ex1;
                        std::string keyInput = oss.str();
                        std::istringstream iss(keyInput);

                        CParser cparser(iss, oss, std::cerr);
                        //For testing, need to read line to get started
                        std::vector<std::string> vopts;
                        std::istream::pos_type next_char;
                        cparser.get_option(vopts, next_char);

                        ex1.read_raw(cparser);

                        struct pp_assemblage *pp_assemblage_ptr = ex1.cxxPPassemblage2pp_assemblage();
                        pp_assemblage_free(&pp_assemblage[i]);
                        pp_assemblage_copy(pp_assemblage_ptr, &pp_assemblage[i], pp_assemblage_ptr->n_user);
                        pp_assemblage_free(pp_assemblage_ptr);
                        free_check_null(pp_assemblage_ptr);

                }

        }
        for (i=0; i < count_kinetics; i++) {
                        std::ostringstream oss;
                        cxxKinetics ex(&(kinetics[i]));
                        ex.dump_raw(oss, 0);
                        std::cerr << oss.str();


                        cxxKinetics ex1;
                        std::string keyInput = oss.str();
                        std::istringstream iss(keyInput);

                        CParser cparser(iss, oss, std::cerr);
                        //For testing, need to read line to get started
                        std::vector<std::string> vopts;
                        std::istream::pos_type next_char;
                        cparser.get_option(vopts, next_char);


                        ex1.read_raw(cparser);

                        struct kinetics *kinetics_ptr = ex1.cxxKinetics2kinetics();
                        kinetics_free(&kinetics[i]);
                        kinetics_copy(kinetics_ptr, &kinetics[i], kinetics_ptr->n_user);
                        kinetics_free(kinetics_ptr);
                        free_check_null(kinetics_ptr);
        }
        for (i=0; i < count_s_s_assemblage; i++) {
                if (s_s_assemblage[i].new_def != TRUE) {
                        std::ostringstream oss;
                        cxxSSassemblage ex(&(s_s_assemblage[i]));
                        ex.dump_raw(oss, 0);
                        std::cerr << oss.str();


                        cxxSSassemblage ex1;
                        std::string keyInput = oss.str();
                        std::istringstream iss(keyInput);

                        CParser cparser(iss, oss, std::cerr);
                        //For testing, need to read line to get started
                        std::vector<std::string> vopts;
                        std::istream::pos_type next_char;
                        cparser.get_option(vopts, next_char);

                        ex1.read_raw(cparser);

                        struct s_s_assemblage *s_s_assemblage_ptr = ex1.cxxSSassemblage2s_s_assemblage();
                        s_s_assemblage_free(&s_s_assemblage[i]);
                        s_s_assemblage_copy(s_s_assemblage_ptr, &s_s_assemblage[i], s_s_assemblage_ptr->n_user);
                        s_s_assemblage_free(s_s_assemblage_ptr);
                        free_check_null(s_s_assemblage_ptr);

                }

        }
        for (i=0; i < count_gas_phase; i++) {
                if (gas_phase[i].new_def != TRUE) {
                        std::ostringstream oss;
                        cxxGasPhase ex(&(gas_phase[i]));
                        ex.dump_raw(oss, 0);
                        std::cerr << oss.str();


                        cxxGasPhase ex1;
                        std::string keyInput = oss.str();
                        std::istringstream iss(keyInput);

                        CParser cparser(iss, oss, std::cerr);
                        //For testing, need to read line to get started
                        std::vector<std::string> vopts;
                        std::istream::pos_type next_char;
                        cparser.get_option(vopts, next_char);

                        ex1.read_raw(cparser);

                        struct gas_phase *gas_phase_ptr = ex1.cxxGasPhase2gas_phase();
                        gas_phase_free(&gas_phase[i]);
                        gas_phase_copy(gas_phase_ptr, &gas_phase[i], gas_phase_ptr->n_user);
                        gas_phase_free(gas_phase_ptr);
                        free_check_null(gas_phase_ptr);

                }

        }
        for (i=0; i < count_irrev; i++) {
                        std::ostringstream oss;
                        cxxReaction ex(&(irrev[i]));
                        ex.dump_raw(oss, 0);
                        std::cerr << oss.str();


                        cxxReaction ex1;
                        std::string keyInput = oss.str();
                        std::istringstream iss(keyInput);

                        CParser cparser(iss, oss, std::cerr);
                        //For testing, need to read line to get started
                        std::vector<std::string> vopts;
                        std::istream::pos_type next_char;
                        cparser.get_option(vopts, next_char);

                        ex1.read_raw(cparser);
                        struct irrev *irrev_ptr = ex1.cxxReaction2irrev();

                        irrev_free(&irrev[i]);
                        irrev_copy(irrev_ptr, &irrev[i], irrev_ptr->n_user);

                        irrev_free(irrev_ptr);
                        free_check_null(irrev_ptr);

        }
        for (i=0; i < count_mix; i++) {
                        std::ostringstream oss;
                        cxxMix ex(&(mix[i]));
                        ex.dump_raw(oss, 0);
                        std::cerr << oss.str();


                        cxxMix ex1;
                        std::string keyInput = oss.str();
                        std::istringstream iss(keyInput);

                        CParser cparser(iss, oss, std::cerr);
                        //For testing, need to read line to get started
                        std::vector<std::string> vopts;
                        std::istream::pos_type next_char;
                        cparser.get_option(vopts, next_char);

                        ex1.read_raw(cparser);
                        struct mix *mix_ptr = ex1.cxxMix2mix();

                        mix_free(&mix[i]);
                        mix_copy(mix_ptr, &mix[i], mix_ptr->n_user);

                        mix_free(mix_ptr);
                        free_check_null(mix_ptr);

        }
        for (i=0; i < count_temperature; i++) {
                        std::ostringstream oss;
                        cxxTemperature ex(&(temperature[i]));
                        ex.dump_raw(oss, 0);
                        std::cerr << oss.str();


                        cxxTemperature ex1;
                        std::string keyInput = oss.str();
                        std::istringstream iss(keyInput);

                        CParser cparser(iss, oss, std::cerr);
                        //For testing, need to read line to get started
                        std::vector<std::string> vopts;
                        std::istream::pos_type next_char;
                        cparser.get_option(vopts, next_char);

                        ex1.read_raw(cparser);
                        struct temperature *temperature_ptr = ex1.cxxTemperature2temperature();

                        temperature_free(&temperature[i]);
                        temperature_copy(temperature_ptr, &temperature[i], temperature_ptr->n_user);

                        temperature_free(temperature_ptr);
                        free_check_null(temperature_ptr);

        }
} 
