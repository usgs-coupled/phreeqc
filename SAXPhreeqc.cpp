// SAXPhreeqc.cpp
#ifdef _DEBUG
#pragma warning(disable : 4786)   // disable truncation warning (Only used by debugger)
#endif

#include <float.h>                         // DBL_DIG
#include <stdio.h>                         // sprintf
#include <wctype.h>                        // iswspace

#include <cassert>                         // assert
#include <strstream>                       // std::ostrstream
#include <iostream>                        // std::cerr

#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#endif

#include <xercesc/framework/StdInInputSource.hpp>
#include <xercesc/framework/MemoryManager.hpp>
#include <xercesc/framework/XMLGrammarPool.hpp>
#include <xercesc/internal/XMLGrammarPoolImpl.hpp>
#include <xercesc/internal/MemoryManagerImpl.hpp>
#include <xercesc/validators/common/Grammar.hpp>
#include <xercesc/parsers/SAXParser.hpp>

#include <xercesc/parsers/SAXParser.hpp>           // SAXParser
#include <xercesc/sax/AttributeList.hpp>           // AttributeList
#include <xercesc/util/XMLUniDefs.hpp>             // Unicode definitions
#include <xercesc/framework/MemBufInputSource.hpp> // MemBufInputSource
#include <xercesc/framework/MemoryManager.hpp>
#include <xercesc/internal/MemoryManagerImpl.hpp>
#include <xercesc/util/PlatformUtils.hpp>          // XMLPlatformUtils::getCurrentMillis
#include <xercesc/util/XMLUni.hpp>

#include "SAXPhreeqc.h"                    // SAX_ functions
#include "SaxPhreeqcHandlers.h"            // SaxPhreeqcHandlers

//XERCES_CPP_NAMESPACE_USE
#define xns XERCES_CPP_NAMESPACE

///static char buffer[300000];
///static std::ostrstream s_oss(buffer, 300000);        // must never go out of scope
static std::ostrstream s_oss;        // must never go out of scope
static bool s_bSysIsOpen = false;    // must never go out of scope
//static SaxPhreeqcHandlers s_handler; // one and only instance
static SaxPhreeqcHandlers *s_handler; // one and only instance

// global.h defines
#define OK 1
#define ERROR 0
#define TRUE 1
#define FALSE 0
#define MAX_LENGTH 200          /* maximum number of characters component name */
#define CONTINUE 0
#define STOP 1

#ifdef USE_LONG_DOUBLE
#define LDBLE long double
#else
#define LDBLE double
#endif

// global.h defs
struct solution {
  int new_def;
  int n_user;
  int n_user_end;
  char *description;
  double tc;
  double ph;
  double solution_pe;
  double mu;
  double ah2o;
  double density;
  LDBLE total_h;
  LDBLE total_o;
  LDBLE cb;
  LDBLE mass_water;
  LDBLE total_alkalinity;
  LDBLE total_co2;
  char *units;
  struct pe_data *pe; 
  int default_pe;
  struct conc *totals;
  struct master_activity *master_activity;
  int count_master_activity;
  int count_isotopes;
  struct isotope *isotopes;
  struct master_activity *species_gamma;
  int count_species_gamma;
};
struct master {                       /* list of name and number of elements in an equation */
  int in;                       /* TRUE if in model, FALSE if out, REWRITE if other mb eq */
  int number;                   /* sequence number in list of masters */
  int last_model;               /* saved to determine if model has changed */
  int type;                     /* AQ or EX */ 
  int primary;                  /* TRUE if master species is primary */
  double coef;                  /* coefficient of element in master species */
  LDBLE total;                 /* total concentration for element or valence state */
  double isotope_ratio;
  double isotope_ratio_uncertainty;
  int isotope;
  double total_primary;
  /*	double la;           */       /* initial guess of master species log activity */
  struct element *elt;          /* element structure */
  double alk;                   /* alkalinity of species */
  double gfw;                   /* default gfw for species */
  char *gfw_formula;            /* formula from which to calcuate gfw */
  struct unknown *unknown;      /* pointer to unknown structure */
  struct species *s;            /* pointer to species structure */
  struct reaction *rxn_primary; /* reaction writes master species in terms of primary 
				   master species */
  struct reaction *rxn_secondary; /* reaction writes master species in terms of secondary 
				     master species */
  struct reaction **pe_rxn;        /* e- written in terms of redox couple (or e-), points
				      to location */
};

struct element {
  char *name;                   /* element name */
  /*	int in; */
  struct master *master;
  struct master *primary;
  double gfw;
};
struct pe_data {
	char *name;
	struct reaction *rxn;
};

// extern routines
extern "C" {
#include "phqalloc.h"
  void *free_check_null(void *ptr);
  int pe_data_store (struct pe_data **pe, const char *token);
  struct phase *phase_bsearch (char *ptr, int *j, int print);
  struct solution *solution_alloc(void);
  struct solution *solution_bsearch(int k, int *n, int print);
  int solution_free (struct solution *solution_ptr);
  void space (void **ptr, int i, int *max, int struct_size);
  char * string_duplicate (const char *token);
  char *string_hsave (const char *str);
  int error_msg (const char *err_str, const int stop);
  struct master *master_bsearch (const char *ptr);
  void malloc_error(void);
}

// globals
extern "C" int input_error;
extern "C" char error_string[10*MAX_LENGTH];
extern "C" struct solution **solution;
extern "C" int count_solution;
extern "C" int max_solution;

class Initializer
{
public:
  Initializer(){ 
	xns::XMLPlatformUtils::Initialize();
	s_handler = new SaxPhreeqcHandlers();
  }
};

static Initializer sInit;  // initialize xerces once

extern "C" void SAX_StartSystem()
{

  assert(!s_bSysIsOpen);     // system already open and has not been closed
#if defined(_DEBUG) 
	int tmpDbgFlag;

	/*
	* Set the debug-heap flag to keep freed blocks in the
	* heap's linked list - This will allow us to catch any
	* inadvertent use of freed memory
	*/
	tmpDbgFlag = _CrtSetDbgFlag(_CRTDBG_REPORT_FLAG);
	tmpDbgFlag |= _CRTDBG_DELAY_FREE_MEM_DF;
	tmpDbgFlag |= _CRTDBG_LEAK_CHECK_DF;
	tmpDbgFlag |= _CRTDBG_CHECK_ALWAYS_DF;
	_CrtSetDbgFlag(tmpDbgFlag);
#endif 
  // init stream
  s_oss.freeze(false);
  s_oss.seekp(0);
  
  // write stream
  // s_oss << "<system>";
s_oss << "<?xml version=\"1.0\" ?>";
s_oss << "<phast_state nx=\"2\">";
s_oss << "  <system system_number=\"1\">";
/*
s_oss << "    <solution ";
s_oss << "        soln_new_def=\"0\"";
s_oss << "	soln_n_user=\"2\" ";
s_oss << "	soln_n_user_end=\"2\" ";
s_oss << "	soln_description=\"Test solution\" ";
s_oss << "	soln_tc=\"25.0\"";
s_oss << "	soln_ph=\"6.3\"";
s_oss << "	soln_solution_pe=\"4.5\"";
s_oss << "	soln_mu=\"0.0011\"";
s_oss << "	soln_ah2o=\"0.995\"";
s_oss << "	soln_density=\"1.0\"";
s_oss << "	soln_total_h=\"111.1\"";
s_oss << "	soln_total_o=\"55.5\"";
s_oss << "	soln_cb=\"1e-15\"";
s_oss << "	soln_mass_water=\"1.0\"";
s_oss << "	soln_total_alkalinity=\"0.001\"";
s_oss << "	soln_total_co2=\"0.002\"";
s_oss << "	soln_units=\"mg/L\"";
s_oss << "	soln_default_pe=\"0\"";
s_oss << "	soln_count_master_activity=\"2\"";
s_oss << "	soln_count_isotopes=\"0\"";
s_oss << "	soln_count_species_gamma=\"0\">";
s_oss << "	<soln_pe soln_pe_name=\"Fe(2)/Fe(3)\"/>";
s_oss << "	<soln_total";
s_oss << "	     conc_description=\"Ca\"";
s_oss << "	     conc_skip=\"0\"";
s_oss << "	     conc_moles=\"0.001\"";
s_oss << "	     conc_input_conc=\"0.001\"";
s_oss << "	     conc_equation_name=\"Fe(2)/Fe(3)\"";
s_oss << "	     conc_phase_si=\"0.0\"";
s_oss << "	     conc_n_pe=\"0\"";
s_oss << "	     conc_as=\"mg/L\"";
s_oss << "	     conc_gfw=\"40.08\"/>";
s_oss << "	<soln_total";
s_oss << "	     conc_description=\"Na\"";
s_oss << "	     conc_skip=\"0\"";
s_oss << "	     conc_moles=\"0.001\"";
s_oss << "	     conc_input_conc=\"0.001\"";
s_oss << "	     conc_equation_name=\"Fe(2)/Fe(3)\"";
s_oss << "	     conc_phase_si=\"0.0\"";
s_oss << "	     conc_n_pe=\"0\"";
s_oss << "	     conc_as=\"mg/L\"";
s_oss << "	     conc_gfw=\"40.08\"/>";
s_oss << "         <soln_master_activity m_a_description=\"Ca+2\" m_a_la=\"-3.0\"/>";
s_oss << "         <soln_master_activity m_a_description=\"Cl-\" m_a_la=\"-3.0\"/>";
s_oss << "         <soln_master_activity m_a_description=\"C\" m_a_la=\"-3.0\"/>";
s_oss << "         <soln_species_gamma m_a_description=\"CO3-2\" m_a_la=\".1\"/>";
s_oss << "         <soln_species_gamma m_a_description=\"HCO3-\" m_a_la=\".2\"/>";
s_oss << "         <soln_species_gamma m_a_description=\"CO2\" m_a_la=\".3\"/>";
s_oss << "         <soln_isotope ";
s_oss << "	      iso_isotope_number=\"2\"";
s_oss << "	      iso_elt_name=\"H\"";
s_oss << "	      iso_isotope_name=\"D\"";
s_oss << "	      iso_total=\"1e-6\"";
s_oss << "	      iso_ratio=\"11.\"";
s_oss << "	      iso_ratio_uncertainty=\"0.5\"";
s_oss << "	      iso_coef=\"2.2\"/>";
s_oss << "         <soln_isotope ";
s_oss << "	      iso_isotope_number=\"3\"";
s_oss << "	      iso_elt_name=\"H\"";
s_oss << "	      iso_isotope_name=\"T\"";
s_oss << "	      iso_total=\"1e-14\"";
s_oss << "	      iso_ratio=\"3.1\"";
s_oss << "	      iso_ratio_uncertainty=\"0.25\"";
s_oss << "	      iso_coef=\"2.1\"/>";
s_oss << "    </solution>";
s_oss << "  </system>";
s_oss << "</phast_state>";
*/
  s_bSysIsOpen = true;

}
char * stringify(char * string) {
	if (string == NULL) return("");
	return(string);
}
extern "C" int SAX_AddSolution(struct solution* solution_ptr)
{
  const char    ERR_MESSAGE[] = "Packing solution message: %s, element not found\n";
  int i;  
  assert(s_bSysIsOpen);      // must call SAX_StartSystem first
	// Solution element and attributes
	s_oss << "    <solution " << std::endl;
	s_oss << "      soln_new_def=\"" << solution_ptr->new_def << "\"" << std::endl;
	s_oss << "      soln_n_user=\"" << solution_ptr->n_user << "\" " << std::endl;
	s_oss << "      soln_n_user_end=\"" << solution_ptr->n_user_end << "\"" << std::endl;
	s_oss << "      soln_description=\"" << solution_ptr->description << "\"" << std::endl;
	s_oss << "      soln_tc=\"" << solution_ptr->tc << "\"" << std::endl;
	s_oss << "      soln_ph=\"" << solution_ptr->ph << "\"" << std::endl;
	s_oss << "      soln_solution_pe=\"" << solution_ptr->solution_pe << "\"" << std::endl;
	s_oss << "      soln_mu=\"" << solution_ptr->mu << "\"" << std::endl;
	s_oss << "      soln_ah2o=\""	 << solution_ptr->ah2o << "\"" << std::endl;
	s_oss << "      soln_density=\"" << solution_ptr->density << "\"" << std::endl;
	s_oss << "      soln_total_h=\"" << solution_ptr->total_h << "\"" << std::endl;
	s_oss << "      soln_total_o=\"" << solution_ptr->total_o << "\"" << std::endl;
	s_oss << "      soln_cb=\"" << solution_ptr->cb << "\"" << std::endl;
	s_oss << "      soln_mass_water=\"" << solution_ptr->mass_water << "\"" << std::endl;
	s_oss << "      soln_total_alkalinity=\"" << solution_ptr->total_alkalinity << "\"" << std::endl;
	s_oss << "      soln_total_co2=\"" << solution_ptr->total_co2 << "\"" << std::endl;
	s_oss << "      soln_units=\"" << solution_ptr->units << "\"" << std::endl;
	s_oss << "      soln_default_pe=\"" << solution_ptr->default_pe << "\"" << std::endl;
	s_oss << "      soln_count_master_activity=\"" << solution_ptr->count_master_activity << "\"" << std::endl;
	s_oss << "      soln_count_isotopes=\"" << solution_ptr->count_isotopes << "\"" << std::endl;
	s_oss << "      soln_count_species_gamma=\"" << solution_ptr->count_species_gamma << "\">" << std::endl;
	// end of solution attributes
	// pe structures
	for (i=0; solution_ptr->pe[i].name != NULL; i++) {
		s_oss << "      <soln_pe soln_pe_name=\"" << solution_ptr->pe[i].name << "\"/>"<< std::endl; 
	}
	// soln_total conc structures
	for (i=0; solution_ptr->totals[i].description != NULL; i++) {
		struct conc *c = &(solution_ptr->totals[i]);
		s_oss << "      <soln_total " << std::endl;
		s_oss << "        conc_description=\"" << c->description << "\"" << std::endl;
		s_oss << "        conc_skip=\"" << c->skip << "\"" << std::endl;
		s_oss << "        conc_moles=\"" << c->moles << "\"" << std::endl;
		s_oss << "        conc_input_conc=\"" << c->input_conc << "\"" << std::endl;
		s_oss << "        conc_equation_name=\"" << stringify(c->equation_name) << "\"" << std::endl;
		s_oss << "        conc_phase_si=\"" << c->phase_si << "\"" << std::endl;
		s_oss << "        conc_n_pe=\"" << c->n_pe << "\"" << std::endl;
		s_oss << "        conc_as=\"" << stringify(c->as) << "\"" << std::endl;
		s_oss << "        conc_gfw=\"" << c->gfw << "\"/>" << std::endl; 
	}
	// master_activity, master_activity structure
	for (i=0; i < solution_ptr->count_master_activity; i++) {
		s_oss << "      <soln_master_activity m_a_description=\"" << stringify(solution_ptr->master_activity[i].description) << "\" m_a_la=\"" << solution_ptr->master_activity[i].la << "\"/>" << std::endl;
	}
	// species_gamma, mater_activity structure
	for (i=0; i < solution_ptr->count_species_gamma; i++) {
		s_oss << "      <soln_species_gamma m_a_description=\"" << stringify(solution_ptr->species_gamma[i].description) << "\" m_a_la=\"" << solution_ptr->species_gamma[i].la << "\"/>" << std::endl;
	}
	// isotopes, isotope structure
	for (i=0; solution_ptr->count_isotopes; i++) {
		s_oss << "      <soln_isotope " << std::endl;
		s_oss << "        iso_isotope_number=\"" << solution_ptr->isotopes[i].isotope_number << "\"" << std::endl;
		s_oss << "        iso_elt_name=\"" << solution_ptr->isotopes[i].elt_name << "\"" << std::endl;
		s_oss << "        iso_isotope_name=\"" << solution_ptr->isotopes[i].isotope_name << "\"" << std::endl;
		s_oss << "        iso_total=\"" << solution_ptr->isotopes[i].total << "\"" << std::endl;
		s_oss << "        iso_ratio=\"" << solution_ptr->isotopes[i].ratio << "\"" << std::endl;
		s_oss << "        iso_ratio_uncertainty=\"" << solution_ptr->isotopes[i].ratio_uncertainty << "\"" << std::endl;
		s_oss << "        iso_coef=\"" << solution_ptr->isotopes[i].coef << "\"" << std::endl;
	}
	// End of solution
	s_oss << "    </solution>" << std::endl;

  return(OK);
}

extern "C" void SAX_EndSystem()
{
  assert(s_bSysIsOpen);      // must call SAX_StartSystem first

  s_oss << "  </system>" << std::endl;
  s_oss << "</phast_state>" << std::endl; 
  s_oss << '\0';
  
  s_bSysIsOpen = false;
  return;
}

extern "C" int SAX_GetXMLLength()
{
  assert(!s_bSysIsOpen);      // must call SAX_EndSystem first
  return s_oss.pcount();
}

extern "C" char* SAX_GetXMLStr()
{
  assert(!s_bSysIsOpen);      // must call SAX_EndSystem first
  return s_oss.str();
}

extern "C" int SAX_UnpackSolutions(void* pvBuffer, int buf_size)
{
  //  Create MemBufferInputSource from the buffer containing the XML
  //  statements.
	xns::MemBufInputSource memBufIS((const XMLByte*)pvBuffer, buf_size, "solution_id", false);
	fprintf(stderr,"%s", (char *) pvBuffer);
  //
  //  Create a SAX parser object.
  //
	xns::SAXParser parser;
	parser.setValidationScheme(xns::SAXParser::Val_Always);
	parser.setDoNamespaces(false);
	parser.setDoSchema(false);

  
  //
  //  Create the handler object and install it as the document and error
  //  handler for the parser. Then parse the MemBufferInputSource and 
  //  catch any exceptions that propogate out
  //
  try
    {
      unsigned long duration = 0;
      parser.setDocumentHandler(s_handler);
      parser.setErrorHandler(s_handler);
      
      int t = 0;
      for (; t < 1; ++t)
	{
		const unsigned long startMillis = xns::XMLPlatformUtils::getCurrentMillis();
		parser.parse(memBufIS);
		const unsigned long endMillis = xns::XMLPlatformUtils::getCurrentMillis();
		duration += endMillis - startMillis;
	}
      std::cerr << "\nSaxParse time = " << (float)duration/t << " millis\n";
    }
  catch (const xns::SAXException& toCatch)
    {
      char* psz = xns::XMLString::transcode(toCatch.getMessage());
      input_error++;
      sprintf(error_string, "SAX_UnpackSolutions: %s\n", psz);
      error_msg(error_string, CONTINUE);
      delete[] psz;
      return ERROR;
    }
  catch (const xns::XMLException& toCatch)
    {
      char* psz = xns::XMLString::transcode(toCatch.getMessage());
      input_error++;
      sprintf(error_string, "SAX_UnpackSolutions: %s\n", psz);
      error_msg(error_string, CONTINUE);
      delete[] psz;
      return ERROR;
    }
  catch (...)
    {
      input_error++;
      sprintf(error_string,"SAX_UnpackSolutions: %s\n", "Unknown error occured.");
      error_msg(error_string, CONTINUE);
      return ERROR;
    }
  
  return(OK);
}


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

// element names

SaxPhreeqcHandlers::SaxPhreeqcHandlers()
  : eltType(typeNULL), attType(attNULL), totals(), acts(), solution_ptr(NULL)
{
	int i;

	int count_elementInfo, count_attributeInfo;
	struct mapElementInfo {char *key; enum elementType type;};
	struct mapAttributeInfo {enum attributeType type; char *key;};
	struct mapElementInfo elementInfo[] = {
		{"phast_state", typePHAST_STATE},
		{"system", typeSYSTEM},
		{"solution", typeSOLUTION},
		{"soln_pe", typeSOLN_PE},
		{"soln_total", typeSOLN_TOTAL},
		{"soln_master_activity", typeSOLN_MASTER_ACTIVITY},
		{"soln_isotope", typeSOLN_ISOTOPE},
		{"soln_species_gamma", typeSOLN_SPECIES_GAMMA}
	};
	count_elementInfo = sizeof(elementInfo)/sizeof(struct mapElementInfo);
	struct mapAttributeInfo	attributeInfo[] = {	
		// Solution structure
	  {attSOLN_new_def, "soln_new_def"},
	  {attSOLN_n_user, "soln_n_user"},
	  {attSOLN_n_user_end, "soln_n_user_end"},
	  {attSOLN_description, "soln_description"},
	  {attSOLN_tc, "soln_tc"},
	  {attSOLN_ph, "soln_ph"},
	  {attSOLN_solution_pe, "soln_solution_pe"},
	  {attSOLN_mu, "soln_mu"},
	  {attSOLN_ah2o, "soln_ah2o"},
	  {attSOLN_density, "soln_density"},
	  {attSOLN_total_h, "soln_total_h"},
	  {attSOLN_total_o, "soln_total_o"},
	  {attSOLN_cb, "soln_cb"},
	  {attSOLN_mass_water, "soln_mass_water"},
	  {attSOLN_total_alkalinity, "soln_total_alkalinity"},
	  {attSOLN_total_co2, "soln_total_co2"},
	  {attSOLN_units, "soln_units"},
	  {attSOLN_default_pe, "soln_default_pe"},
	  {attSOLN_count_master_activity, "soln_count_master_activity"},
	  {attSOLN_count_isotopes, "soln_count_isotopes"},
	  {attSOLN_count_species_gamma, "soln_count_species_gamma"},
	  {attSOLN_PE_name, "soln_pe_name"},
		// master_activity structure
	  {attM_A_description, "m_a_description"},
	  {attM_A_la, "m_a_la"},
		// isotope structure
	  {attISO_isotope_number, "iso_isotope_number"},
	  {attISO_elt_name, "iso_elt_name"},
	  {attISO_isotope_name, "iso_isotope_name"},
	  {attISO_total, "iso_total"},
	  {attISO_ratio, "iso_ratio"},
	  {attISO_ratio_uncertainty, "iso_ratio_uncertainty"},
	  {attISO_x_ratio_uncertainty, "iso_x_ratio_uncertainty"},
	  {attISO_coef, "iso_coef"},
		// conc structure
	  {attCONC_description, "conc_description"},
	  {attCONC_skip, "conc_skip"},
	  {attCONC_moles, "conc_moles"},
	  {attCONC_input_conc, "conc_input_conc"},
	  {attCONC_equation_name, "conc_equation_name"},
	  {attCONC_phase_si, "conc_phase_si"},
	  {attCONC_n_pe, "conc_n_pe"},
	  {attCONC_as, "conc_as"},
	  {attCONC_gfw, "conc_gfw"},

	};
	count_attributeInfo = sizeof(attributeInfo)/sizeof(struct mapAttributeInfo);
	for (i = 0; i < count_elementInfo; i++) {
		this->mapXMLCh2Type[xns::XMLString::transcode(elementInfo[i].key)] = elementInfo[i].type;
	}
	for (i = 0; i < count_attributeInfo; i++) {
		this->mapXMLCh2AttType[xns::XMLString::transcode(attributeInfo[i].key)] = attributeInfo[i].type;
	}
	this->totals.reserve(50);
	this->acts.reserve(50);
}

SaxPhreeqcHandlers::~SaxPhreeqcHandlers()
{
	std::map<const XMLCh*, elementType, XMLCH_LESS>::iterator it = this->mapXMLCh2Type.begin();
	for (; it != this->mapXMLCh2Type.end(); ++it)
	{
		XMLCh* p = (XMLCh*)it->first;
		xns::XMLString::release(&p);
	}

}

// -----------------------------------------------------------------------
//  Implementations of the SAX DocumentHandler interface
// -----------------------------------------------------------------------

void SaxPhreeqcHandlers::endDocument()
{
}

void SaxPhreeqcHandlers::endElement(const XMLCh* const name)
{
  switch (this->mapXMLCh2Type[name])
    {
	case typeSOLUTION:
		// solution is finished now copy into solutions array
		{
			int n;
			// copy vector of conc's to solution
			this->solution_ptr->totals = (struct conc*) PHRQ_realloc(this->solution_ptr->totals, (size_t) (this->totals.size() + 1) * sizeof(struct conc));
			std::copy(this->totals.begin(), this->totals.end(), this->solution_ptr->totals);
			this->solution_ptr->totals[this->totals.size()].description=NULL;
			this->totals.clear();
			assert(this->totals.size() == 0);

			// copy vector of master_activities's to solution
			this->solution_ptr->master_activity = (struct master_activity *) free_check_null(this->solution_ptr->master_activity);
			if (this->acts.size() > 0) {
				this->solution_ptr->master_activity = (struct master_activity *) PHRQ_realloc(this->solution_ptr->master_activity, (size_t) (this->acts.size()) * sizeof(struct master_activity));
				std::copy(this->acts.begin(), this->acts.end(), this->solution_ptr->master_activity);
			}
			this->solution_ptr->count_master_activity = this->acts.size();
			this->acts.clear();
			assert(this->acts.size() == 0);

			// copy vector of s_gamma's to solution
			this->solution_ptr->species_gamma = (struct master_activity *) free_check_null(this->solution_ptr->species_gamma);
			if (this->s_gammas.size() > 0) {
				this->solution_ptr->species_gamma = (struct master_activity *) PHRQ_realloc(this->solution_ptr->species_gamma, (size_t) (this->s_gammas.size()) * sizeof(struct master_activity));
				std::copy(this->s_gammas.begin(), this->s_gammas.end(), this->solution_ptr->species_gamma);
			}
			this->solution_ptr->count_species_gamma = this->s_gammas.size();
			this->s_gammas.clear();
			assert(this->s_gammas.size() == 0);

			// copy vector of isotopes's to solution
			this->solution_ptr->isotopes = (struct isotope *) free_check_null(this->solution_ptr->isotopes);
			if (this->isotopes.size() > 0) {
				this->solution_ptr->isotopes = (struct isotope *) PHRQ_realloc(this->solution_ptr->isotopes, (size_t) (this->isotopes.size()) * sizeof(struct isotope));
				std::copy(this->isotopes.begin(), this->isotopes.end(), this->solution_ptr->isotopes);
			}
			this->solution_ptr->count_isotopes = this->isotopes.size();
			this->isotopes.clear();
			assert(this->isotopes.size() == 0);


			// store solution for now
			if (solution_bsearch(this->solution_ptr->n_user, &n, FALSE) != NULL) {
				solution_free(solution[n]);
				solution[n] = this->solution_ptr;
			} else {
				n = count_solution++;
				if (count_solution >= max_solution)
				{
					space ((void **) &(solution), count_solution, &max_solution, sizeof (struct solution *) );
				}
				solution[n] = this->solution_ptr;
			}
			this->solution_ptr = NULL;
		}
		break;
 /*     
    case typeTOTALS:
      // have all totals now allocate and copy into solution_ptr
      m_solution_ptr->totals = (struct conc*) PHRQ_realloc(m_solution_ptr->totals, (size_t) (m_totals.size() + 1) * sizeof(struct conc));
      std::copy(m_totals.begin(), m_totals.end(), m_solution_ptr->totals);
      m_totals.clear();
      assert(m_totals.size() == 0);
      break;
      
    case typeACTS:
      // have all acts now allocate and copy into solution_ptr
      m_solution_ptr->master_activity = (struct master_activity*) PHRQ_realloc(m_solution_ptr->master_activity, (size_t) (m_acts.size() + 1) * sizeof(struct master_activity));
      std::copy(m_acts.begin(), m_acts.end(), m_solution_ptr->master_activity);
      m_acts.clear();
      assert(m_acts.size() == 0);
      break;
*/      
    default:
      break;
    }
}

void SaxPhreeqcHandlers::characters(const XMLCh* const chars, const unsigned int length)
{
  // skip whitespace
  XMLCh* pChar = (XMLCh*)chars;
/*
  while(pChar && iswspace(*pChar)) ++pChar;
  if (*pChar)
    {
      switch(m_type)
	{

	  case typeTC:
	  m_solution_ptr->tc = wcstod((wchar_t*)pChar, NULL);
	  break;
	case typePH:
	  m_solution_ptr->ph = wcstod((wchar_t*)pChar, NULL);
	  break;
	case typeSOLUTION_PE:
	  m_solution_ptr->solution_pe = wcstod((wchar_t*)pChar, NULL);
	  break;
	case typeMU:
	  m_solution_ptr->mu = wcstod((wchar_t*)pChar, NULL);
	  break;
	case typeAH2O:
	  m_solution_ptr->ah2o = wcstod((wchar_t*)pChar, NULL);
	  break;
	case typeTOTAL_H:
	  m_solution_ptr->total_h = wcstod((wchar_t*)pChar, NULL);
	  break;
	case typeTOTAL_O:
	  m_solution_ptr->total_o = wcstod((wchar_t*)pChar, NULL);
	  break;
	case typeCB:
	  m_solution_ptr->cb = wcstod((wchar_t*)pChar, NULL);
	  break;
	case typeMASS_WATER:
	  m_solution_ptr->mass_water = wcstod((wchar_t*)pChar, NULL);
	  break;
	case typeTOTAL:
	  // description stored in startElement
	  m_solution_ptr->totals[0].moles = wcstod((wchar_t*)pChar, NULL);
	  m_totals.push_back(m_solution_ptr->totals[0]);
	  break;
	case typeACT:
	  // description stored in startElement
	  m_solution_ptr->master_activity[0].la = wcstod((wchar_t*)pChar, NULL);
	  m_acts.push_back(m_solution_ptr->master_activity[0]);
	  break;

	default:
	  break;
	}
 	}
	  */
}

void SaxPhreeqcHandlers::ignorableWhitespace(const XMLCh* const chars, const unsigned int length)
{
}

void SaxPhreeqcHandlers::processingInstruction(const XMLCh* const target, const XMLCh* const data)
{
}

void SaxPhreeqcHandlers::startDocument()
{
}
void SaxPhreeqcHandlers::startElement(const XMLCh* const name, xns::AttributeList& attributes)
{
  const char ERR_MSG[] = "Unpacking solution message: %s, element not found\n";
  char *string;
  int i;

  string = xns::XMLString::transcode(name);
  this->eltType = this->mapXMLCh2Type[name];
  xns::XMLString::release(&string);
  switch (this->eltType)
    {
    case typePHAST_STATE:
		i =  wcstol((wchar_t*)attributes.getValue(xns::XMLString::transcode("nx")), NULL, 10);

	// m_solution_ptr->n_user = (int) wcstol((wchar_t*)attributes.getValue(s_strN_USER), NULL, 10);
      break;
    case typeSOLUTION:
      assert(this->solution_ptr == NULL);
      assert(this->totals.size() == 0);
      assert(this->acts.size() == 0);
      
      // allocate space for solution
      this->solution_ptr = solution_alloc();
      
      // process attributes for solution
	  processSolutionAttributes(attributes);
	  break;
    case typeSOLN_PE:
		assert(this->solution_ptr->pe != NULL);
		// store pe, no need to clean up at end of solution
		string = xns::XMLString::transcode(attributes.getName(0));
		if ((attributes.getLength() >= 1) && (this->mapXMLCh2AttType[attributes.getName(0)] == attSOLN_PE_name)){
			string = xns::XMLString::transcode(attributes.getValue((unsigned int) 0));
			pe_data_store(&(this->solution_ptr->pe), string);
		} else {
			++input_error;
			sprintf(error_string, "No attribute data for SOLN_PE.\n");
			error_msg(error_string, CONTINUE);
		}
		break;
    case typeSOLN_TOTAL:
		{
			// store in c, push_back on totals
			// need to copy and clean up at end of </solution>
			struct conc *c;
			c = processSolutionTotalAttributes(attributes);
			this->totals.push_back(*c);
		}
		break;  
    case typeSOLN_MASTER_ACTIVITY:
		{
			// store in ma, push_back on acts
			// need to copy and clean up at end of </solution>
			struct master_activity *ma;
			ma = processMasterActivityAttributes(attributes);
			this->acts.push_back(*ma);
		}
    case typeSOLN_SPECIES_GAMMA:
		{
			// store in ma, push_back on s_gammas
			// need to copy and clean up at end of </solution>
			struct master_activity *ma;
			ma = processMasterActivityAttributes(attributes);
			this->s_gammas.push_back(*ma);
			break;
		}
    case typeSOLN_ISOTOPE:
		{
			// store in iso, push_back on isotopes
			// need to copy and clean up at end of </solution>
			struct isotope *iso;
			iso = processIsotopeAttributes(attributes);
			this->isotopes.push_back(*iso);
			break;
		}
		break;  
				
				
		/*     
    case typeTOTAL:
      {
	m_solution_ptr->totals[0].description = XMLString::transcode(attributes.getValue(s_strNAME));
	struct master *master_ptr = master_bsearch(m_solution_ptr->totals[0].description);
	if (master_ptr == NULL)
	  {
	    ++input_error;
	    sprintf(error_string, ERR_MSG, m_solution_ptr->totals[0].description);
	    error_msg(error_string, CONTINUE);
	    delete[] m_solution_ptr->totals[0].description;
	  }
	else
	  {
	    delete[] m_solution_ptr->totals[0].description;
	    m_solution_ptr->totals[0].description = master_ptr->elt->name;
	  }
      }
      break;
      
    case typeACT:
      {
	m_solution_ptr->master_activity[0].description = XMLString::transcode(attributes.getValue(s_strNAME));
	struct master *master_ptr = master_bsearch(m_solution_ptr->master_activity[0].description);
	if (master_ptr == NULL)
	  {
	    ++input_error;
	    sprintf(error_string, ERR_MSG, m_solution_ptr->master_activity[0].description);
	    error_msg(error_string, CONTINUE);
	    delete[] m_solution_ptr->master_activity[0].description;
	  }
	else
	  {
	    delete[] m_solution_ptr->master_activity[0].description;
	    m_solution_ptr->master_activity[0].description = master_ptr->elt->name;
	  }
      }
      break;
      */
    default:
      break;
    }
}
int SaxPhreeqcHandlers::processSolutionAttributes(xns::AttributeList& attributes)
{
	const char ERR_MSG[] = "Unpacking solution attributes: %s, attribute not found\n";	
	unsigned int i;
	char *string;
	attributeType attType;
	assert(this->eltType == typeSOLUTION);
	assert(this->solution_ptr != NULL);

	// Get attribute name, map to attribute type, process

	for (i = 0; i < attributes.getLength(); i++) {
		attType = this->mapXMLCh2AttType[attributes.getName(i)];
		switch (attType) {
		case attSOLN_new_def:
			this->solution_ptr->new_def = strtol(string = xns::XMLString::transcode(attributes.getValue(i)), NULL, 10);
			xns::XMLString::release(&string);
			break;
		case attSOLN_n_user:
			this->solution_ptr->n_user = strtol(xns::XMLString::transcode(attributes.getValue(i)), NULL, 10);
			break;
		case attSOLN_n_user_end:
			this->solution_ptr->n_user_end = strtol(xns::XMLString::transcode(attributes.getValue(i)), NULL, 10);
			break;
		case attSOLN_description:
			this->solution_ptr->description = string_duplicate(xns::XMLString::transcode(attributes.getValue(i)));
			break;
		case attSOLN_tc:
			this->solution_ptr->tc = strtod(xns::XMLString::transcode(attributes.getValue(i)), NULL);
			break;
		case attSOLN_ph:
			this->solution_ptr->ph = strtod(xns::XMLString::transcode(attributes.getValue(i)), NULL);
			break;
		case attSOLN_solution_pe:
			this->solution_ptr->solution_pe = strtod(xns::XMLString::transcode(attributes.getValue(i)), NULL);
			break;
		case attSOLN_mu:
			this->solution_ptr->mu = strtod(xns::XMLString::transcode(attributes.getValue(i)), NULL);
			break;
		case attSOLN_ah2o:
			this->solution_ptr->ah2o = strtod(xns::XMLString::transcode(attributes.getValue(i)), NULL);
			break;
		case attSOLN_density:
			this->solution_ptr->density = strtod(xns::XMLString::transcode(attributes.getValue(i)), NULL);
			break;
		case attSOLN_total_h:
			this->solution_ptr->total_h = strtod(xns::XMLString::transcode(attributes.getValue(i)), NULL);
			break;
		case attSOLN_total_o:
			this->solution_ptr->total_o = strtod(xns::XMLString::transcode(attributes.getValue(i)), NULL);
			break;
		case attSOLN_cb:
			this->solution_ptr->cb = strtod(xns::XMLString::transcode(attributes.getValue(i)), NULL);
			break;
		case attSOLN_mass_water:
			this->solution_ptr->mass_water = strtod(xns::XMLString::transcode(attributes.getValue(i)), NULL);
			break;
		case attSOLN_total_alkalinity:
			this->solution_ptr->total_alkalinity = strtod(xns::XMLString::transcode(attributes.getValue(i)), NULL);
			break;
		case attSOLN_total_co2:
			this->solution_ptr->total_co2 = strtod(xns::XMLString::transcode(attributes.getValue(i)), NULL);
			break;
		case attSOLN_units:
			this->solution_ptr->units = string_hsave(xns::XMLString::transcode(attributes.getValue(i)));
			break;
		case attSOLN_default_pe:
			this->solution_ptr->default_pe = strtol(xns::XMLString::transcode(attributes.getValue(i)), NULL, 10);
			break;
			break;
		case attSOLN_count_master_activity:
			this->solution_ptr->count_master_activity = strtol(xns::XMLString::transcode(attributes.getValue(i)), NULL, 10);
			break;
		case attSOLN_count_isotopes:
			this->solution_ptr->count_isotopes = strtol(xns::XMLString::transcode(attributes.getValue(i)), NULL, 10);
			break;
		case attSOLN_count_species_gamma:
			this->solution_ptr->count_species_gamma = strtol(xns::XMLString::transcode(attributes.getValue(i)), NULL, 10);
			break;
		default:
			++input_error;
			sprintf(error_string, ERR_MSG, xns::XMLString::transcode(attributes.getName(i)));
			error_msg(error_string, CONTINUE);		
			break;
		}
	}
	return 0;
}
struct conc *SaxPhreeqcHandlers::processSolutionTotalAttributes(xns::AttributeList& attributes)
{
	const char ERR_MSG[] = "Unpacking solution totals attributes: %s, attribute not found\n";	
	unsigned int i;
	int l;
	struct conc *c = new conc();
	attributeType attType;
	assert(this->eltType == typeSOLN_TOTAL);


	// Get attribute name, map to attribute type, process

	for (i = 0; i < attributes.getLength(); i++) {
		attType = this->mapXMLCh2AttType[attributes.getName(i)];
		switch (attType) {
		case attCONC_description:
			c->description = string_hsave(xns::XMLString::transcode(attributes.getValue(i)));
			break;
		case attCONC_skip:
			c->skip = strtol(xns::XMLString::transcode(attributes.getValue(i)), NULL, 10);
			break;
		case attCONC_moles:
			c->moles = strtod(xns::XMLString::transcode(attributes.getValue(i)), NULL);
			break;
		case attCONC_input_conc:
			c->input_conc = strtod(xns::XMLString::transcode(attributes.getValue(i)), NULL);
			break;		
		case attCONC_equation_name:
			c->equation_name = string_hsave(xns::XMLString::transcode(attributes.getValue(i)));
			c->phase = phase_bsearch(c->equation_name, &l, FALSE);
			break;
		case attCONC_phase_si:
			c->phase_si = strtod(xns::XMLString::transcode(attributes.getValue(i)), NULL);
			break;		
		case attCONC_n_pe:
			c->n_pe = strtol(xns::XMLString::transcode(attributes.getValue(i)), NULL, 10);
			break;
		case attCONC_as:
			c->as = string_hsave(xns::XMLString::transcode(attributes.getValue(i)));
			break;
		case attCONC_gfw:
			c->gfw = strtod(xns::XMLString::transcode(attributes.getValue(i)), NULL);
			break;		
	
		default:
			++input_error;
			sprintf(error_string, ERR_MSG, xns::XMLString::transcode(attributes.getName(i)));
			error_msg(error_string, CONTINUE);		
			break;
		}
	}
	return c;
}
struct master_activity *SaxPhreeqcHandlers::processMasterActivityAttributes(xns::AttributeList& attributes)
{
	int i;
	const char ERR_MSG[] = "Unpacking master activity attributes: %s, attribute not found\n";	
	struct master_activity *ma = new master_activity();
	attributeType attType;
	for (i = 0; i < attributes.getLength(); i++) {
		attType = this->mapXMLCh2AttType[attributes.getName(i)];
		switch (attType) {
			case attM_A_description:
				ma->description = string_hsave(xns::XMLString::transcode(attributes.getValue(i)));
				break;
			case attM_A_la:
				ma->la = strtod(xns::XMLString::transcode(attributes.getValue(i)), NULL);
				break;
			default:
				++input_error;
				sprintf(error_string, ERR_MSG, xns::XMLString::transcode(attributes.getName(i)));
				error_msg(error_string, CONTINUE);		
				break;
		}	
	}
	return ma;
}
struct isotope *SaxPhreeqcHandlers::processIsotopeAttributes(xns::AttributeList& attributes)
{
	int i;
	const char ERR_MSG[] = "Unpacking isotope attributes: %s, attribute not found\n";	
	struct isotope *iso = new isotope();

	iso->primary = iso->master = NULL;
	attributeType attType;
	for (i = 0; i < attributes.getLength(); i++) {
		attType = this->mapXMLCh2AttType[attributes.getName(i)];
		switch (attType) {
			case attISO_isotope_number:
				iso->isotope_number = strtol(xns::XMLString::transcode(attributes.getValue(i)), NULL, 10);
				break;
			case attISO_elt_name:
				iso->elt_name = string_hsave(xns::XMLString::transcode(attributes.getValue(i)));
				break;
			case attISO_isotope_name:
				iso->isotope_name = string_hsave(xns::XMLString::transcode(attributes.getValue(i)));
				break;					
			case attISO_total:
				iso->total = strtod(xns::XMLString::transcode(attributes.getValue(i)), NULL);
				break;
			case attISO_ratio:
				iso->ratio = strtod(xns::XMLString::transcode(attributes.getValue(i)), NULL);
				break;
			case attISO_ratio_uncertainty:
				iso->ratio_uncertainty = strtod(xns::XMLString::transcode(attributes.getValue(i)), NULL);
				break;
			case attISO_coef:
				iso->coef = strtod(xns::XMLString::transcode(attributes.getValue(i)), NULL);
				break;

			default:
				++input_error;
				sprintf(error_string, ERR_MSG, xns::XMLString::transcode(attributes.getName(i)));
				error_msg(error_string, CONTINUE);		
				break;
		}	
	}
	return iso;
}