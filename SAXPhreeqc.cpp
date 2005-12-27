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

#include <xercesc/parsers/SAXParser.hpp>           // SAXParser
#include <xercesc/sax/AttributeList.hpp>           // AttributeList
#include <xercesc/util/XMLUniDefs.hpp>             // Unicode definitions
#include <xercesc/framework/MemBufInputSource.hpp> // MemBufInputSource
#include <xercesc/util/PlatformUtils.hpp>          // XMLPlatformUtils::getCurrentMillis

#include "SAXPhreeqc.h"                    // SAX_ functions
#include "SaxPhreeqcHandlers.h"            // SaxPhreeqcHandlers

///static char buffer[300000];
///static std::ostrstream s_oss(buffer, 300000);        // must never go out of scope
static std::ostrstream s_oss;        // must never go out of scope
static bool s_bSysIsOpen = false;    // must never go out of scope
static SaxPhreeqcHandlers s_handler; // one and only instance

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
  int count_isotopes;
  struct isotope *isotopes;
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


// extern routines
extern "C" {
#include "phqalloc.h"
  struct solution *solution_alloc(void);
  struct solution *solution_bsearch(int k, int *n, int print);
  int solution_free (struct solution *solution_ptr);
  void space (void **ptr, int i, int *max, int struct_size);
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
  Initializer(){ XMLPlatformUtils::Initialize();}
};

static Initializer sInit;  // initialize xerces once

extern "C" void SAX_StartSystem()
{
  assert(!s_bSysIsOpen);     // system already open and has not been closed
  
  // init stream
  s_oss.freeze(false);
  s_oss.seekp(0);
  
  // write stream
  s_oss << "<system>";
  s_bSysIsOpen = true;
}

extern "C" int SAX_AddSolution(struct solution* solution_ptr)
{
  const char    ERR_MESSAGE[] = "Packing solution message: %s, element not found\n";
  
  assert(s_bSysIsOpen);      // must call SAX_StartSystem first

  // <solution num="1">
  s_oss << "<solution num=\"" << solution_ptr->n_user << "\">";
  
  // set output precision
  s_oss.precision(DBL_DIG + 1);
  //s_oss << std::ios::scientific;
  
  // <temp_c>25.0</temp_c>
  s_oss << "<temp_c>" << solution_ptr->tc << "</temp_c>";
  
  // <pH>7.0</pH>
  s_oss << "<pH>" << solution_ptr->ph << "</pH>";
  
  // <solution_pe>4.0</solution_pe>
  s_oss << "<solution_pe>" << solution_ptr->solution_pe << "</solution_pe>";
  
  // <mu>1e-007</mu>
  s_oss << "<mu>" << solution_ptr->mu << "</mu>";
  
  // <ah2o>1</ah2o>
  s_oss << "<ah2o>" << solution_ptr->ah2o << "</ah2o>";
  
  // <total_h>111.0147</total_h>
  s_oss << "<total_h>" << solution_ptr->total_h << "</total_h>";
  
  // <total_o>55.63047</total_o>
  s_oss << "<total_o>" << solution_ptr->total_o << "</total_o>";
  
  // <cb>0</cb>
  s_oss << "<cb>" << solution_ptr->cb << "</cb>";
  
  // <mass_water>1</mass_water>
  s_oss << "<mass_water>" << solution_ptr->mass_water << "</mass_water>";
  
  // <totals>
  s_oss << "<totals>";
  
  struct master *master_ptr;
  for (int i = 0; solution_ptr->totals[i].description != NULL; ++i)
    {
      master_ptr = master_bsearch(solution_ptr->totals[i].description);
      if (master_ptr == NULL)
	{
	  ++input_error;
	  sprintf(error_string, ERR_MESSAGE, solution_ptr->totals[i].description);
	  error_msg(error_string, CONTINUE);
	}
      
      // <total name="Alkalinity">2.406e-003</total>
      s_oss << "<total name=\"" << solution_ptr->totals[i].description << "\">";
      s_oss << solution_ptr->totals[i].moles;
      s_oss << "</total>";
    }
  
  // </totals>
  s_oss << "</totals>";
  
  // <acts>
  s_oss << "<acts>";
  
  for (int j = 0; solution_ptr->master_activity[j].description != NULL; ++j)
    {
      master_ptr = master_bsearch(solution_ptr->master_activity[j].description);
      if (master_ptr == NULL)
	{
	  ++input_error;
	  sprintf(error_string, ERR_MESSAGE, solution_ptr->master_activity[j].description);
	  error_msg(error_string, CONTINUE);
	}
      
      // <act name="OH-">-5.788</act>
      s_oss << "<act name=\"" << solution_ptr->master_activity[j].description << "\">";
      s_oss << solution_ptr->master_activity[j].la;
      s_oss << "</act>";
    }
  
  // </acts>
  s_oss << "</acts>";
  
  // </solution>
  s_oss << "</solution>";
  
  return(OK);
}

extern "C" void SAX_EndSystem()
{
  assert(s_bSysIsOpen);      // must call SAX_StartSystem first

  s_oss << "</system>";
  
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
  MemBufInputSource memBufIS((const XMLByte*)pvBuffer, buf_size, "solution_id", false);

  //
  //  Create a SAX parser object.
  //
  SAXParser parser;
  parser.setValidationScheme(SAXParser::Val_Always);
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
      parser.setDocumentHandler(&s_handler);
      parser.setErrorHandler(&s_handler);
      
      int t = 0;
      for (; t < 10; ++t)
	{
	  const unsigned long startMillis = XMLPlatformUtils::getCurrentMillis();
	  parser.parse(memBufIS);
	  const unsigned long endMillis = XMLPlatformUtils::getCurrentMillis();
	  duration += endMillis - startMillis;
	}
      std::cerr << "\nSaxParse time = " << (float)duration/t << " millis\n";
    }
  catch (const SAXException& toCatch)
    {
      char* psz = XMLString::transcode(toCatch.getMessage());
      input_error++;
      sprintf(error_string, "SAX_UnpackSolutions: %s\n", psz);
      error_msg(error_string, CONTINUE);
      delete[] psz;
      return ERROR;
    }
  catch (const XMLException& toCatch)
    {
      char* psz = XMLString::transcode(toCatch.getMessage());
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
const XMLCh SaxPhreeqcHandlers::s_strSYSTEM[] 	    =  // "system"
  {chLatin_s, chLatin_y, chLatin_s, chLatin_t, chLatin_e, chLatin_m, 0};

const XMLCh SaxPhreeqcHandlers::s_strSOLUTION[] 	=  // "solution"
  {chLatin_s, chLatin_o, chLatin_l, chLatin_u, chLatin_t, chLatin_i, chLatin_o, chLatin_n, 0};

const XMLCh SaxPhreeqcHandlers::s_strTC[] 			=  // "temp_c"
  {chLatin_t, chLatin_e, chLatin_m, chLatin_p, chUnderscore, chLatin_c, 0};

const XMLCh SaxPhreeqcHandlers::s_strPH[] 			=  // "pH"
  {chLatin_p, chLatin_H, 0};

const XMLCh SaxPhreeqcHandlers::s_strSOLUTION_PE[] 	=  // "solution_pe"
  {chLatin_s, chLatin_o, chLatin_l, chLatin_u, chLatin_t, chLatin_i, chLatin_o, chLatin_n, chUnderscore, chLatin_p, chLatin_e, 0};

const XMLCh SaxPhreeqcHandlers::s_strMU[] 			=  // "mu"
  {chLatin_m, chLatin_u, 0};

const XMLCh SaxPhreeqcHandlers::s_strAH2O[] 		=  // "ah2o"
  {chLatin_a, chLatin_h, chDigit_2, chLatin_o, 0};

const XMLCh SaxPhreeqcHandlers::s_strTOTAL_H[] 		=  // "total_h"
  {chLatin_t, chLatin_o, chLatin_t, chLatin_a, chLatin_l, chUnderscore, chLatin_h, 0};

const XMLCh SaxPhreeqcHandlers::s_strTOTAL_O[] 		=  // "total_o"
  {chLatin_t, chLatin_o, chLatin_t, chLatin_a, chLatin_l, chUnderscore, chLatin_o, 0};

const XMLCh SaxPhreeqcHandlers::s_strCB[] 			=  // "cb"
  {chLatin_c, chLatin_b, 0};

const XMLCh SaxPhreeqcHandlers::s_strMASS_WATER[] 	=  // "mass_water"
  {chLatin_m, chLatin_a, chLatin_s, chLatin_s, chUnderscore, chLatin_w, chLatin_a, chLatin_t, chLatin_e, chLatin_r, 0};

const XMLCh SaxPhreeqcHandlers::s_strTOTALS[] 		=  // "totals"
  {chLatin_t, chLatin_o, chLatin_t, chLatin_a, chLatin_l, chLatin_s, 0};

const XMLCh SaxPhreeqcHandlers::s_strTOTAL[] 		=  // "total"
  {chLatin_t, chLatin_o, chLatin_t, chLatin_a, chLatin_l, 0};

const XMLCh SaxPhreeqcHandlers::s_strACTS[] 		=  // "acts"
  {chLatin_a, chLatin_c, chLatin_t, chLatin_s, 0};

const XMLCh SaxPhreeqcHandlers::s_strACT[]  		=  // "act"
  {chLatin_a, chLatin_c, chLatin_t, 0};


// attribute names
const XMLCh SaxPhreeqcHandlers::s_strN_USER[] 		=  // "num"
  {chLatin_n, chLatin_u, chLatin_m, 0};

const XMLCh SaxPhreeqcHandlers::s_strNAME[] 		=  // "name"
  {chLatin_n, chLatin_a, chLatin_m, chLatin_e, 0};




SaxPhreeqcHandlers::SaxPhreeqcHandlers()
  : m_type(typeNULL), m_totals(), m_acts(), m_solution_ptr(NULL)
{
  m_mapXMLCh2Type[s_strSYSTEM]      = typeSYSTEM;
  m_mapXMLCh2Type[s_strSOLUTION]    = typeSOLUTION;
  m_mapXMLCh2Type[s_strTC]          = typeTC;
  m_mapXMLCh2Type[s_strPH]          = typePH;
  m_mapXMLCh2Type[s_strSOLUTION_PE] = typeSOLUTION_PE;
  m_mapXMLCh2Type[s_strMU]          = typeMU;
  m_mapXMLCh2Type[s_strAH2O]        = typeAH2O;
  m_mapXMLCh2Type[s_strTOTAL_H]     = typeTOTAL_H;
  m_mapXMLCh2Type[s_strTOTAL_O]     = typeTOTAL_O;
  m_mapXMLCh2Type[s_strCB]          = typeCB;
  m_mapXMLCh2Type[s_strMASS_WATER]  = typeMASS_WATER;
  m_mapXMLCh2Type[s_strTOTALS]      = typeTOTALS;
  m_mapXMLCh2Type[s_strTOTAL]       = typeTOTAL;
  m_mapXMLCh2Type[s_strACTS]        = typeACTS;
  m_mapXMLCh2Type[s_strACT]         = typeACT;
  
  m_totals.reserve(50);
  m_acts.reserve(50);
}

SaxPhreeqcHandlers::~SaxPhreeqcHandlers()
{
}

// -----------------------------------------------------------------------
//  Implementations of the SAX DocumentHandler interface
// -----------------------------------------------------------------------

void SaxPhreeqcHandlers::endDocument()
{
}

void SaxPhreeqcHandlers::endElement(const XMLCh* const name)
{
  switch (m_mapXMLCh2Type[name])
    {
    case typeSOLUTION:
      // solution is finished now copy into solutions array
      {
	int n;
	if (solution_bsearch(m_solution_ptr->n_user, &n, FALSE) != NULL)
	  {
	    solution_free(solution[n]);
	    solution[n] = m_solution_ptr;
	  }
	else
	  {
	    n = count_solution++;
	    if (count_solution >= max_solution)
	      {
		space ((void **) &(solution), count_solution, &max_solution, sizeof (struct solution *) );
	      }
	    solution[n] = m_solution_ptr;
	  }
	m_solution_ptr = NULL;
      }
      break;
      
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
      
    default:
      break;
    }
}

void SaxPhreeqcHandlers::characters(const XMLCh* const chars, const unsigned int length)
{
  // skip whitespace
  XMLCh* pChar = (XMLCh*)chars;
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

void SaxPhreeqcHandlers::startElement(const XMLCh* const name, AttributeList& attributes)
{
  const char ERR_MSG[] = "Unpacking solution message: %s, element not found\n";
  m_type = m_mapXMLCh2Type[name];
  
  switch (m_type)
    {
    case typeSOLUTION:
      assert(m_solution_ptr == NULL);
      assert(m_totals.size() == 0);
      assert(m_acts.size() == 0);
      
      // allocate space for solution
      m_solution_ptr = solution_alloc();
      
      // get n_user
      m_solution_ptr->n_user = (int) wcstol((wchar_t*)attributes.getValue(s_strN_USER), NULL, 10);
      break;
      
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
      
    default:
      break;
    }
}
