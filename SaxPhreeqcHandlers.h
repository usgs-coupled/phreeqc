// SaxPhreeqcHandlers.h: interface for the SaxPhreeqcHandlers class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_SAXPHREEQCHANDLERS_H__4A69D5F5_2E57_4001_911D_4ABF6F2C0A0B__INCLUDED_)
#define AFX_SAXPHREEQCHANDLERS_H__4A69D5F5_2E57_4001_911D_4ABF6F2C0A0B__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <wchar.h>                         // iswspace sprintf
#include <xercesc/sax/HandlerBase.hpp>
/**
#include <hash_map>
**/
#include <vector>
#include <map>

#define xns XERCES_CPP_NAMESPACE

#ifdef USE_LONG_DOUBLE
#define LDBLE long double
#else
#define LDBLE double
#endif

struct conc {
	char *description;
	int skip;
	double moles;
	double input_conc;
	char *units;
	char *equation_name;
	struct phase *phase;
	double phase_si;
	int n_pe;
	char *as;
	double gfw;
};
struct isotope {
	LDBLE isotope_number;
	char *elt_name;
	char *isotope_name;
	LDBLE total;
	LDBLE ratio;
	LDBLE ratio_uncertainty;
	LDBLE x_ratio_uncertainty;
	struct master *master;
	struct master *primary;
	LDBLE coef;                    /* coefficient of element in phase */
};
struct master_activity {
	char *description;
	LDBLE la;
};

///XERCES_CPP_NAMESPACE_USE

struct XMLCH_LESS : std::binary_function<const XMLCh*, const XMLCh*, bool> {
bool operator()(const XMLCh* _X, const XMLCh* _Y) const
{return xns::XMLString::compareString( _X, _Y) < 0;}
//	{return wcscmp((const wchar_t*) _X, (const wchar_t*) _Y) < 0; }

};

/**
struct XMLCH_EQUALS : std::binary_function<const XMLCh*, const XMLCh*, bool> {
bool operator()(const XMLCh* _X, const XMLCh* _Y) const
	{return wcscmp((const wchar_t*) _X, (const wchar_t*) _Y) == 0; }
};
**/


class SaxPhreeqcHandlers : public xns::HandlerBase
{
public:
	SaxPhreeqcHandlers();
	virtual ~SaxPhreeqcHandlers();

    // -----------------------------------------------------------------------
    //  Implementations of the SAX DocumentHandler interface
    // -----------------------------------------------------------------------
    void endDocument();

    void endElement(const XMLCh* const name);

    void characters(const XMLCh* const chars, const unsigned int length);

    void ignorableWhitespace(const XMLCh* const chars, const unsigned int length);

    void processingInstruction(const XMLCh* const target, const XMLCh* const data);

    void startDocument();

    void startElement(const XMLCh* const name, xns::AttributeList& attributes);

	// element types
	enum elementType
	{
		typeNULL,
		typePHAST_STATE,
		typeSYSTEM,
		typeSOLUTION,
		typeSOLN_PE,
		typeSOLN_TOTAL,
		typeSOLN_MASTER_ACTIVITY,
		typeSOLN_ISOTOPE,
		typeSOLN_SPECIES_GAMMA
	} eltType;
	enum attributeType
	{
		attNULL,
		attSOLN_new_def,
		attSOLN_n_user,
		attSOLN_n_user_end,
		attSOLN_description,
		attSOLN_tc,
		attSOLN_ph,
		attSOLN_solution_pe,
		attSOLN_mu,
		attSOLN_ah2o,
		attSOLN_density,
		attSOLN_total_h,
		attSOLN_total_o,
		attSOLN_cb,
		attSOLN_mass_water,
		attSOLN_total_alkalinity,
		attSOLN_total_co2,
		attSOLN_units,
		attSOLN_count_pe,
		attSOLN_default_pe,
		attSOLN_count_totals,
		attSOLN_count_master_activity,
		attSOLN_count_isotopes,
		attSOLN_count_species_gamma,
		attSOLN_PE_name,
		attM_A_description,
		attM_A_la,
		attISO_isotope_number,
		attISO_elt_name,
		attISO_isotope_name,
		attISO_total,
		attISO_ratio,		
		attISO_ratio_uncertainty,			
		attISO_x_ratio_uncertainty,			
		attISO_coef,	
		attCONC_description,
		attCONC_skip,
		attCONC_moles,
		attCONC_input_conc,
		attCONC_equation_name,
		attCONC_phase_si,
		attCONC_n_pe,
		attCONC_as,
		attCONC_gfw,
	} attType;

	int  processSolutionAttributes(xns::AttributeList& attributes);
	struct conc *processSolutionTotalAttributes(xns::AttributeList& attributes);
	struct master_activity *processMasterActivityAttributes(xns::AttributeList& attributes);
	struct isotope *processIsotopeAttributes(xns::AttributeList& attributes);

protected:
	std::vector<conc> totals;
	std::vector<master_activity> acts, s_gammas;
	std::vector<isotope> isotopes;
	std::map<const XMLCh*, elementType, XMLCH_LESS> mapXMLCh2Type;
	std::map<const XMLCh*, attributeType, XMLCH_LESS> mapXMLCh2AttType;
	/**
	std::hash_map<const XMLCh*, ElementType, hash<const XMLCh*>, XMLCH_EQUALS> m_hashmapXMLCh2Type;
	**/
	struct solution* solution_ptr;
};

#endif // !defined(AFX_SAXPHREEQCHANDLERS_H__4A69D5F5_2E57_4001_911D_4ABF6F2C0A0B__INCLUDED_)
