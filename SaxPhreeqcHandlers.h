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
struct master_activity {
	char *description;
	LDBLE la;
};

struct XMLCH_LESS : std::binary_function<const XMLCh*, const XMLCh*, bool> {
bool operator()(const XMLCh* _X, const XMLCh* _Y) const
	{return wcscmp((const wchar_t*) _X, (const wchar_t*) _Y) < 0; }
};

/**
struct XMLCH_EQUALS : std::binary_function<const XMLCh*, const XMLCh*, bool> {
bool operator()(const XMLCh* _X, const XMLCh* _Y) const
	{return wcscmp((const wchar_t*) _X, (const wchar_t*) _Y) == 0; }
};
**/
XERCES_CPP_NAMESPACE_USE

class SaxPhreeqcHandlers : public HandlerBase
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

    void startElement(const XMLCh* const name, AttributeList& attributes);

	// element names
	static const XMLCh s_strSYSTEM[];
	static const XMLCh s_strSOLUTION[];
	static const XMLCh s_strTC[];
	static const XMLCh s_strPH[];
	static const XMLCh s_strSOLUTION_PE[];
	static const XMLCh s_strMU[];
	static const XMLCh s_strAH2O[];
	static const XMLCh s_strTOTAL_H[];
	static const XMLCh s_strTOTAL_O[];
	static const XMLCh s_strCB[];
	static const XMLCh s_strMASS_WATER[];
	static const XMLCh s_strTOTALS[];
	static const XMLCh s_strTOTAL[];
	static const XMLCh s_strACTS[];
	static const XMLCh s_strACT[];

	// attribute names
	static const XMLCh s_strN_USER[];
	static const XMLCh s_strNAME[];


	// element types
	enum ElementType
	{
		typeNULL,
		typePHAST_STATE,
		typeSYSTEM,
		typeSOLUTION,
		typeSOLN_PE,
		typeREACTION,
		typeSOLN_TOTAL,
		typeSOLN_MASTER_ACTIVITY,

		typeTC,
		typePH,
		typeSOLUTION_PE,
		typeMU,
		typeAH2O,
		typeTOTAL_H,
		typeTOTAL_O,
		typeCB,
		typeMASS_WATER,
		typeTOTALS,
		typeTOTAL,
		typeNAME,
		typeMOLES,
		typeACTS,
		typeACT,
		typeLA,
	} m_type;
	XMLCh *Keyword[10];

private:
	std::vector<conc> m_totals;
	std::vector<master_activity> m_acts;
	std::map<const XMLCh*, ElementType, XMLCH_LESS> m_mapXMLCh2Type;
	/**
	std::hash_map<const XMLCh*, ElementType, hash<const XMLCh*>, XMLCH_EQUALS> m_hashmapXMLCh2Type;
	**/
	struct solution* m_solution_ptr;
};

#endif // !defined(AFX_SAXPHREEQCHANDLERS_H__4A69D5F5_2E57_4001_911D_4ABF6F2C0A0B__INCLUDED_)
