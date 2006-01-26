// Parser.cpp: implementation of the CParser class.
//
//////////////////////////////////////////////////////////////////////

#include "Parser.h"
#include "utilities.h"
#include <algorithm>    // std::transform
#include <map>          // std::map
#include <cassert>      // assert
#include <iostream>     // std::cout std::cerr

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CParser::CParser(std::istream& input)
: m_input_stream(input), m_output_stream(std::cout), m_error_stream(std::cerr)
, m_input_error(0), m_next_keyword(KT_NONE)
{
	m_line_save.reserve(80);
	m_line.reserve(80);
}

CParser::CParser(std::istream& input, std::ostream& output)
: m_input_stream(input), m_output_stream(output), m_error_stream(std::cerr)
, m_input_error(0), m_next_keyword(KT_NONE)
{
	m_line_save.reserve(80);
	m_line.reserve(80);
}

CParser::CParser(std::istream& input, std::ostream& output, std::ostream& error)
: m_input_stream(input), m_output_stream(output), m_error_stream(error)
, m_input_error(0), m_next_keyword(KT_NONE)
{
	m_line_save.reserve(80);
	m_line.reserve(80);
}

CParser::~CParser()
{
}

CParser::LINE_TYPE CParser::check_line(const std::string& str, bool allow_empty, bool allow_eof, bool allow_keyword, bool print)
{
	LINE_TYPE i;

	// Get line
	do {
		i = get_line();
		// reset iss
		m_line_iss.str(m_line);
		m_line_iss.seekg(0, std::ios_base::beg);
		m_line_iss.clear();


		if (true) // pr.echo_input == TRUE
		{
			if ((print && i != LT_EOF) || i == LT_KEYWORD)
			{
				get_output() << "\t" << m_line_save << "\n";
			}
		}

	} while (i == LT_EMPTY && allow_empty == false);

	// Check eof
	if (i == LT_EOF && allow_eof == false)
	{
		std::ostringstream msg;
		msg << "Unexpected eof while reading " << str << "\nExecution terminated.\n";
		error_msg(msg, OT_STOP);
	}

	// Check keyword
	if (i == LT_KEYWORD && allow_keyword == false)
	{
		std::ostringstream msg;
		msg << "Expected data for " << str << ", but got a keyword ending data block.";
		error_msg(msg, OT_CONTINUE);
		incr_input_error();
	}
	return i;
}

CParser::LINE_TYPE CParser::get_line()
{
	CParser::LINE_TYPE return_value = LT_EMPTY;
	while (return_value == LT_EMPTY)
	{
		//
		// Eliminate all characters after # sign as a comment
		//

		//
		// Get line, check for eof
		//
		if (get_logical_line() == LT_EOF)
		{
			if (!m_input_stream.eof())
			{
				error_msg("Reading input file.", OT_CONTINUE);
				error_msg("istream::get() returned an error.", OT_STOP);
			}
			else
			{
				//{{MOD
				m_line.erase(m_line.begin(), m_line.end()); // m_line.clear();
				//}}MOD
				m_next_keyword = KT_EOF;
				return LT_EOF;
			}
		}

		//
		// Get long lines
		//
		bool empty = true;
		m_line = m_line_save.substr(0, m_line_save.find_first_of('#'));
		for (unsigned int i = 0; i < m_line.size(); ++i)
		{
			if (!::isspace(m_line[i]))
			{
				empty = false;
				break;
			}
		}

		//
		// New line character encountered 
		//
		return_value = (empty ? LT_EMPTY : LT_OK);
	}

	//
	// Determine return_value
	// 
	if (return_value == LT_OK)
	{
		if (check_key(m_line.begin(), m_line.end()))
		{
			return_value = LT_KEYWORD;
		}
		else
		{
			std::string::iterator beg = m_line.begin();
			std::string::iterator end = m_line.end();
			std::string token;
			copy_token(token, beg, end);

			if (token.size() > 1 && token[0] == '-' && ::isalpha(token[1]))
			{
				return_value = LT_OPTION;
			}			
		}
	}
	return return_value;
}

/**
	Reads input stream until end of line, ";", or eof
	stores characters in line_save

	returns:
		EOF on empty line on end of file or
		OK otherwise
*/
CParser::LINE_TYPE CParser::get_logical_line()
{
	int j;
	unsigned int pos;
	char c;

	m_line_save.erase(m_line_save.begin(), m_line_save.end()); // m_line_save.clear();

	while ((j = m_input_stream.get()) != std::char_traits<char>::eof()) {
		c = (char) j;
		if (c == '#') {
			// ignore all chars after # until newline
			do {
				c = (char) j;
				if (c == '\n') {
					break;
				}
				m_line_save += c;
			} while ((j = m_input_stream.get()) != std::char_traits<char>::eof());
		}
		if (c == ';') break;
		if (c == '\n') {
			break;
		}
		if (c == '\\') {
			pos = m_line_save.size();
			m_line_save += c;
			while ((j = m_input_stream.get()) != std::char_traits<char>::eof()) {
				c = (char) j;
				if (c == '\\') {
					pos = m_line_save.size();
					m_line_save += c;
					continue;
				}
				if (c == '\n') {
					// remove '\\'
					for (; pos < m_line_save.size(); pos++) {
						m_line_save[pos] = m_line_save[pos+1];
					}
					m_line_save.erase(m_line_save.size() - 1, 1);
					break;
				}
				m_line_save += c;
				if (!::isspace(j)) break;
			}
		} else {
			m_line_save += c;
		}
	}
 	if (j == std::char_traits<char>::eof() && m_line_save.size() == 0) {
		return(LT_EOF);
	}
	return(LT_OK);
}


//bool CParser::check_key(const std::string::iterator ptr)
bool CParser::check_key(std::string::iterator begin, std::string::iterator end)
{
	static std::map<std::string, KEY_TYPE> s_keyword_map;
	if (s_keyword_map.size() == 0)
	{
		s_keyword_map.insert(std::map<std::string, KEY_TYPE>::value_type("solution", KT_SOLUTION));
		s_keyword_map.insert(std::map<std::string, KEY_TYPE>::value_type("end",      KT_END));
	}

	std::string lowercase;
	copy_token(lowercase, begin, end);
	std::transform(lowercase.begin(), lowercase.end(), lowercase.begin(), tolower);

	m_next_keyword = KT_NONE;
	std::map<std::string, KEY_TYPE>::iterator map_iter = s_keyword_map.find(lowercase);
	if (map_iter == s_keyword_map.end())
		return false;
	m_next_keyword = (*map_iter).second;
	return true;
}

CParser::STATUS_TYPE CParser::check_units(std::string& tot_units, bool alkalinity, bool check_compatibility,
				const std::string& default_units, bool print)
{
/*
 *   Check if legitimate units
 *   Input:
 *           tot_units           character string to check,
 *           alkalinity          true if alkalinity, false if any other total,
 *           check_compatibility true check alk and default units, false otherwise
 *           default_units       character string of default units (check /L, /kg, etc)
 *           print               true print warning messages
 *   Output:
 *           tot_units           standard form for unit
 */
	using utilities::str_tolower;
	using utilities::replace;
	using utilities::squeeze_white;

	static const char *units[] = {
		"Mol/l",       /* 0 */
		"mMol/l",      /* 1 */
		"uMol/l",      /* 2 */
		"g/l",         /* 3 */
		"mg/l",        /* 4 */
		"ug/l",        /* 5 */
		"Mol/kgs",     /* 6 */
		"mMol/kgs",    /* 7 */
		"uMol/kgs",    /* 8 */
		"g/kgs",       /* 9 = ppt */
		"mg/kgs",      /* 10 = ppm */
		"ug/kgs",      /* 11 = ppb */
		"Mol/kgw",     /* 12 = mol/kg H2O */
		"mMol/kgw",    /* 13 = mmol/kg H2O */
		"uMol/kgw",    /* 14 = umol/kg H2O */
		"g/kgw",       /* 15 = mol/kg H2O */
		"mg/kgw",      /* 16 = mmol/kg H2O */
		"ug/kgw",      /* 17 = umol/kg H2O */
		"eq/l",        /* 18 */
		"meq/l",       /* 19 */
		"ueq/l",       /* 20 */
		"eq/kgs",      /* 21 */
		"meq/kgs",     /* 22 */
		"ueq/kgs",     /* 23 */
		"eq/kgw",      /* 24 */
		"meq/kgw",     /* 25 */
		"ueq/kgw",     /* 26 */
	};

	squeeze_white(tot_units);
	str_tolower(tot_units);
	replace("milli",       "m",      tot_units);
	replace("micro",       "u",      tot_units);
	replace("grams",       "g",      tot_units);
	replace("gram",        "g",      tot_units);
	replace("moles",       "Mol",    tot_units);
	replace("mole",        "Mol",    tot_units);
	replace("mol",         "Mol",    tot_units);
	replace("liter",       "l",      tot_units);
	replace("kgh",         "kgw",    tot_units);
	replace("ppt",         "g/kgs",  tot_units);
	replace("ppm",         "mg/kgs", tot_units);
	replace("ppb",         "ug/kgs", tot_units);
	replace("equivalents", "eq",     tot_units);
	replace("equivalent",  "eq",     tot_units);
	replace("equiv",       "eq",     tot_units);

	std::string::size_type end;
	if ((end = tot_units.find("/l")) != std::string::npos) {
		tot_units.resize(end + 2);
	}
	if ((end = tot_units.find("/kgs")) != std::string::npos) {
		tot_units.resize(end + 4);
	}
	if ((end = tot_units.find("/kgw")) != std::string::npos) {
		tot_units.resize(end + 4);
	}

	//
	//  Check if unit in list
	//
	bool found = false;
	for (unsigned int i = 0; i < sizeof(units) / sizeof(char *); ++i) {
		if (tot_units.compare(units[i]) == 0) {
			found = true;
			break;
		}
	}
	if (!found) {
		if (print) {
			std::ostringstream err;
			err << "Unknown unit, " << tot_units;
			error_msg(err, OT_CONTINUE);
		}
		return ERROR;
	}

	//
	//   Check if units are compatible with default_units
	//
	if (check_compatibility == false) return OK;

	//
	//   Special cases for alkalinity
	//
	if (alkalinity == true && tot_units.find("Mol") != std::string::npos) {
		if (print) {
			warning_msg("Alkalinity given in moles, assumed to be equivalents.");
		}
		replace("Mol", "eq", tot_units);
	}
	if (alkalinity == false && tot_units.find("eq") != std::string::npos) {
		if (print) {
			error_msg("Only alkalinity can be entered in equivalents.", OT_CONTINUE);
		}
		return ERROR;
	}

	//
	//  See if default_units are compatible with tot_units
	//
	if (default_units.find("/l")   != std::string::npos && tot_units.find("/l")   != std::string::npos) return OK;
	if (default_units.find("/kgs") != std::string::npos && tot_units.find("/kgs") != std::string::npos) return OK;
	if (default_units.find("/kgw") != std::string::npos && tot_units.find("/kgw") != std::string::npos) return OK;

	std::string str = default_units;
	replace("kgs", "kg solution", str);
	replace("kgs", "kg solution", tot_units);
	replace("kgw", "kg water",    str);
	replace("kgw", "kg water",    tot_units);
	replace("/l",  "/L",          str);
	replace("Mol", "mol",         str);
	replace("/l",  "/L",          tot_units);
	replace("Mol", "mol",         tot_units);

	if (print) {
		std::ostringstream err;
		err << "Units for master species, " << tot_units << ", are not compatible with default units, " << str << ".";
		error_msg(err, OT_CONTINUE);
	}
	return ERROR;
}

CParser::TOKEN_TYPE CParser::token_type(const std::string& token)
{
	if (!token.empty()) {
		if (::isupper(token[0])) {
			return CParser::TT_UPPER;
		} else if (::islower(token[0])) {
			return CParser::TT_LOWER;
		} else if (::isdigit(token[0]) || token[0] == '.' || token[0] == '-') {
			return CParser::TT_DIGIT;
		} else {
			assert(!::isspace(token[0]));
			return CParser::TT_UNKNOWN;
		}
	}
	else {
		return CParser::TT_EMPTY;
	}
}

CParser::TOKEN_TYPE CParser::peek_token()
{
	std::istringstream::pos_type pos = m_line_iss.tellg();
	std::string token;
	m_line_iss >> token;
	m_line_iss.seekg(pos);
	return token_type(token);
}

CParser::TOKEN_TYPE CParser::copy_token(std::string& token, std::string::iterator& begin, std::string::iterator& end)
{
	if (begin != end)
	{
		std::string::iterator b = begin;
		for (; b < end && ::isspace(*b); ++b);

		begin = b;
		for (; begin < end && !::isspace(*begin); ++begin);

		token.assign(b, begin);
	}
	else
	{
		token.resize(0);
	}

	return token_type(token);
}

CParser::TOKEN_TYPE CParser::copy_token(std::string& token, std::istream& is)
{
	is >> token;
	return token_type(token);
}

CParser::TOKEN_TYPE CParser::copy_token(std::string& token, std::istream::pos_type& pos)
{
	m_line_iss.seekg(pos);	
	// m_line_iss >> token;
	if( !(m_line_iss >> token)) {
		token.erase(token.begin(), token.end()); // token.clear();
	}
	pos = m_line_iss.tellg();
	return token_type(token);
}

CParser::FIND_TYPE CParser::find_option(const std::string& item, int *n, const std::vector<std::string>& list, bool exact)
{
	std::string token(item);
	std::transform(token.begin(), token.end(), token.begin(), tolower);
	for (unsigned int i = 0; i < list.size(); i++)
	{
		if (exact == true)
		{
			if (list[i].compare(token) == 0)
			{
				*n = i;
				return FT_OK;
			}
		}
		else
		{
			if(list[i].find(token) == 0)
			{
				*n = i;
				return FT_OK;
			}
		}
	}

	*n = -1;
	return FT_ERROR;
}

// OPTION_TYPE get_option(const char **opt_list, int count_opt_list, char **next_char)
// OPTION_TYPE CParser::get_option(const std::vector<std::string>& opt_list, std::string::iterator& next_char)
int CParser::get_option(const std::vector<std::string>& opt_list, std::string::iterator& next_char)
{
	//
	// Read a line and check for options
	//
	int j;	
	int /*  opt_l,  */ opt;
	//char *opt_ptr;
	std::string::iterator opt_ptr;

	// char option[MAX_LENGTH];
	std::string option;

	//
	// Read line
	//
	LINE_TYPE lt = check_line("get_option", false, true, true, false);
	if (lt == LT_EOF)
	{
		j = OPTION_EOF;
	}
	else if (lt == LT_KEYWORD)
	{
		j = OPTION_KEYWORD;
	}
	else if (lt == LT_OPTION)
	{
		opt_ptr = m_line.begin();
		std::string::iterator end = m_line.end();
		copy_token(option, opt_ptr, end);
		if (find_option(option, &opt, opt_list, false) == CParser::FT_OK)
		{
			j = opt;
			m_line_save.replace(m_line_save.find(option), option.size(), opt_list[opt]);
			m_line.replace(m_line.find(option), option.size(), opt_list[opt]);
			opt_ptr = m_line.begin();
			std::string::iterator end = m_line.end();
			copy_token(option, opt_ptr, end);
			next_char = opt_ptr;
			if (true) // pr.echo_input == TRUE
			{
				if (true) // database_file == NULL
				{
					get_output() << "\t" << m_line_save  << "\n";
				}
			}
		}
		else
		{
			if (true) // (database_file == NULL)
			{
				get_output() << "\t" << m_line_save  << "\n";
			}
			//////std::istringstream err_msg;
			//////err_msg << "Unknown option.";
			//////err_msg << line_save;
			//////error_msg(const std::string& msg, ONERROR_TYPE);

			// error_msg("Unknown option.", CONTINUE);
			// error_msg(line_save, CONTINUE);
			// input_error++;
			std::cerr << "Unknown option." << "\n";
			std::cerr << m_line_save << "\n";

			j = OPTION_ERROR;
			next_char = m_line.begin();
		}
	}
	else
	{
		opt_ptr = m_line.begin();
		std::string::iterator end = m_line.end();
		copy_token(option, opt_ptr, end);
		if (find_option(option, &opt, opt_list, true) == FT_OK)
		{
			j = opt;
			next_char = opt_ptr;
		}
		else
		{
			j = OPTION_DEFAULT;
			next_char = m_line.begin();
		}
		if (true) // pr.echo_input == TRUE
		{
			if (true) // database_file == NULL
			{
				std::cout << "\t" << m_line_save  << "\n";
			}
		}
	}
	return (j);
}

int CParser::get_option(const std::vector<std::string>& opt_list, std::istream::pos_type& next_pos)
{
	//
	// Read a line and check for options
	//
	int j;	
	int opt;
	std::istream::pos_type pos_ptr;
	std::string option;

	//
	// Read line
	//
	LINE_TYPE lt = check_line("get_option", false, true, true, false);
	if (lt == LT_EOF)
	{
		j = OPTION_EOF;
	}
	else if (lt == LT_KEYWORD)
	{
		j = OPTION_KEYWORD;
	}
	else if (lt == LT_OPTION)
	{
		std::string::iterator opt_ptr = m_line.begin();
		std::string::iterator end = m_line.end();
		copy_token(option, opt_ptr, end);
		if (find_option(option.substr(1), &opt, opt_list, false) == FT_OK)
		{
			// replace -option with option
			j = opt;
			m_line_save.replace(m_line_save.find(option), option.size(), opt_list[opt]);
			m_line.replace(m_line.find(option), option.size(), opt_list[opt]);

			// reset iss
			m_line_iss.str(m_line);
			m_line_iss.seekg(0, std::ios_base::beg);
			m_line_iss.clear();

			pos_ptr = 0;
			copy_token(option, pos_ptr);
			next_pos = pos_ptr;
			//{{
			////  m_line_iss.clear();
			//}}
			if (true) // pr.echo_input == TRUE
			{
				if (true) // database_file == NULL
				{
					get_output() << "\t" << m_line_save << "\n";
				}
			}
		}
		else
		{
			if (true) // (database_file == NULL)
			{
				get_output() << "\t" << m_line_save << "\n";
			}
			error_msg("Unknown option.", OT_CONTINUE);
			error_msg(m_line_save.c_str(), OT_CONTINUE);
			incr_input_error();
			j = OPTION_ERROR;
			next_pos = pos_ptr;
		}
	}
	else
	{
		pos_ptr = 0;
		copy_token(option, pos_ptr);
		if (find_option(option, &opt, opt_list, true) == FT_OK)
		{
			j = opt;
			next_pos = pos_ptr;
		}
		else
		{
			j = OPTION_DEFAULT;
			next_pos = 0;
		}
		if (true) // pr.echo_input == TRUE
		{
			if (true) // database_file == NULL
			{
				get_output() << "\t" << m_line_save  << "\n";
			}
		}
	}
	return (j);
}

int CParser::error_msg(const char *err_str, ONERROR_TYPE ot)
{
	m_error_stream  << "ERROR: " << err_str << "\n";
	m_error_stream.flush();

	m_output_stream << "ERROR: " << err_str << "\n";
	m_output_stream.flush();

	if (ot == OT_STOP)
	{
		exit(1);
	}
	return 0;
}

int CParser::warning_msg(const char *err_str)
{
	m_error_stream  << "WARNING: " << err_str << "\n";
	m_error_stream.flush();

	m_output_stream << "WARNING: " << err_str << "\n";
	m_output_stream.flush();

	return 0;
}

CParser::STATUS_TYPE CParser::get_elt(std::string::iterator& begin, const std::string::iterator end, std::string& element)
{
	element.erase(element.begin(), element.end()); // element.clear();

	if (begin == end) {
		error_msg("Empty string in get_elt.  Expected an element name.", OT_CONTINUE);
		return ERROR;
	}

	//
	// Load name into char array element
	//
	char c = *begin;
	++begin;
	element.insert(element.end(), c);  // element.push_back(c);
	if (c == '[') {
		while ( (c = *begin) != ']' ) {
			element.insert(element.end(), c);  // element.push_back(c);
			++begin;
			if ( (c = *begin) == ']' ) {
				element.insert(element.end(), c);  // element.push_back(c);
				++begin;
				break;
			} else if (begin == end) {
				error_msg("No ending bracket (]) for element name", OT_CONTINUE);
				incr_input_error();
				return ERROR;
			}
		}
		while (::islower(c = *begin) || c == '_') {
			element.insert(element.end(), c);  // element.push_back(c);
			++begin;
			if (begin == end) break;
		}
	} else {
		while (::islower(c = *begin) || c == '_') {
			element.insert(element.end(), c);  // element.push_back(c);
			++begin;
			if (begin == end) break;
		}
	}
	return OK;
}

CParser::STATUS_TYPE CParser::parse_couple(std::string& token)
{
	// Parse couple puts redox couples in standard form
	// "+" is removed and couples are rewritten in sort
	// order.

	if (utilities::strcmp_nocase_arg1(token.c_str(), "pe") == 0) {
		utilities::str_tolower(token);
		return OK;
	}

	while ( utilities::replace("+", "", token) );

	std::string::iterator ptr = token.begin();
	std::string elt1;
	get_elt(ptr, token.end(), elt1);

	if (*ptr != '(') {
		std::ostringstream err_msg;
		err_msg << "Element name must be followed by " <<
			"parentheses in redox couple, " << token << ".";
		error_msg(err_msg, OT_CONTINUE);
		incr_input_error();
		return ERROR;
	}

	int paren_count = 1;
	std::string paren1 = "(";
	while ( ptr != token.end() ) {
		++ptr;
		if (*ptr == '/' || ptr == token.end()) {
			std::ostringstream err_msg;
			err_msg << "End of line or  ""/"" encountered before end of parentheses, " <<
				token << ".";
			error_msg(err_msg, OT_CONTINUE);
			return ERROR;
		}
		paren1.insert(paren1.end(), *ptr);  // element.push_back(c);
		if (*ptr == '(') ++paren_count;
		if (*ptr == ')') --paren_count;
		if (paren_count == 0) break;
	}

	++ptr;
	if (ptr == token.end() || *ptr != '/') {
		std::ostringstream err_msg;
		err_msg << " ""/"" must follow parentheses " <<
			"ending first half of redox couple, " << token << ".";
		error_msg(err_msg, OT_CONTINUE);
		return ERROR;
	}
	++ptr;
	std::string elt2;
	get_elt(ptr, token.end(), elt2);
	if (elt1.compare(elt2) != 0) {
		std::ostringstream err_msg;
		err_msg << "Redox couple must be two redox states " <<
			"of the same element, " << token << ".";
		error_msg(err_msg, OT_CONTINUE);
		return ERROR;
	}
	if (*ptr != '(') {
		std::ostringstream err_msg;
		err_msg << "Element name must be followed by "
			"parentheses in redox couple, " << token << ".";
		error_msg(err_msg, OT_CONTINUE);
		incr_input_error();
		return ERROR;
	}
	std::string paren2 = "(";
	paren_count = 1;

	while ( ptr != token.end() ) {
		++ptr;
		if (*ptr == '/' || ptr == token.end()) {
			std::ostringstream err_msg;
			err_msg << "End of line or  ""/"" encountered before end of parentheses, " <<
				token << ".";
			error_msg(err_msg, OT_CONTINUE);
			return ERROR;
		}
		paren2.insert(paren2.end(), *ptr);  // element.push_back(c);
		if (*ptr == '(') ++paren_count;
		if (*ptr == ')') --paren_count;
		if (paren_count == 0) break;
	}
	if (paren1.compare(paren2) <  0) {
		token = elt1 + paren1 + std::string("/") + elt2 + paren2;
	} else if (paren1.compare(paren2) >  0) {
		token = elt2 + paren2 + std::string("/") + elt1 + paren1;
	} else {
		std::ostringstream err_msg;
		err_msg << "Both parts of redox couple are the same, " <<
			token << ".";
		error_msg(err_msg, OT_CONTINUE);
		return ERROR;
	}
	return OK;
}
