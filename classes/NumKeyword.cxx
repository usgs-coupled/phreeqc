// NumKeyword.cxx: implementation of the CNumKeyword class.
//
//////////////////////////////////////////////////////////////////////

#include "NumKeyword.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CNumKeyword::CNumKeyword()
{
}

CNumKeyword::~CNumKeyword()
{
}

void CNumKeyword::dump_xml(std::ostream& os, unsigned int indent)const
{
	unsigned int i;

	for(i = 0; i < indent + 1; ++i) os << "  ";		
	os << "<n_user>" << m_n_user << "</n_user>" << "\n";

	for(i = 0; i < indent + 1; ++i) os << "  ";		
	os << "<n_user_end>" << m_n_user_end << "</n_user_end>" << "\n";

	for(i = 0; i < indent + 1; ++i) os << "  ";		
	os << "<Description>" << m_description << "</Description>" << "\n";
}

void CNumKeyword::read_number_description(CParser& parser)
{
	// skip keyword
	std::string keyword;
	std::istream::pos_type ptr;
	parser.copy_token(keyword, ptr);

	std::istream::pos_type ptr1 = ptr;
	std::string::size_type pos;
	std::string token;
	if (parser.copy_token(token, ptr) != CParser::TT_DIGIT)
	{
		m_n_user = 1;
		m_n_user_end = 1;
	}
	else if ( (pos = token.find_first_of("-")) != std::string::npos )
	{
		token.replace(pos, 1, " ");
		std::istringstream iss(token);
		if (!(iss >> m_n_user >> m_n_user_end))
		{
			std::ostringstream err_oss;
			if (parser.next_keyword() >= 0)
			{
				err_oss << "Reading number range for " << keyword << ".";
			}
			else
			{
				err_oss << "Reading number range for keyword.";
			}
			parser.error_msg(err_oss, CParser::OT_CONTINUE);
			parser.incr_input_error();
		}
		ptr1 = ptr;
	}
	else
	{
		std::istringstream iss(token);
		iss >> m_n_user;
		m_n_user_end = m_n_user;
		ptr1 = ptr;
	}

	// reset get position
	parser.get_iss().seekg(ptr1);

	// skip whitespace
	while (::isspace(parser.get_iss().peek())) parser.get_iss().ignore();

	// copy description
	std::getline(parser.get_iss(), m_description);
}

void CNumKeyword::read_number_description(std::istream& is)
{
	// KEYWORD [[1[-20]] [This is the description]]

	// eat keyword
	std::string token;
	is >> token;

	// skip whitespace
	while (::isspace(is.peek())) is.ignore();

	if (::isdigit(is.peek()))
	{
		is >> m_n_user;
		char ch = is.peek();
		if (ch == '-')
		{
			is >> ch;  // eat '-'
			is >> m_n_user_end;
		}
		else
		{
			m_n_user_end = m_n_user;
		}
	}
	else
	{
		m_n_user = m_n_user_end = 1;
	}

	while (::isspace(is.peek())) is.ignore();

	std::getline(is, m_description);
}
