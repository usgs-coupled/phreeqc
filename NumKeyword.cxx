#ifdef _DEBUG
#pragma warning(disable : 4786)   // disable truncation warning (Only used by debugger)
#endif
// NumKeyword.cxx: implementation of the cxxNumKeyword class.
//
//////////////////////////////////////////////////////////////////////

#include "NumKeyword.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cxxNumKeyword::cxxNumKeyword()
{
}

cxxNumKeyword::~cxxNumKeyword()
{
}

void cxxNumKeyword::dump_xml(std::ostream& os, unsigned int indent)const
{
        unsigned int i;

        for(i = 0; i < indent + 1; ++i) os << "  ";             
        os << "<n_user>" << this->n_user << "</n_user>" << "\n";

        for(i = 0; i < indent + 1; ++i) os << "  ";             
        os << "<n_user_end>" << this->n_user_end << "</n_user_end>" << "\n";

        for(i = 0; i < indent + 1; ++i) os << "  ";             
        os << "<Description>" << this->description << "</Description>" << "\n";
}

void cxxNumKeyword::read_number_description(CParser& parser)
{
        std::string keyword;
        std::istream::pos_type ptr;

        // skip keyword
        parser.copy_token(keyword, ptr);

        std::istream::pos_type ptr1 = ptr;
        std::string::size_type pos;
        std::string token;
        if (parser.copy_token(token, ptr) != CParser::TT_DIGIT)
        {
                this->n_user = 1;
                this->n_user_end = 1;
        } else {
                std::istringstream iss(token);
                iss >> this->n_user;
                this->n_user_end = this->n_user;
                std::string token1;
                iss >> token1;
                if ( (pos = token1.find_first_of("-")) != std::string::npos ) {
                        token1.replace(pos, 1, " ");
                        std::istringstream iss1(token1);
                        iss1 >> this->n_user_end;
                        //      ptr1 = ptr;
                }
        }
        /*
        else if ( (pos = token.find_first_of("-")) != std::string::npos )
        {
                token.replace(pos, 1, " ");
                std::istringstream iss(token);
                if (!(iss >> this->n_user >> this->n_user_end))
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
                iss >> this->n_user;
                this->n_user_end = this->n_user;
                ptr1 = ptr;
        }
        */
        // reset get position
        //parser.get_iss().seekg(ptr1);

        // skip whitespace
        while (::isspace(parser.get_iss().peek())) parser.get_iss().ignore();

        // copy description
        std::getline(parser.get_iss(), this->description);
}


void cxxNumKeyword::read_number_description(std::istream& is)
{
        // KEYWORD [[1[-20]] [This is the description]]

        // eat keyword
        std::string token;
        is >> token;

        // skip whitespace
        while (::isspace(is.peek())) is.ignore();

        if (::isdigit(is.peek()))
        {
                is >> this->n_user;
                char ch = is.peek();
                if (ch == '-')
                {
                        is >> ch;  // eat '-'
                        is >> this->n_user_end;
                }
                else
                {
                        this->n_user_end = this->n_user;
                }
        }
        else
        {
                this->n_user = this->n_user_end = 1;
        }

        while (::isspace(is.peek())) is.ignore();

        std::getline(is, this->description);
}

