#ifdef _DEBUG
#pragma warning(disable : 4786)   // disable truncation warning (Only used by debugger)
#endif
#include <stdlib.h>  // ::tolower
#include <ctype.h>   // ::tolower
#include <algorithm> //std::transform
#include <iostream>     // std::cout std::cerr

#include "Utils.h"
#include "Parser.h"
#include "output.h"
////////////////////////////////////////////////////////////////////////////
int Utilities::strcmp_nocase_arg1(const char *str1, const char *str2)
////////////////////////////////////////////////////////////////////////////
{
        //
        // Compare two strings disregarding case
        //
        int c1, c2;
        while ((c1 = ::tolower(*str1++)) == (c2 = *str2++)) {
                if (c1 == '\0') return(0);
        }
        if (c1 < c2) return(-1);
        return(1);
}

////////////////////////////////////////////////////////////////////////////
int Utilities::strcmp_nocase(const char *str1, const char *str2)
////////////////////////////////////////////////////////////////////////////
{
        //
        // Compare two strings disregarding case
        //
        int c1, c2;
        while ((c1 = ::tolower(*str1++)) == (c2 = ::tolower(*str2++))) {
                if (c1 == '\0') return(0);
        }
        if (c1 < c2) return(-1);
        return(1);
}

////////////////////////////////////////////////////////////////////////////
void Utilities::str_tolower(std::string& str)
////////////////////////////////////////////////////////////////////////////
{
        std::transform(str.begin(), str.end(), str.begin(), tolower);
}



////////////////////////////////////////////////////////////////////////////
bool Utilities::replace(const char* str1, const char* str2, std::string& str)
////////////////////////////////////////////////////////////////////////////
{
        std::string::size_type n = str.find(str1, 0);
        if (n == std::string::npos) return false;

        str.replace(n, ::strlen(str1), str2);
        return true;
}
////////////////////////////////////////////////////////////////////////////
void Utilities::squeeze_white(std::string& s_l)
////////////////////////////////////////////////////////////////////////////
{
        std::string str;
        std::string::iterator beg = s_l.begin();
        std::string::iterator end = s_l.end();
        CParser::copy_token(str, beg, end);
        s_l = str;
}
#ifdef SKIP
////////////////////////////////////////////////////////////////////////////
void Utilities::error_msg (const std::string& err_str, const int stop)
////////////////////////////////////////////////////////////////////////////
{
        //std::cerr << err_str << std::endl;
        output_message(OUTPUT_ERROR, err_str, stop, "", args);
}
#endif
