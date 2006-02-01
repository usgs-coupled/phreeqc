#include "Parser.h"
#include "Solution.h"
#define EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "output.h"
#include "phrqproto.h"

extern int check_line(const char *string, int allow_empty, int allow_eof, int allow_keyword,
		      int print);

/* ---------------------------------------------------------------------- */
int read_solution_raw (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads SOLUTION_RAW data block
 *
 *      Arguments:
 *         none
 *
 *      Returns:
 *         KEYWORD if keyword encountered, input_error may be incremented if
 *                    a keyword is encountered in an unexpected position
 *         EOF     if eof encountered while reading mass balance concentrations
 *         ERROR   if error occurred reading data
 *
 */
	char *ptr, *ptr1;
	int l, title_x_length, line_length;
	int return_value;
	char token[MAX_LENGTH];
	/*
	 *  Accumulate lines in std string
	 */
	std::string keywordLines("");

	keywordLines.append(line);
	keywordLines.append("\n");
/*
 *   Read additonal lines
 */
	for (;;) {
		return_value = check_line("solution_raw",TRUE,TRUE,TRUE,TRUE);
                                               /* empty, eof, keyword, print */
		if (return_value == EOF || return_value == KEYWORD ) break;
		keywordLines.append(line);
		keywordLines.append("\n");
	}

	std::ostringstream oss_out;
	std::ostringstream oss_err;
	std::istringstream iss_in(keywordLines);

	CParser parser(iss_in, oss_out, oss_err);

	return(return_value);
}
