#ifdef _DEBUG
#pragma warning(disable : 4786)   // disable truncation warning (Only used by debugger)
#endif
#include "Parser.h"
#include "Solution.h"
#include "Exchange.h"
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
	int return_value;
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

	std::istringstream iss_in(keywordLines);
	std::ostringstream oss_out;
	std::ostringstream oss_err;

	CParser parser(iss_in, oss_out, oss_err);
	//For testing, need to read line to get started
	std::vector<std::string> vopts;
	std::istream::pos_type next_char;
	parser.get_option(vopts, next_char);


	cxxSolution sol;
	sol.read_raw(parser);
	struct solution *soln_ptr = sol.cxxSolution2solution();
	int n;

	/*
	 *  This is not quite right, may not produce sort order
	 */

	if (solution_bsearch(soln_ptr->n_user, &n, FALSE) != NULL) {
		solution_free(solution[n]);
	} else {
		n=count_solution++;
		if (count_solution >= max_solution) {
			space ((void **) ((void *) &(solution)), count_solution, &max_solution, sizeof (struct solution *) );
		}
	}
	solution[n] = soln_ptr;


	return(return_value);
}
/* ---------------------------------------------------------------------- */
int read_exchange_raw (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads EXCHANGE_RAW data block
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
	int return_value;
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
		return_value = check_line("exchange_raw",TRUE,TRUE,TRUE,TRUE);
                                               /* empty, eof, keyword, print */
		if (return_value == EOF || return_value == KEYWORD ) break;
		keywordLines.append(line);
		keywordLines.append("\n");
	}

	std::istringstream iss_in(keywordLines);
	std::ostringstream oss_out;
	std::ostringstream oss_err;

	CParser parser(iss_in, oss_out, oss_err);
	//For testing, need to read line to get started
	std::vector<std::string> vopts;
	std::istream::pos_type next_char;
	parser.get_option(vopts, next_char);

	cxxExchange ex;
	ex.read_raw(parser);
	struct exchange *exchange_ptr = ex.cxxExchange2exchange();
	int n;

	/*
	 *  This is not quite right, may not produce sort order
	 */

	if (exchange_bsearch(exchange_ptr->n_user, &n) != NULL) {
		exchange_free(&exchange[n]);
	} else {
		n=count_exchange++;
		if (count_exchange >= max_exchange) {
			space ((void **) ((void *) &(exchange)), count_exchange, &max_exchange, sizeof (struct exchange *) );
		}
	}
	exchange_copy(exchange_ptr, &exchange[n], exchange_ptr->n_user);
	exchange_free(exchange_ptr);
	free_check_null(exchange_ptr);
	return(return_value);
}
