#ifdef _DEBUG
#pragma warning(disable : 4786)   // disable truncation warning (Only used by debugger)
#endif
#include "Parser.h"
#include "Solution.h"
#include "Exchange.h"
#include "Surface.h"
#include "PPassemblage.h"
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
/* ---------------------------------------------------------------------- */
int read_surface_raw (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads SURFACE_RAW data block
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
		return_value = check_line("surface_raw",TRUE,TRUE,TRUE,TRUE);
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

	cxxSurface ex;
	ex.read_raw(parser);
	struct surface *surface_ptr = ex.cxxSurface2surface();
	int n;

	/*
	 *  This is not quite right, may not produce sort order
	 */

	if (surface_bsearch(surface_ptr->n_user, &n) != NULL) {
		surface_free(&surface[n]);
	} else {
		n=count_surface++;
		if (count_surface >= max_surface) {
			space ((void **) ((void *) &(surface)), count_surface, &max_surface, sizeof (struct surface *) );
		}
	}
	surface_copy(surface_ptr, &surface[n], surface_ptr->n_user);
	surface_free(surface_ptr);
	free_check_null(surface_ptr);
	return(return_value);
}
/* ---------------------------------------------------------------------- */
int read_equilibrium_phases_raw (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads EQUILIBRIUM_PHASES_RAW data block
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
		return_value = check_line("equilibrium_phases_raw",TRUE,TRUE,TRUE,TRUE);
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

	cxxPPassemblage ex;
	ex.read_raw(parser);
	struct pp_assemblage *pp_assemblage_ptr = ex.cxxPPassemblage2pp_assemblage();
	int n;

	/*
	 *  This is not quite right, may not produce sort order
	 */

	if (pp_assemblage_bsearch(pp_assemblage_ptr->n_user, &n) != NULL) {
		pp_assemblage_free(&pp_assemblage[n]);
	} else {
		n=count_pp_assemblage++;
		if (count_pp_assemblage >= max_pp_assemblage) {
			space ((void **) ((void *) &(pp_assemblage)), count_pp_assemblage, &max_pp_assemblage, sizeof (struct pp_assemblage *) );
		}
	}
	pp_assemblage_copy(pp_assemblage_ptr, &pp_assemblage[n], pp_assemblage_ptr->n_user);
	pp_assemblage_free(pp_assemblage_ptr);
	free_check_null(pp_assemblage_ptr);
	return(return_value);
}
