#include <cmath>
extern "C" {
#define EXTERNAL
#define PHREEQC_IDENT
#include "global.h"
#include "output.h"
#include "phrqproto.h"
#include "input.h"
#undef PHREEQC_IDENT
}

#include <iostream>
#include <fstream>

int cpp_handler(const int action, const int type, const char *err_str, const int stop, void *cookie, const char *, va_list args);
int get_callback(void *cookie);


// files
std::ifstream *input_file    = NULL;
std::ifstream *database_file = NULL;
std::ostream  *error_file    = NULL;
std::ofstream *output        = NULL;
std::ofstream *log_file      = NULL;
std::ofstream *echo_file     = NULL;
std::ofstream *punch_file    = NULL;
std::ofstream *dump_file     = NULL;

int process_file_names(void);

int main(int argc, char *argv[])
{
	int errors;

	error_file = &std::cerr;

/*
 *   Add callbacks for error_msg and warning_msg
 */
	if (add_output_callback(cpp_handler, NULL) != 1) {
		std::cerr << "ERROR: NULL pointer returned from malloc or realloc.\n";
		std::cerr << "ERROR: Program terminating.\n";
		return -1;
	}

/*
 *   Open input files
 */

/*
 *   Open input/output files
 */
	errors = process_file_names();
	if (errors != 0) {
		clean_up();
		return errors;
	}

/*
 *   Initialize arrays
 */
	errors = do_initialize();
	if (errors != 0) {
		clean_up();
		return errors;
	}

/*
 *   Load database into memory
 */
	errors = read_database(get_callback, (void*)database_file);
	if (errors != 0) {
		clean_up();
		return errors;
	}

/*
 *   Read input data for simulation
 */
	errors = run_simulations(get_callback, (void*)input_file);
	if (errors != 0) {
		clean_up();
		return errors;
	}

/*
 *   Display successful status
 */
	errors = do_status();
	if (errors != 0) {
		clean_up();
		return errors;
	}

	clean_up();


	return 0;
}

/* ---------------------------------------------------------------------- */
int cpp_handler(const int action, const int type, const char *err_str, const int stop, void *cookie, const char *format, va_list args)
/* ---------------------------------------------------------------------- */
{
	char big_buffer[5000];

	if (get_forward_output_to_log()) {
	}

	switch (type) {
	case OUTPUT_ERROR:
		if (status_on == TRUE) {
			if (error_file != NULL) {
				(*error_file) << "\n";
			}
			status_on = FALSE;
		}		
		if (error_file != NULL) {
			(*error_file) << "ERROR: " << err_str << "\n";
		}
		if (output != NULL) {
			(*output) << "ERROR: " << err_str << "\n";
		}
		if (stop == STOP) {
			if (error_file != NULL) {
				(*error_file) << "Stopping.\n";
				(*error_file).flush();
			}
			if (output != NULL) {
				(*output) << "Stopping.\n";
				(*output).flush();
			}
		}
		break;

	case OUTPUT_WARNING:
		if (pr.logfile == TRUE && log_file != NULL) {
			(*log_file) << "WARNING: " << err_str << "\n";
			(*log_file).flush();
		}
		if (state == TRANSPORT && transport_warnings == FALSE) return(OK);
		if (state == ADVECTION && advection_warnings == FALSE) return(OK);
		if (pr.warnings >= 0) {
			if (count_warnings > pr.warnings) return(OK);
		}
		if (status_on == TRUE) {
			if (error_file != NULL) {
				(*error_file) << "\n";
			}
#ifndef DOS
			status_on = FALSE;
#endif
		}		
		if (error_file != NULL) {
			(*error_file) << "WARNING: " << err_str << "\n";
			(*error_file).flush();
		}
		if (phast == TRUE && echo_file != NULL) {
			(*echo_file) << "WARNING: " << err_str << "\n";
		}
		if (output != NULL) {
			(*output) << "WARNING: " << err_str << "\n";
		}
		break;
	case OUTPUT_MESSAGE:
	case OUTPUT_BASIC:
		if (output != NULL) {
			vsprintf(big_buffer, format, args);
			(*output) << big_buffer;
		}
		break;
	case OUTPUT_PUNCH:
		if (punch_file != NULL) {
			vsprintf(big_buffer, format, args);
			(*punch_file) << big_buffer;
		}
		break;
// COMMENT: {10/22/2004 3:10:57 PM}	case OUTPUT_ECHO:
// COMMENT: {9/27/2004 8:00:11 PM}		if (echo_file != NULL) {
// COMMENT: {9/27/2004 8:00:11 PM}			vfprintf(echo_file, format, args);
// COMMENT: {9/27/2004 8:00:11 PM}			if (flush) fflush(echo_file);
// COMMENT: {9/27/2004 8:00:11 PM}		}
		break;
	case OUTPUT_GUI_ERROR:
// COMMENT: {9/27/2004 8:00:11 PM}		if (error_file != NULL) {
// COMMENT: {9/27/2004 8:00:11 PM}			vfprintf(error_file, format, args);
// COMMENT: {9/27/2004 8:00:11 PM}			if (flush) fflush(error_file);
// COMMENT: {9/27/2004 8:00:11 PM}		}
		break;
	case OUTPUT_LOG:
// COMMENT: {9/27/2004 8:00:11 PM}		if (pr.logfile == TRUE && log_file != NULL) {
// COMMENT: {9/27/2004 8:00:11 PM}			vfprintf(log_file, format, args);
// COMMENT: {9/27/2004 8:00:11 PM}			if (flush) fflush(error_file);
// COMMENT: {9/27/2004 8:00:11 PM}		}
		break;
	case OUTPUT_SCREEN:
// COMMENT: {9/27/2004 8:00:11 PM}		if (error_file != NULL) {
// COMMENT: {9/27/2004 8:00:11 PM}			vfprintf(error_file, format, args);
// COMMENT: {9/27/2004 8:00:11 PM}			if (flush) fflush(error_file);
// COMMENT: {9/27/2004 8:00:11 PM}		}
		break;
	case OUTPUT_STDERR:
	case OUTPUT_CVODE:
// COMMENT: {9/27/2004 8:00:11 PM}		if (stderr != NULL) {
// COMMENT: {9/27/2004 8:00:11 PM}			vfprintf(stderr, format, args);
// COMMENT: {9/27/2004 8:00:11 PM}			fflush(stderr);
// COMMENT: {9/27/2004 8:00:11 PM}		}
		break;
	case OUTPUT_DUMP:
// COMMENT: {9/27/2004 8:00:11 PM}		if (dump_file != NULL) {
// COMMENT: {9/27/2004 8:00:11 PM}			vfprintf(dump_file, format, args);
// COMMENT: {9/27/2004 8:00:11 PM}			if (flush) fflush(dump_file);
// COMMENT: {9/27/2004 8:00:11 PM}		}
		break;
	default:
		break;
	}

	if (get_forward_output_to_log()) {
	}
	return OK;
}

/* ---------------------------------------------------------------------- */
int open_output_file(const int type, const char *file_name)
/* ---------------------------------------------------------------------- */
{
	switch (type) {


	case OUTPUT_PUNCH:
		if (punch_file != NULL) {
			delete punch_file;
			punch_file = NULL;
		}
		punch_file = new std::ofstream(file_name);
		if (punch_file == NULL || !(*punch_file).is_open()) {
			delete punch_file;
			punch_file = NULL;
			return ERROR;
		}
		break;

	case OUTPUT_DUMP:
		if (dump_file != NULL) {
			delete dump_file;
			dump_file = NULL;
		}
		dump_file = new std::ofstream(file_name);
		if (dump_file == NULL || !(*dump_file).is_open()) {
			delete dump_file;
			dump_file = NULL;
			return ERROR;
		}
		break;

	default:
		assert(FALSE); // unexpected type
		break;
	}


	return OK;
}

int get_callback(void *cookie)
{
	std::ifstream *pifs = (std::ifstream *)cookie;
	return pifs->get();
}

int process_file_names(void)
{
	char db[]  = "C:\\Program Files\\USGS\\Phreeqc Interactive 2.8\\phreeqc.dat";
	char in[]  = "test2";
	char out[] = "test2.out";
	char log[] = "phreeqc.log";

	database_file = new std::ifstream(db);
	if (!database_file->is_open()) {
		std::cerr << "Can't open file, " << db << ".\n";
		return -1;
	}
	input_file = new std::ifstream(in);
	if (!input_file->is_open()) {
		std::cerr << "Can't open file, " << in << ".\n";
		return -1;
	}
	output = new std::ofstream(out);
	if (!output->is_open()) {
		std::cerr << "Can't open file, " << out << ".\n";
		return -1;
	}
	log_file = new std::ofstream(log);
	if (!log_file->is_open()) {
		error_msg ("Can't open log file, phreeqc.log.", STOP);
	}

	output_msg(OUTPUT_MESSAGE, "   Input file: %s\n", in);
	output_msg(OUTPUT_MESSAGE, "  Output file: %s\n", out);
	output_msg(OUTPUT_MESSAGE, "Database file: %s\n\n", db);

	return 0;
}
