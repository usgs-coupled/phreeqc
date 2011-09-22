#ifndef USE_OLD_IO
#include <assert.h>
#if !defined(PHREEQC_CLASS)
#define EXTERNAL extern
#include "global.h"
#include <istream>
#include <fstream>
#else
#include "Phreeqc.h"
#endif
#include <setjmp.h>
////#include "output.h"
#ifndef PHREEQC_CLASS
#include "phrqproto.h"
#endif
#include "phqalloc.h"
////static char const svnid[] =
////	"$Id: output.cpp 5605 2011-08-24 16:23:37Z dlpark $";
////#if !defined(PHREEQC_CLASS)
////#define MAX_CALLBACKS 10
////static struct output_callback output_callbacks[MAX_CALLBACKS];
////static size_t count_output_callback = 0;
////static int forward_output_to_log = 0;
////#endif
#if !defined(PHREEQC_CLASS)
static int forward_output_to_log = 0;
#include "input.h"
#endif
/* ---------------------------------------------------------------------- */
int CLASS_QUALIFIER
error_msg(const char *err_str, const int stop, ...)
/* ---------------------------------------------------------------------- */
{
#if !defined(PHREEQC_CLASS)
	extern jmp_buf mark;
#endif

	if (input_error <= 0)
		input_error = 1;
#if defined MULTICHART
	if (stop)
		chart_handler.End_timer(PHREEQC_THIS);
#endif
	va_list args;
	va_start(args, err_str);
	bool bstop = (stop == STOP);
	if (status_on == TRUE)
	{
		phrq_io.output_string(PHRQ_io::OUTPUT_ERROR, "\n");
#ifndef DOS
		status_on = FALSE;
#endif
	}
	phrq_io.phreeqc_handler(PHRQ_io::ACTION_OUTPUT, PHRQ_io::OUTPUT_ERROR, err_str, bstop, "", args);
	//phrq_io.output_handler(PHRQ_io::OUTPUT_ERROR, err_str, bstop, "", args);
	va_end(args);

	if (stop == STOP)
	{
		longjmp(mark, input_error);
	}
	return OK;
}

/* ---------------------------------------------------------------------- */
int CLASS_QUALIFIER
warning_msg(const char *err_str, ...)
/* ---------------------------------------------------------------------- */
{
	va_list args;
	va_start(args, err_str);
	if (count_warnings <= pr.warnings)
	{
		if (status_on == TRUE)
		{
			phrq_io.output_string(PHRQ_io::OUTPUT_ERROR, "\n");
#ifndef DOS
			status_on = FALSE;
#endif
		}
		phrq_io.phreeqc_handler(PHRQ_io::ACTION_OUTPUT, PHRQ_io::OUTPUT_WARNING, err_str, false, "", args);
		//phrq_io.output_handler(PHRQ_io::OUTPUT_WARNING, err_str, false, "", args);
	}
	va_end(args);
	count_warnings++;
	return OK;
}
#ifdef SKIP
/* ---------------------------------------------------------------------- */
int CLASS_QUALIFIER
output_msg(const int type, const char *format, ...)
/* ---------------------------------------------------------------------- */
{
		//if (phast == TRUE)
		//{
		//	if (type == OUTPUT_CHECKLINE && pr.echo_input == TRUE && phreeqc_mpi_myself == 0)
		//	{
		//		va_list args2;
		//		va_start(args2, format);
		//		(output_callbacks[i].callback) (ACTION_OUTPUT, OUTPUT_ECHO, NULL, CONTINUE,
		//										output_callbacks[i].cookie, format,
		//										args2);
		//		va_end(args2);
		//	}
		//}
	if (type == OUTPUT_LOG && pr.logfile != TRUE)
		return OK;

	va_list args1;
	va_start(args1, format);
	
	phrq_io.phreeqc_handler(PHRQ_io::ACTION_OUTPUT, type, NULL, false, format, args1);
	//phrq_io.output_handler(type, NULL, false, format, args1);
	va_end(args1);
	return OK;
}
#endif
/* ---------------------------------------------------------------------- */
int CLASS_QUALIFIER
output_msg(const int type, const char *format, ...)
/* ---------------------------------------------------------------------- */
{
		//if (phast == TRUE)
		//{
		//	if (type == OUTPUT_CHECKLINE && pr.echo_input == TRUE && phreeqc_mpi_myself == 0)
		//	{
		//		va_list args2;
		//		va_start(args2, format);
		//		(output_callbacks[i].callback) (ACTION_OUTPUT, OUTPUT_ECHO, NULL, CONTINUE,
		//										output_callbacks[i].cookie, format,
		//										args2);
		//		va_end(args2);
		//	}
		//}

	va_list args1;
	va_start(args1, format);
	
	//if (get_forward_output_to_log())
	//{
	//	save_output = output_file;
	//	output_file = log_file;
	//}

	switch (type)
	{

	case PHRQ_io::OUTPUT_ERROR:
		if (status_on == TRUE)
		{
			phrq_io.output_string(PHRQ_io::OUTPUT_ERROR, "\n");
#ifndef DOS
			status_on = FALSE;
#endif
		}
		phrq_io.output_handler(type, NULL, false, format, args1);
#if defined MULTICHART
		chart_handler.End_timer(PHREEQC_THIS);
#endif
		break;

	case PHRQ_io::OUTPUT_WARNING:
		if (state == TRANSPORT && transport_warnings == FALSE)
			return (OK);
		if (state == ADVECTION && advection_warnings == FALSE)
			return (OK);
		if (pr.warnings >= 0)
		{
			if (count_warnings > pr.warnings)
				return (OK);
		}
		if (status_on == TRUE)
		{
			phrq_io.output_string(PHRQ_io::OUTPUT_ERROR, "\n");
#ifndef DOS
			status_on = FALSE;
#endif
		}
		phrq_io.output_handler(type, NULL, false, format, args1);
		break;
	case PHRQ_io::OUTPUT_CHECKLINE:
		if (pr.echo_input == TRUE)
		{
			phrq_io.output_handler(type, NULL, false, format, args1);
		}
		break;
	case PHRQ_io::OUTPUT_MESSAGE:
	case PHRQ_io::OUTPUT_BASIC:
	case PHRQ_io::OUTPUT_PUNCH:
	case PHRQ_io::OUTPUT_LOG:
	case PHRQ_io::OUTPUT_SCREEN:
	case PHRQ_io::OUTPUT_STDERR:
	case PHRQ_io::OUTPUT_CVODE:
	case PHRQ_io::OUTPUT_DUMP:
		phrq_io.output_handler(type, NULL, false, format, args1);
		break;
	}

	//if (get_forward_output_to_log())
	//{
	//	output_file = save_output;
	//}
	//return (OK);
	va_end(args1);
	return OK;
}
/* ---------------------------------------------------------------------- */
void CLASS_QUALIFIER
set_forward_output_to_log(int value)
/* ---------------------------------------------------------------------- */
{
	forward_output_to_log = value;
}

/* ---------------------------------------------------------------------- */
int CLASS_QUALIFIER
get_forward_output_to_log(void)
/* ---------------------------------------------------------------------- */
{
	return forward_output_to_log;
}

/* ---------------------------------------------------------------------- */
int CLASS_QUALIFIER
output_fflush(const int type, ...)
/* ---------------------------------------------------------------------- */
{
	int check = OK;

	va_list args;
	va_start(args, type);
	check = phrq_io.phreeqc_handler(PHRQ_io::ACTION_FLUSH, type, NULL, false, NULL, args);	
	//check = phrq_io.fileop_handler(type, fflush);
	va_end(args);

	if (check != OK)
		return (ERROR);

	return (OK);
}

/* ---------------------------------------------------------------------- */
int CLASS_QUALIFIER
output_rewind(const int type, ...)
/* ---------------------------------------------------------------------- */
{
	int check;

	check = OK;

	va_list args;
	va_start(args, type);
	check = phrq_io.phreeqc_handler(PHRQ_io::ACTION_REWIND, type, NULL, false, NULL, args);	
	//check = phrq_io.fileop_handler(type, &PHRQ_io::rewind_wrapper);
	va_end(args);

	if (check != OK)
		return (ERROR);
	return (OK);
}

/* ---------------------------------------------------------------------- */
int CLASS_QUALIFIER
output_close(const int type, ...)
/* ---------------------------------------------------------------------- */
{
	int check = OK;

	va_list args;
	va_start(args, type);
	check = phrq_io.phreeqc_handler(PHRQ_io::ACTION_CLOSE, type, NULL, false, NULL, args);
	va_end(args);

	if (check != OK)
		return (ERROR);

	return (OK);
}

/* ---------------------------------------------------------------------- */
int CLASS_QUALIFIER
output_open(const int type, const char *file_name, ...)
/* ---------------------------------------------------------------------- */
{
	int check = OK;
	assert(file_name && strlen(file_name));

	va_list args;
	va_start(args, file_name);
	check = phrq_io.phreeqc_handler(PHRQ_io::ACTION_OPEN, type, file_name, false, NULL, args);
	va_end(args);

	if (check != OK)
		return (ERROR);

	return (OK);
}

#if defined(HDF5_CREATE)
extern void HDFWriteHyperSlabV(const char *name, const char *format,
							   va_list argptr);
#endif

#if defined(USE_MPI) && defined(HDF5_CREATE) && defined(MERGE_FILES)
extern int Merge_fpunchf(const int length, const char *format,
						 va_list argptr);
#endif

int CLASS_QUALIFIER
fpunchf(const char *name, const char *format, ...)
{

	va_list args;
	va_start(args, format);
	phrq_io.phreeqc_handler(PHRQ_io::ACTION_OUTPUT, PHRQ_io::OUTPUT_PUNCH, name, false, format, args);
	va_end(args);

	return OK;
}

int CLASS_QUALIFIER
fpunchf_user(int user_index, const char *format, ...)
{
	static int s_warning = 0;
	static char buffer[80];
	char *name;

	// check headings
	if (user_index < user_punch_count_headings)
	{
		name = user_punch_headings[user_index];
	}
	else
	{
		if (s_warning == 0)
		{
			sprintf(error_string,
					"USER_PUNCH: Headings count doesn't match number of calls to PUNCH.\n");
			warning_msg(error_string);
			s_warning = 1;
		}
		sprintf(buffer, "no_heading_%d",
				(user_index - user_punch_count_headings) + 1);
		name = buffer;
	}

	va_list args;
	va_start(args, format);
	phrq_io.phreeqc_handler(PHRQ_io::ACTION_OUTPUT, PHRQ_io::OUTPUT_PUNCH, name, false, format, args);
	va_end(args);

	return OK;
}

int CLASS_QUALIFIER
fpunchf_end_row(const char *format, ...)
{

	va_list args;
	va_start(args, format);
	phrq_io.phreeqc_handler(PHRQ_io::ACTION_OUTPUT, PHRQ_io::OUTPUT_PUNCH_END_ROW, "", false, format, args);
	va_end(args);

	return OK;
}
/* ---------------------------------------------------------------------- */
int CLASS_QUALIFIER
process_file_names(int argc, char *argv[], void **db_cookie,
				   void **input_cookie, int log)
/* ---------------------------------------------------------------------- */
{
	int l;
	char token[2 * MAX_LENGTH], default_name[2 * MAX_LENGTH];
	char query[2 * MAX_LENGTH];
	char in_file[2 * MAX_LENGTH], out_file[2 * MAX_LENGTH], db_file[2 * MAX_LENGTH];
	char *env_ptr;
	char *ptr;
	int errors;
	ENTRY item, *found_item;

/*
 *   Prepare error handling
 */
	errors = setjmp(mark);
	if (errors != 0)
	{
		return errors;
	}

/*
 *   Prep for get_line
 */
	max_line = MAX_LINE;
	space((void **) ((void *) &line), INIT, &max_line, sizeof(char));
	space((void **) ((void *) &line_save), INIT, &max_line, sizeof(char));
	hcreate_multi(5, &strings_hash_table);
	hcreate_multi(2, &keyword_hash_table);

/*
 *   Initialize hash table
 */
	keyword_hash = (struct key *) PHRQ_malloc(sizeof(struct key));
	if (keyword_hash == NULL)
	{
		malloc_error();
	}
	else
	{
		keyword_hash->name = string_hsave("database");
		keyword_hash->keycount = 0;
		item.key = keyword_hash->name;
		item.data = keyword_hash;
		found_item = hsearch_multi(keyword_hash_table, item, ENTER);
		if (found_item == NULL)
		{
			sprintf(error_string,
					"Hash table error in keyword initialization.");
			error_msg(error_string, STOP);
		}
	}

/*
 *   Open file for error output
 */
	if (argc > 4)
	{
		if (!phrq_io.open_handler(PHRQ_io::OUTPUT_ERROR, argv[4]))
		{
			sprintf(error_string, "Error opening file, %s.", argv[4]);
			warning_msg(error_string);
		}
	}
	else
	{
		phrq_io.open_handler(PHRQ_io::OUTPUT_ERROR, NULL);
	}

/*
 *   Open user-input file
 */
	strcpy(query, "Name of input file?");
	FILE * local_input_file;
	if (argc <= 1)
	{
		default_name[0] = '\0';
		local_input_file = file_open(query, default_name, "r", FALSE);
	}
	else
	{
		strcpy(default_name, argv[1]);
		local_input_file = file_open(query, default_name, "r", TRUE);
	}
	phrq_io.Set_input_file(local_input_file);
	output_msg(PHRQ_io::OUTPUT_SCREEN, "Input file: %s\n\n", default_name);
	output_msg(PHRQ_io::OUTPUT_SEND_MESSAGE, "Input file: %s\r\n\r\n", default_name);
	strcpy(in_file, default_name);
/*
 *   Open file for output
 */
	strcpy(query, "Name of output file?");

	ptr = default_name;
	copy_token(token, &ptr, &l);
	strcat(token, ".out");
	FILE * local_output_file;
	if (argc <= 1)
	{
		local_output_file = file_open(query, token, "w", FALSE);
	}
	else if (argc == 2)
	{
		local_output_file = file_open(query, token, "w", TRUE);
	}
	else if (argc >= 3)
	{
		strcpy(token, argv[2]);
		local_output_file = file_open(query, token, "w", TRUE);
	}
	phrq_io.Set_output_file(local_output_file);
	output_msg(PHRQ_io::OUTPUT_SCREEN, "Output file: %s\n\n", token);
	output_msg(PHRQ_io::OUTPUT_SEND_MESSAGE, "Output file: %s\r\n\r\n", token);
	strcpy(out_file, token);
/*
 *   Open file for errors
 */
	if (log == TRUE)
	{
		//if ((phrq_io.log_file = fopen("phreeqc.log", "w")) == NULL)
		if (!phrq_io.open_handler(PHRQ_io::OUTPUT_LOG, "phreeqc.log"))
		{
			error_msg("Can't open log file, phreeqc.log.", STOP);
		}
	}
	/*
	 *  Read input file for DATABASE keyword
	 */
	std::ifstream * temp_input = new std::ifstream(in_file, std::ifstream::in);
	set_cookie(temp_input);
	if (get_line(PHRQ_io::istream_getc, temp_input) == KEYWORD)
	{
		ptr = line;
		copy_token(token, &ptr, &l);
		if (strcmp_nocase(token, "database") == 0)
		{
#ifdef PHREEQ98
			user_database = string_duplicate(prefix_database_dir(ptr));
#else
			user_database = string_duplicate(ptr);
#endif
			if (string_trim(user_database) == EMPTY)
			{
				warning_msg("DATABASE file name is missing; default database will be used.");
				user_database = (char *) free_check_null(user_database);
			}
		}
	}
	phrq_io.close_input();

	pop_cookie();

	if ((local_input_file = fopen(in_file, "r")) == NULL)
	{;
		error_msg("Can't reopen input file.", STOP);
	}
	else
	{
		phrq_io.Set_input_file(local_input_file);
	}
/*
 *   Open data base
 */
	strcpy(query, "Name of database file?");
	env_ptr = getenv("PHREEQC_DATABASE");
	if (user_database != NULL)
	{
		strcpy(token, user_database);
	}
	else if (env_ptr != NULL)
	{
		strcpy(token, env_ptr);
	}
	else
	{
		strcpy(token, default_data_base);
	}

	FILE * local_database_file;
	if (argc <= 1)
	{
		local_database_file = file_open(query, token, "r", FALSE);
	}
	else if (argc < 4)
	{
		local_database_file = file_open(query, token, "r", TRUE);
	}
	else if (argc >= 4)
	{
		if (user_database == NULL)
		{
			strcpy(token, argv[3]);
		}
		else
		{
#ifndef PHREEQCI_GUI
			warning_msg
				("Database file from DATABASE keyword is used; command line argument ignored.");
#endif
		}
		local_database_file = file_open(query, token, "r", TRUE);
	}
	if (local_database_file != NULL)
	{
		phrq_io.Set_database_file(local_database_file);
	}
	output_msg(PHRQ_io::OUTPUT_SCREEN, "Database file: %s\n\n", token);
	output_msg(PHRQ_io::OUTPUT_SEND_MESSAGE, "Database file: %s\r\n\r\n", token);
	strcpy(db_file, token);

	output_msg(PHRQ_io::OUTPUT_MESSAGE, "   Input file: %s\n", in_file);
	output_msg(PHRQ_io::OUTPUT_MESSAGE, "  Output file: %s\n", out_file);
	output_msg(PHRQ_io::OUTPUT_MESSAGE, "Database file: %s\n\n", token);
/*
 *   local cleanup
 */
	user_database = (char *) free_check_null(user_database);
	line = (char *) free_check_null(line);
	line_save = (char *) free_check_null(line_save);

	hdestroy_multi(keyword_hash_table);
	keyword_hash = (struct key *) free_check_null(keyword_hash);
	keyword_hash_table = NULL;

	free_hash_strings(strings_hash_table);
	hdestroy_multi(strings_hash_table);
	strings_hash_table = NULL;

	db_stream = new std::ifstream(db_file, std::ifstream::in);
	in_stream = new std::ifstream(in_file, std::ifstream::in);
	*db_cookie = db_stream;
	*input_cookie = in_stream;

	return 0;
}
#endif /* USE_OLD_IO */