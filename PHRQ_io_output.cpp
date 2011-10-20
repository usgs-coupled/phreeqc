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
#ifndef PHREEQC_CLASS
#include "phrqproto.h"
#endif
#include "phqalloc.h"
#if !defined(PHREEQC_CLASS)
static int forward_output_to_log = 0;
#include "input.h"
#endif


/* ---------------------------------------------------------------------- */
int CLASS_QUALIFIER
warning_msg(const char *err_str)
/* ---------------------------------------------------------------------- */
{
	if (state == TRANSPORT && transport_warnings == FALSE)
		return (OK);
	if (state == ADVECTION && advection_warnings == FALSE)
		return (OK);
	count_warnings++;
	if (pr.warnings >= 0)
	{
		if (count_warnings > pr.warnings)
			return (OK);
	}
	phrq_io->warning_msg(err_str);
	
	return OK;
}
/* ---------------------------------------------------------------------- */
void CLASS_QUALIFIER
screen_msg(const char *err_str)
/* ---------------------------------------------------------------------- */
{
	fprintf(stderr, "%s", err_str);
}
/* ---------------------------------------------------------------------- */
void CLASS_QUALIFIER
echo_msg(const char *str)
/* ---------------------------------------------------------------------- */
{
	if (pr.echo_input == TRUE)
	{
		phrq_io->output_temp_msg(str);
	}
}

///* ---------------------------------------------------------------------- */
//int CLASS_QUALIFIER
//output_msg(const int type, const char *format, ...)
///* ---------------------------------------------------------------------- */
//{
//
//	va_list args;
//	va_start(args, format);
//	int type_local = type;
//	
//	//if (get_forward_output_to_log() && type == PHRQ_io::OUTPUT_MESSAGE)
//	//{
//	//	//type_local = PHRQ_io::OUTPUT_LOG;
//	//	char temp_buff[1000];
//	//	vsnprintf(temp_buff, 999, format, args);
//	//	phrq_io->log_msg(temp_buff);
//	//	return 1;
//	//}
//
//	//switch (type_local)
//	//{
//	//default:
//	//	break;
//	//}
//
//	phrq_io->output_msg(type_local, format, args);
//	va_end(args);
//
//	return OK;
//}
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

///* ---------------------------------------------------------------------- */
//int CLASS_QUALIFIER
//output_fflush(const int type)
///* ---------------------------------------------------------------------- */
//{
//
//	phrq_io->output_fflush(type);
//	return OK;
//}
//
///* ---------------------------------------------------------------------- */
//int CLASS_QUALIFIER
//output_rewind(const int type)
///* ---------------------------------------------------------------------- */
//{
//	phrq_io->output_rewind(type);
//	return OK;
//}
//
//int CLASS_QUALIFIER
//output_close(const int type)
///* ---------------------------------------------------------------------- */
//{
//	phrq_io->output_close(type);
//	return OK;
//}
//
///* ---------------------------------------------------------------------- */
//int CLASS_QUALIFIER
//output_open(const int type, const char *file_name)
///* ---------------------------------------------------------------------- */
//{
//	assert(file_name && strlen(file_name));
//	return phrq_io->output_open(type, file_name);
//}

#if defined(HDF5_CREATE)
extern void HDFWriteHyperSlabV(const char *name, const char *format,
							   va_list argptr);
#endif

#if defined(USE_MPI) && defined(HDF5_CREATE) && defined(MERGE_FILES)
extern int Merge_fpunchf(const int length, const char *format,
						 va_list argptr);
#endif
void CLASS_QUALIFIER
fpunchf_heading(const char *name)
{
	if (pr.punch == TRUE && punch.in == TRUE)
	{
		punch_msg(name);
	}
}
//int CLASS_QUALIFIER
//fpunchf(const char *name, const char *format, ...)
//{
//
//	va_list args;
//	va_start(args, format);
//	phrq_io->fpunchf(name, format, args);
//	va_end(args);
//
//	return OK;
//}
void CLASS_QUALIFIER
fpunchf(const char *name, const char *format, double d)
{
	phrq_io->fpunchf(name, format, d);
}
void CLASS_QUALIFIER
fpunchf(const char *name, const char *format, char * s)
{
	phrq_io->fpunchf(name, format, s);
}
void CLASS_QUALIFIER
fpunchf(const char *name, const char *format, int d)
{
	phrq_io->fpunchf(name, format, d);
}

void CLASS_QUALIFIER
fpunchf_user(int user_index, const char *format, double d)
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
	phrq_io->fpunchf(name, format, d);
}

void CLASS_QUALIFIER
fpunchf_user(int user_index, const char *format, char * d)
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
	phrq_io->fpunchf(name, format, d);
}

int CLASS_QUALIFIER
fpunchf_end_row(const char *format)
{
	//NOOP for Phreeqc
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
	ENTRY item, *found_item;

/*
 *   Prepare error handling
 */
#ifdef PHREEQC_CLASS
	try {
#else
	int errors = setjmp(mark);
	if (errors != 0)
	{
		return errors;
	}
#endif

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
		if (!phrq_io->error_open(argv[4]))
		{
			sprintf(error_string, "Error opening file, %s.", argv[4]);
			warning_msg(error_string);
		}
	}
	else
	{
		phrq_io->error_open(NULL);
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
	phrq_io->Set_input_file(local_input_file);
	screen_msg(sformatf("Input file: %s\n\n", default_name));
	//output_msg(PHRQ_io::OUTPUT_SEND_MESSAGE, "Input file: %s\r\n\r\n", default_name);
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
	phrq_io->Set_output_file(local_output_file);
	screen_msg(sformatf("Output file: %s\n\n", token));
	//output_msg(PHRQ_io::OUTPUT_SEND_MESSAGE, "Output file: %s\r\n\r\n", token);
	strcpy(out_file, token);
/*
 *   Open file for errors
 */
	if (log == TRUE)
	{
		if (!phrq_io->log_open("phreeqc.log"))
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

	phrq_io->close_input();

	pop_cookie();

	if ((local_input_file = fopen(in_file, "r")) == NULL)
	{;
		error_msg("Can't reopen input file.", STOP);
	}
	else
	{
		phrq_io->Set_input_file(local_input_file);
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
		phrq_io->Set_database_file(local_database_file);
	}
	screen_msg(sformatf("Database file: %s\n\n", token));
	//output_msg(PHRQ_io::OUTPUT_SEND_MESSAGE, "Database file: %s\r\n\r\n", token);
	strcpy(db_file, token);

	output_temp_msg(sformatf("   Input file: %s\n", in_file));
	output_temp_msg(sformatf("  Output file: %s\n", out_file));
	output_temp_msg(sformatf("Database file: %s\n\n", token));
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
#ifdef PHREEQC_CLASS
	}
	catch (PhreeqcStop e)
	{
		return get_input_errors();
	}
#endif
	return 0;
}
// ---------------------------------------------------------------------- */
// dump file methods
// ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
bool CLASS_QUALIFIER
dump_open(const char *file_name)
/* ---------------------------------------------------------------------- */
{
	return this->phrq_io->dump_open(file_name);
}
/* ---------------------------------------------------------------------- */
void CLASS_QUALIFIER
dump_fflush(void)
/* ---------------------------------------------------------------------- */
{
	this->phrq_io->dump_fflush();
}
/* ---------------------------------------------------------------------- */
void CLASS_QUALIFIER
dump_close(void)
/* ---------------------------------------------------------------------- */
{
	this->phrq_io->dump_close();
}
/* ---------------------------------------------------------------------- */
void CLASS_QUALIFIER
dump_rewind(void)
/* ---------------------------------------------------------------------- */
{
	this->phrq_io->dump_rewind();
}
/* ---------------------------------------------------------------------- */
bool CLASS_QUALIFIER
dump_isopen(void)
/* ---------------------------------------------------------------------- */
{
	return this->phrq_io->dump_isopen();
}
/* ---------------------------------------------------------------------- */
void CLASS_QUALIFIER
dump_msg(const char * str)
/* ---------------------------------------------------------------------- */
{
	this->phrq_io->dump_msg(str);
}
// ---------------------------------------------------------------------- */
// error file methods
// ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
bool CLASS_QUALIFIER
error_open(const char *file_name)
/* ---------------------------------------------------------------------- */
{
	return this->phrq_io->error_open(file_name);
}
/* ---------------------------------------------------------------------- */
void CLASS_QUALIFIER
error_fflush(void)
/* ---------------------------------------------------------------------- */
{
	this->phrq_io->error_fflush();
}
/* ---------------------------------------------------------------------- */
void CLASS_QUALIFIER
error_close(void)
/* ---------------------------------------------------------------------- */
{
	this->phrq_io->error_close();
}
/* ---------------------------------------------------------------------- */
void CLASS_QUALIFIER
error_rewind(void)
/* ---------------------------------------------------------------------- */
{
	this->phrq_io->error_rewind();
}
/* ---------------------------------------------------------------------- */
bool CLASS_QUALIFIER
error_isopen(void)
/* ---------------------------------------------------------------------- */
{
	return this->phrq_io->error_isopen();
}
/* ---------------------------------------------------------------------- */
void CLASS_QUALIFIER
error_msg(const char *err_str, bool stop)
/* ---------------------------------------------------------------------- */
{
	if (get_input_errors() <= 0)
		input_error = 1;

#if defined MULTICHART
	if (stop)
		chart_handler.End_timer();
#endif

	phrq_io->error_msg(err_str, stop!=0);

	if (stop)
	{
		// TODO, End_time should be moved to catch
#if defined MULTICHART
		chart_handler.End_timer();
#endif
		throw PhreeqcStop();
	}
}
// ---------------------------------------------------------------------- */
// log file methods
// ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
bool CLASS_QUALIFIER
log_open(const char *file_name)
/* ---------------------------------------------------------------------- */
{
	return this->phrq_io->log_open(file_name);
}
/* ---------------------------------------------------------------------- */
void CLASS_QUALIFIER
log_fflush(void)
/* ---------------------------------------------------------------------- */
{
	this->phrq_io->log_fflush();
}
/* ---------------------------------------------------------------------- */
void CLASS_QUALIFIER
log_close(void)
/* ---------------------------------------------------------------------- */
{
	this->phrq_io->log_close();
}
/* ---------------------------------------------------------------------- */
void CLASS_QUALIFIER
log_rewind(void)
/* ---------------------------------------------------------------------- */
{
	this->phrq_io->log_rewind();
}
/* ---------------------------------------------------------------------- */
bool CLASS_QUALIFIER
log_isopen(void)
/* ---------------------------------------------------------------------- */
{
	return this->phrq_io->log_isopen();
}
/* ---------------------------------------------------------------------- */
void CLASS_QUALIFIER
log_msg(const char * str)
/* ---------------------------------------------------------------------- */
{
	return this->phrq_io->log_msg(str);
}
// ---------------------------------------------------------------------- */
// output_temp file methods
// ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
bool CLASS_QUALIFIER
output_temp_open(const char *file_name)
/* ---------------------------------------------------------------------- */
{
	return this->phrq_io->output_temp_open(file_name);
}
/* ---------------------------------------------------------------------- */
void CLASS_QUALIFIER
output_temp_fflush(void)
/* ---------------------------------------------------------------------- */
{
	this->phrq_io->output_temp_fflush();
}
/* ---------------------------------------------------------------------- */
void CLASS_QUALIFIER
output_temp_close(void)
/* ---------------------------------------------------------------------- */
{
	this->phrq_io->output_temp_close();
}
/* ---------------------------------------------------------------------- */
void CLASS_QUALIFIER
output_temp_rewind(void)
/* ---------------------------------------------------------------------- */
{
	this->phrq_io->output_temp_rewind();
}
/* ---------------------------------------------------------------------- */
bool CLASS_QUALIFIER
output_temp_isopen(void)
/* ---------------------------------------------------------------------- */
{
	return this->phrq_io->output_temp_isopen();
}
/* ---------------------------------------------------------------------- */
void CLASS_QUALIFIER
output_temp_msg(const char * str)
/* ---------------------------------------------------------------------- */
{
	if (get_forward_output_to_log())
	{
		phrq_io->log_msg(str);
	}
	else
	{
		phrq_io->output_temp_msg(str);
	}
}
// ---------------------------------------------------------------------- */
// punch file methods
// ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
bool CLASS_QUALIFIER
punch_open(const char *file_name)
/* ---------------------------------------------------------------------- */
{
	return this->phrq_io->punch_open(file_name);
}
/* ---------------------------------------------------------------------- */
void CLASS_QUALIFIER
punch_fflush(void)
/* ---------------------------------------------------------------------- */
{
	this->phrq_io->punch_fflush();
}
/* ---------------------------------------------------------------------- */
void CLASS_QUALIFIER
punch_close(void)
/* ---------------------------------------------------------------------- */
{
	this->phrq_io->punch_close();
}
/* ---------------------------------------------------------------------- */
void CLASS_QUALIFIER
punch_rewind(void)
/* ---------------------------------------------------------------------- */
{
	this->phrq_io->punch_rewind();
}
/* ---------------------------------------------------------------------- */
bool CLASS_QUALIFIER
punch_isopen(void)
/* ---------------------------------------------------------------------- */
{
	return this->phrq_io->punch_isopen();
}
/* ---------------------------------------------------------------------- */
void CLASS_QUALIFIER
punch_msg(const char * str)
/* ---------------------------------------------------------------------- */
{
	return this->phrq_io->punch_msg(str);
}
#endif /* USE_OLD_IO */