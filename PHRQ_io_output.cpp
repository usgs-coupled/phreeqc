#include <assert.h>
#include "Phreeqc.h"
#include "phqalloc.h"

/* ---------------------------------------------------------------------- */
int Phreeqc::
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
	if (phrq_io) phrq_io->warning_msg(err_str);
	
	return OK;
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
screen_msg(const char *err_str)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) phrq_io->screen_msg(err_str);
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
echo_msg(const char *str)
/* ---------------------------------------------------------------------- */
{
	if (pr.echo_input == TRUE)
	{
		if (phrq_io) phrq_io->echo_msg(str);
	}
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
set_forward_output_to_log(int value)
/* ---------------------------------------------------------------------- */
{
	forward_output_to_log = value;
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
get_forward_output_to_log(void)
/* ---------------------------------------------------------------------- */
{
	return forward_output_to_log;
}

#if defined(HDF5_CREATE)
extern void HDFWriteHyperSlabV(const char *name, const char *format,
							   va_list argptr);
#endif

#if defined(USE_MPI) && defined(HDF5_CREATE) && defined(MERGE_FILES)
extern int Merge_fpunchf(const int length, const char *format,
						 va_list argptr);
#endif
void Phreeqc::
fpunchf_heading(const char *name)
{
	if (pr.punch == TRUE && punch.in == TRUE)
	{
		punch_msg(name);
	}
}
void Phreeqc::
fpunchf(const char *name, const char *format, double d)
{
	if (phrq_io) phrq_io->fpunchf(name, format, d);
}
void Phreeqc::
fpunchf(const char *name, const char *format, char * s)
{
	if (phrq_io) phrq_io->fpunchf(name, format, s);
}
void Phreeqc::
fpunchf(const char *name, const char *format, int d)
{
	if (phrq_io) phrq_io->fpunchf(name, format, d);
}

void Phreeqc::
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
	if (phrq_io) phrq_io->fpunchf(name, format, d);
}

void Phreeqc::
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
	if (phrq_io) phrq_io->fpunchf(name, format, d);
}

int Phreeqc::
fpunchf_end_row(const char *format)
{
	if (phrq_io) 
	{
		phrq_io->fpunchf_end_row(format);
	}
	return OK;
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
process_file_names(int argc, char *argv[], std::istream **db_cookie,
				   std::istream **input_cookie, int log)
/* ---------------------------------------------------------------------- */
{
	int l;
	char token[2 * MAX_LENGTH], default_name[2 * MAX_LENGTH];
	char query[2 * MAX_LENGTH];
	char in_file[2 * MAX_LENGTH], out_file[2 * MAX_LENGTH], db_file[2 * MAX_LENGTH];
	char *env_ptr;
	char *ptr;
/*
 *   Prepare error handling
 */
	try {
		if (phrq_io == NULL) 
		{
			std::cerr << "No PHRQ_io output handler defined in process_file_names" << std::endl;
		}
/*
 *   Prep for get_line
 */
		max_line = MAX_LINE;
		space((void **) ((void *) &line), INIT, &max_line, sizeof(char));
		space((void **) ((void *) &line_save), INIT, &max_line, sizeof(char));
		hcreate_multi(5, &strings_hash_table);
/*
 *   Open file for error output
 */
		if (phrq_io)
		{
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
		}
		else 
		{
			std::cerr << "No PHRQ_io output handler defined in process_file_names" << std::endl;
		}

		/*
		*   Open user-input file
		*/
		strcpy(query, "Name of input file?");
		FILE * local_input_file = NULL;
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
		if (phrq_io) phrq_io->Set_input_file(local_input_file);
		screen_msg(sformatf("Input file: %s\n\n", default_name));
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
		if (phrq_io) phrq_io->Set_output_file(local_output_file);
		screen_msg(sformatf("Output file: %s\n\n", token));
		strcpy(out_file, token);
/*
 *   Open file for errors
 */
		if (log == TRUE)
		{
			if (phrq_io)
			{
				if (!phrq_io->log_open("phreeqc.log"))
				{
					error_msg("Can't open log file, phreeqc.log.", STOP);
				}
			}
		}
	/*
	 *  Read input file for DATABASE keyword
	 */
		std::ifstream * temp_input = new std::ifstream(in_file, std::ifstream::in);
		push_istream(temp_input);
		if (get_line() == KEYWORD)
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
		if (phrq_io) phrq_io->close_input();

		pop_istream();

		if ((local_input_file = fopen(in_file, "r")) == NULL)
		{;
		error_msg("Can't reopen input file.", STOP);
		}
		else
		{
			if (phrq_io) phrq_io->Set_input_file(local_input_file);
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
			if (phrq_io) phrq_io->Set_database_file(local_database_file);
		}
		screen_msg(sformatf("Database file: %s\n\n", token));
		strcpy(db_file, token);

		output_msg(sformatf("   Input file: %s\n", in_file));
		output_msg(sformatf("  Output file: %s\n", out_file));
		output_msg(sformatf("Database file: %s\n\n", token));
		/*
		*   local cleanup
		*/
		user_database = (char *) free_check_null(user_database);
		line = (char *) free_check_null(line);
		line_save = (char *) free_check_null(line_save);

		free_hash_strings(strings_hash_table);
		hdestroy_multi(strings_hash_table);
		strings_hash_table = NULL;

		*db_cookie = new std::ifstream(db_file, std::ifstream::in);
		*input_cookie = new std::ifstream(in_file, std::ifstream::in);
	}
	catch (PhreeqcStop e)
	{
		return get_input_errors();
	}
	return 0;
}
// ---------------------------------------------------------------------- */
// dump file methods
// ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
bool Phreeqc::
dump_open(const char *file_name)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io)
		return this->phrq_io->dump_open(file_name);
	return false;
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
dump_fflush(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->dump_fflush();
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
dump_close(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->dump_close();
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
dump_rewind(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->dump_rewind();
}
/* ---------------------------------------------------------------------- */
bool Phreeqc::
dump_isopen(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io)
		return this->phrq_io->dump_isopen();
	return false;
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
dump_msg(const char * str)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->dump_msg(str);
}
// ---------------------------------------------------------------------- */
// error file methods
// ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
bool Phreeqc::
error_open(const char *file_name)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io)
		return this->phrq_io->error_open(file_name);
	return false;
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
error_fflush(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->error_fflush();
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
error_close(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->error_close();
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
error_rewind(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->error_rewind();
}
/* ---------------------------------------------------------------------- */
bool Phreeqc::
error_isopen(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io)
		return this->phrq_io->error_isopen();
	return false;
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
error_msg(const char *err_str, bool stop)
/* ---------------------------------------------------------------------- */
{
	if (get_input_errors() <= 0)
		input_error = 1;
	std::ostringstream msg;
	msg << "ERROR: " << err_str << std::endl;
	if (phrq_io)
	{
		phrq_io->output_msg(msg.str().c_str());
		phrq_io->log_msg(msg.str().c_str());

// COMMENT: {11/23/2011 3:51:53 PM}		phrq_io->error_msg("\n");
		if (status_on)
		{
			phrq_io->screen_msg("\n");
		}
		phrq_io->error_msg(msg.str().c_str(), stop);
	}

	if (stop)
	{
		throw PhreeqcStop();
	}
}
// ---------------------------------------------------------------------- */
// log file methods
// ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
bool Phreeqc::
log_open(const char *file_name)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io)
		return this->phrq_io->log_open(file_name);
	return false;
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
log_fflush(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->log_fflush();
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
log_close(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->log_close();
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
log_rewind(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->log_rewind();
}
/* ---------------------------------------------------------------------- */
bool Phreeqc::
log_isopen(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) 
		return this->phrq_io->log_isopen();
	return false;
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
log_msg(const char * str)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->log_msg(str);
}
// ---------------------------------------------------------------------- */
// output_temp file methods
// ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
bool Phreeqc::
output_open(const char *file_name)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) 
		return this->phrq_io->output_open(file_name);
	return false;
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
output_fflush(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->output_fflush();
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
output_close(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->output_close();
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
output_rewind(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->output_rewind();
}
/* ---------------------------------------------------------------------- */
bool Phreeqc::
output_isopen(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io)
		return this->phrq_io->output_isopen();
	return false;
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
output_msg(const char * str)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io)
	{
		if (get_forward_output_to_log())
		{
			phrq_io->log_msg(str);
		}
		else
		{
			phrq_io->output_msg(str);
		}
	}
}
// ---------------------------------------------------------------------- */
// punch file methods
// ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
bool Phreeqc::
punch_open(const char *file_name)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io)
		return this->phrq_io->punch_open(file_name);
	return false;
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
punch_fflush(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->punch_fflush();
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
punch_close(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->punch_close();
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
punch_rewind(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->punch_rewind();
}
/* ---------------------------------------------------------------------- */
bool Phreeqc::
punch_isopen(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io)
		return this->phrq_io->punch_isopen();
	return false;
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
punch_msg(const char * str)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->punch_msg(str);
}
