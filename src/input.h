#ifndef _INC_INPUT_H
#define _INC_INPUT_H

#ifdef PHREEQC_IDENT
static char const svnidinput[] = "$Id$";
#endif
typedef int (*PFN_READ_CALLBACK)(void *cookie);

struct read_callback {
	PFN_READ_CALLBACK callback;
	void *cookie;
	int database;
};

int add_char_to_line(int *i, char c);
int check_line_impl(PFN_READ_CALLBACK pfn, void *cookie, const char *string, int allow_empty, int allow_eof, int allow_keyword, int print);
int get_line(PFN_READ_CALLBACK pfn, void *cookie);
int get_logical_line(PFN_READ_CALLBACK pfn, void *cookie, int *l);
int read_database(PFN_READ_CALLBACK pfn, void *cookie);
int run_simulations(PFN_READ_CALLBACK pfn, void *cookie);
int set_read_callback(PFN_READ_CALLBACK pfn, void *cookie, int database);


#endif /* _INC_INPUT_H */
