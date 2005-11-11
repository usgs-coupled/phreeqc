#define EXTERNAL
#include "global.h"
#include "output.h"
#include "phrqproto.h"
#include "input.h"


static char const svnid[] = "$Id$";

#ifdef DOS
static int write_banner(void);
#endif

/* ----------------------------------------------------------------------
 *   MAIN
 * ---------------------------------------------------------------------- */
int main(int argc, char *argv[])
/*
 *   Main program for PHREEQC
 */
{
	int errors;
	void *db_cookie = NULL;
	void *input_cookie = NULL;
	if (svnid == NULL) fprintf(stderr," ");
	phast = FALSE;
/*
 *   Add callbacks for error_msg and warning_msg
 */
	if (add_output_callback(phreeqc_handler, NULL) != OK) {
		fprintf(stderr, "ERROR: %s\n", "NULL pointer returned from malloc or realloc.");
		fprintf(stderr, "ERROR: %s\n", "Program terminating.");
		return -1;
	}

/*
 *   Open input/output files
 */
	errors = process_file_names(argc, argv, &db_cookie, &input_cookie, TRUE);
	if (errors != 0) {
		clean_up();
		return errors;
	}
#ifdef DOS
	write_banner(); 
#endif

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
	errors = read_database(getc_callback, db_cookie);
	if (errors != 0) {
		clean_up();
		return errors;
	}

/*
 *   Read input data for simulation
 */
	errors = run_simulations(getc_callback, input_cookie);
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
	close_input_files();
	close_output_files();

	return 0;
}
/* ---------------------------------------------------------------------- */
int write_banner(void)
/* ---------------------------------------------------------------------- */
{
output_msg(OUTPUT_SCREEN, "              ����������������������������������������������\n");
output_msg(OUTPUT_SCREEN, "              �                                            �\n");
output_msg(OUTPUT_SCREEN, "              �              * PHREEQC-@VV@ *              �\n");
output_msg(OUTPUT_SCREEN, "              �                                            �\n");
output_msg(OUTPUT_SCREEN, "              �      A hydrogeochemical transport model    �\n");
output_msg(OUTPUT_SCREEN, "              �                                            �\n");
output_msg(OUTPUT_SCREEN, "              �                    by                      �\n");
output_msg(OUTPUT_SCREEN, "              �       D.L. Parkhurst and C.A.J. Appelo     �\n");
output_msg(OUTPUT_SCREEN, "              �                                            �\n");
output_msg(OUTPUT_SCREEN, "              �                 @V_FIXDATE@                �\n");
output_msg(OUTPUT_SCREEN, "              ����������������������������������������������\n\n");

return 0;
}
