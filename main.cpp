#define EXTERNAL
#include "global.h"
#include "output.h"
#include "phrqproto.h"
#include "input.h"
#include <istream>
#include <fstream>
/*#define PHREEQC_XML*/
#ifdef PHREEQC_XML
#include "SAXPhreeqc.h"
extern void SAX_cleanup(void);
#endif

static char const svnid[] = "$Id$";

#ifdef DOS
static int write_banner(void);
#endif

/* ----------------------------------------------------------------------
 *   MAIN
 * ---------------------------------------------------------------------- */
int
main(int argc, char *argv[])
/*
 *   Main program for PHREEQC
 */
{

	int errors;
	void *db_cookie = NULL;
	void *input_cookie = NULL;
#if defined(WIN32_MEMORY_DEBUG)
	int tmpDbgFlag;

	/*
	 * Set the debug-heap flag to keep freed blocks in the
	 * heap's linked list - This will allow us to catch any
	 * inadvertent use of freed memory
	 */
	tmpDbgFlag = _CrtSetDbgFlag(_CRTDBG_REPORT_FLAG);
	/*tmpDbgFlag |= _CRTDBG_DELAY_FREE_MEM_DF;*/
	tmpDbgFlag |= _CRTDBG_LEAK_CHECK_DF;
	/*tmpDbgFlag |= _CRTDBG_CHECK_ALWAYS_DF;*/
	_CrtSetDbgFlag(tmpDbgFlag);
	/*_crtBreakAlloc = 9482;*/
#endif

	if (svnid == NULL)
		fprintf(stderr, " ");
	phast = FALSE;
/*
 *   Add callbacks for error_msg and warning_msg
 */
#ifdef USE_OLD_IO
	if (add_output_callback(phreeqc_handler, NULL) != OK)
	{
		fprintf(stderr, "ERROR: %s\n",
				"NULL pointer returned from malloc or realloc.");
		fprintf(stderr, "ERROR: %s\n", "Program terminating.");
		return -1;
	}
#else
	phrq_io = new PHRQ_io;
#endif

/*
 *   Open input/output files
 */
	errors = process_file_names(argc, argv, &db_cookie, &input_cookie, TRUE);
	if (errors != 0)
	{
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
	if (errors != 0)
	{
		clean_up();
		return errors;
	}

/*
 *   Load database into memory
 */
#if defined(MERGE_INCLUDE_FILES) 
	//errors = read_database(phrq_io.istream_getc, db_cookie);
	set_cookie((std::ifstream *) db_cookie);
	errors = read_database(PHRQ_io::istream_getc, db_cookie);
	clear_cookie();
#else
	errors = read_database(getc_callback, db_cookie);
#endif
	if (errors != 0)
	{
		clean_up();
		return errors;
	}

/*
 *   Read input data for simulation
 */

#if defined(MERGE_INCLUDE_FILES) 
	//errors = run_simulations(phrq_io.istream_getc, input_cookie);
	set_cookie((std::ifstream *)input_cookie);
	errors = run_simulations(PHRQ_io::istream_getc, input_cookie);
	clear_cookie();
#else
	errors = run_simulations(getc_callback, input_cookie);
#endif
	if (errors != 0)
	{
		clean_up();
		return errors;
	}

/*
 *   Display successful status
 */
	errors = do_status();
	if (errors != 0)
	{
		clean_up();
		return errors;
	}
#ifdef PHREEQC_XML
	{
		int n;
		SAX_StartSystem();
		for (n = 0; n < count_solution; ++n)
		{
			SAX_AddSolution(solution[n]);
		}
		SAX_EndSystem();
		SAX_UnpackSolutions(SAX_GetXMLStr(), SAX_GetXMLLength());
	}
#endif
	clean_up();
#ifdef PHREEQC_CPP
	phrq_io->close_input_files();
	phrq_io->close_output_files();
#else
	close_input_files();
	close_output_files();
#endif
#ifdef PHREEQC_XML
	SAX_cleanup();
#endif
	return 0;
}

/* ---------------------------------------------------------------------- */
int
write_banner(void)
/* ---------------------------------------------------------------------- */
{
	char buffer[80];
	int len, indent;
	screen_msg(
			   "              €ﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂ€\n");
	screen_msg(
			   "              ∫                                            ∫\n");

	/* version */
	len = sprintf(buffer, "* PHREEQC-%s *", "@VERSION@");
	indent = (44 - len) / 2;
	screen_msg(sformatf( "%14c∫%*c%s%*c∫\n", ' ', indent, ' ', buffer,
			   44 - indent - len, ' ').c_str());

	screen_msg(,
			   "              ∫                                            ∫\n");
	screen_msg(
			   "              ∫      A hydrogeochemical transport model    ∫\n");
	screen_msg(
			   "              ∫                                            ∫\n");
	screen_msg(
			   "              ∫                    by                      ∫\n");
	screen_msg(
			   "              ∫       D.L. Parkhurst and C.A.J. Appelo     ∫\n");
	screen_msg(
			   "              ∫                                            ∫\n");


	/* date */
	len = sprintf(buffer, "%s", "@VER_DATE@");
	indent = (44 - len) / 2;
	screen_msg(sformatf("%14c∫%*c%s%*c∫\n", ' ', indent, ' ', buffer,
			   44 - indent - len, ' ').c_str());

	screen_msg(
			   "              €‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹€\n\n");

	return 0;
}
