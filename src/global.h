#if !defined(_INC_PHREEQC_H)
#define CLASS_QUALIFIER 
#define STATIC static
#define CLASS_STATIC
#endif
#ifdef PHREEQC_IDENT
static char const svnidglobal[] =
	"$Id$";
#endif
#ifndef _INC_GLOBAL_H
#define _INC_GLOBAL_H
#include "phrqtype.h"
typedef enum { kcal, cal, kjoules, joules } DELTA_H_UNIT;
#include "global_structures.h"
/* #define NO_DOS */
/* #define PHREEQ98 *//* PHREEQ98: code for graphical user interface */
#if defined (PHREEQ98) || defined (_MSC_VER) 
#define isnan _isnan
#define isfinite _finite
#else 
#if defined (DJGPP)
/* #define isnan(x) (x != x) */
/* #define isfinite(x) (!(x != x)) */
#define isfinite(x) finite(x)
#endif
#endif

/*
 * uncomment following line, to use default DOS file name for
 * output file
 */
/* #define DOS */

/* ----------------------------------------------------------------------
 *   INCLUDE FILES
 * ---------------------------------------------------------------------- */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <errno.h>
#include <float.h>
#if !defined(_INC_PHREEQC_H)
#include <assert.h>
#endif
#include <setjmp.h>


/* must be defined here and in cl.c */
/* #include <nan.h> */
#ifndef NAN
#   define NAN -99999999
#endif
#define MISSING -9999.999
/* search.h -- declarations for POSIX/SVID-compatible search functions */

/* HSEARCH(3C) */
typedef struct entry
{
	char *key;
	void *data;
} ENTRY;
typedef enum
{ FIND, ENTER } ACTION;

/* TSEARCH(3C) */
typedef enum
{ preorder, postorder, endorder, leaf } VISIT;

/* ----------------------------------------------------------------------
 *   DEFINITIONS
 * ---------------------------------------------------------------------- */
#define F_C_MOL 96493.5			/* C/mol or joule/volt-eq */
#define F_KJ_V_EQ  96.4935		/* kJ/volt-eq */
#define F_KCAL_V_EQ 23.0623		/* kcal/volt-eq */
#define R_LITER_ATM 0.0820597	/* L-atm/deg-mol */
#define R_KCAL_DEG_MOL 0.00198726	/* kcal/deg-mol */
#define R_KJ_DEG_MOL 0.00831470	/* kJ/deg-mol */
#define EPSILON 78.5			/* dialectric constant, dimensionless */
#define EPSILON_ZERO 8.854e-12	/* permittivity of free space, C/V-m = C**2/m-J */
#define JOULES_PER_CALORIE 4.1840
#define AVOGADRO 6.02252e23		/* atoms / mole */


#define TRUE 1
#define FALSE 0
#define OK 1
#define ERROR 0
#define STOP 1
#define CONTINUE 0

#define DISP 2
#define STAG 3
#define NOMIX 4

#define CONVERGED 2
#define MASS_BALANCE 3
/*
  #define OSCILLATE 4
  #define H2O_LIMITS 5
*/
#define REWRITE 2
#define INIT -1

/* check_line values, plus EMPTY, EOF, OK */
#define KEYWORD 3

/* copy_token values */
#define EMPTY 2
#define UPPER 4
#define LOWER 5
#define DIGIT 6
#define UNKNOWN 7
#define OPTION 8

/* species types */
#define AQ 0
#define HPLUS 1
#define H2O 2
#define EMINUS 3
#define SOLID 4
#define EX 5
#define SURF 6
#define SURF_PSI 7
#define SURF_PSI1 8
#define SURF_PSI2 9

/* unknown types */
#define MB 10
#define ALK 11
#define CB 12
#define SOLUTION_PHASE_BOUNDARY 13
#define MU 14
#define AH2O 15
#define MH 16
#define MH2O 17
#define PP 18
#define EXCH 19
#define SURFACE 20
#define SURFACE_CB 21
#define SURFACE_CB1 22
#define SURFACE_CB2 23
#define GAS_MOLES 24
#define S_S_MOLES 25
#define PITZER_GAMMA 26
/* state */
#define INITIALIZE	       0
#define INITIAL_SOLUTION   1
#define INITIAL_EXCHANGE   2
#define INITIAL_SURFACE 3
#define INITIAL_GAS_PHASE  4
#define REACTION		   5
#define INVERSE		 6
#define ADVECTION		 7
#define TRANSPORT		 8
#define PHAST		     9

/* constaints in mass balance */
#define EITHER 0
#define DISSOLVE 1
#define PRECIPITATE -1

/* gas phase type */
#define PRESSURE 1
#define VOLUME 2

#define MAX_PP_ASSEMBLAGE 10	/* default estimate of the number of phase assemblages */
#define MAX_ADD_EQUATIONS 20	/* maximum number of equations added together to reduce eqn to
								   master species */
#define MAX_ELEMENTS 50			/* default estimate of the number of elements */
#define MAX_LENGTH 256			/* maximum number of characters component name */
#define MAX_LINE 80				/* estimate of maximum line length */
#define MAX_LM 3.0				/* maximum log molality allowed in intermediate iterations */
#define MIN_LM -30.0			/* minimum log molality allowed before molality set to zero */
#define MAX_MASS_BALANCE 10		/* initial guess of number mass balance equations for a solution */
#define MAX_MASTER 50			/* default estimate of the number of master species */
#define MAX_ELTS 15				/* default estimate for maximum number of times elements occur in
								   an equation */
#define MAX_PHASES 500			/* initial guess of number of phases defined */
#define MAX_SOLUTION 10			/* The maximum number of solutions allowed */
#define MAX_S 500				/* default estimate for maximum number of species in aqueous model */
#define MAX_STRINGS 3000
#define MAX_SUM_JACOB0 50		/* list used to calculate jacobian */
#define MAX_SUM_JACOB1 500		/* list used to calculate jacobian */
#define MAX_SUM_JACOB2 500		/* list used to calculate jacobian */
#define MAX_SUM_MB 500			/* list used to calculate mass balance sums */
#define MAX_TRXN 16				/* default estimate for maximum number of components in an eqn */
#define MAX_UNKNOWNS 15			/* default estimate for maximum number of unknowns in model */
#define TOL 1e-9				/* tolerance for comparisons of double numbers */
#define LOG_ZERO_MOLALITY -30	/* molalities <= LOG_ZERO_MOLALITY are considered equal to zero */
#define MIN_TOTAL 1e-25
#define MIN_TOTAL_SS MIN_TOTAL
#define MIN_RELATED_SURFACE MIN_TOTAL*100
#define MIN_RELATED_LOG_ACTIVITY -30
/* ----------------------------------------------------------------------
 *   STRUCTURES
 * ---------------------------------------------------------------------- */
enum SURFACE_TYPE
{ UNKNOWN_DL, NO_EDL, DDL, CD_MUSIC, CCM };
enum DIFFUSE_LAYER_TYPE
{ NO_DL, BORKOVEK_DL, DONNAN_DL };
enum SITES_UNITS
{ SITES_ABSOLUTE, SITES_DENSITY };

EXTERNAL struct model last_model;
EXTERNAL int same_model;
EXTERNAL int same_temperature;


EXTERNAL struct punch punch;
/* ----------------------------------------------------------------------
 *   Temperatures
 * ---------------------------------------------------------------------- */

EXTERNAL struct temperature *temperature;
EXTERNAL int count_temperature;
/* ----------------------------------------------------------------------
 *   Surface
 * --------------------------------------------------------------------- */

EXTERNAL int g_iterations;
EXTERNAL LDBLE G_TOL;
EXTERNAL struct surface *surface;
EXTERNAL struct surface *dbg_surface;
EXTERNAL int count_surface;
EXTERNAL int max_surface;
EXTERNAL struct Charge_Group *charge_group;
EXTERNAL int change_surf_count;
EXTERNAL struct Change_Surf *change_surf;

/* ----------------------------------------------------------------------
 *   Exchange
 * ---------------------------------------------------------------------- */

EXTERNAL struct exchange *exchange;
EXTERNAL struct exchange *dbg_exchange;
EXTERNAL int count_exchange;
EXTERNAL int max_exchange;
/* ----------------------------------------------------------------------
 *   Kinetics
 * ---------------------------------------------------------------------- */

EXTERNAL struct kinetics *kinetics;
EXTERNAL struct kinetics *dbg_kinetics;
EXTERNAL int count_kinetics;
EXTERNAL int max_kinetics;

EXTERNAL int count_save_values;
EXTERNAL struct save_values *save_values;

#ifdef SKIP
struct kin_exch
{
	char *exch_name;
	char *phase_name;
	LDBLE phase_proportion;
};
EXTERNAL struct kin_exch *kin_exch;
EXTERNAL int count_kin_exch;
struct kin_surf
{
	char *surf_name;
	char *phase_name;
	LDBLE phase_proportion;
};
EXTERNAL struct kin_surf *kin_surf;
EXTERNAL int count_kin_surf;
#endif
/*----------------------------------------------------------------------
 *   Save
 *---------------------------------------------------------------------- */

EXTERNAL struct save save;
/*----------------------------------------------------------------------
 *   Use
 *---------------------------------------------------------------------- */

EXTERNAL struct Use use;
EXTERNAL struct Use *dbg_use;
/*----------------------------------------------------------------------
 *   Copy
 *---------------------------------------------------------------------- */

EXTERNAL struct copier copy_solution;
EXTERNAL struct copier copy_pp_assemblage;
EXTERNAL struct copier copy_exchange;
EXTERNAL struct copier copy_surface;
EXTERNAL struct copier copy_s_s_assemblage;
EXTERNAL struct copier copy_gas_phase;
EXTERNAL struct copier copy_kinetics;
EXTERNAL struct copier copy_mix;
EXTERNAL struct copier copy_irrev;
EXTERNAL struct copier copy_temperature;


/*----------------------------------------------------------------------
 *   Inverse
 *---------------------------------------------------------------------- */

EXTERNAL struct inverse *inverse;
EXTERNAL int count_inverse;

/*----------------------------------------------------------------------
 *   Mix
 *---------------------------------------------------------------------- */

EXTERNAL struct mix *mix;
EXTERNAL struct mix *dbg_mix;
EXTERNAL int count_mix;
/*----------------------------------------------------------------------
 *   Irreversible reaction
 *---------------------------------------------------------------------- */

EXTERNAL struct irrev *irrev;
EXTERNAL struct irrev *dbg_irrev;
EXTERNAL int count_irrev;
/*----------------------------------------------------------------------
 *   Gas phase
 *---------------------------------------------------------------------- */

EXTERNAL int count_gas_phase;
EXTERNAL int max_gas_phase;
EXTERNAL struct gas_phase *gas_phase;
/*----------------------------------------------------------------------
 *   Solid solution
 *---------------------------------------------------------------------- */

EXTERNAL int count_s_s_assemblage;
EXTERNAL int max_s_s_assemblage;
EXTERNAL struct s_s_assemblage *s_s_assemblage;
/*----------------------------------------------------------------------
 *   Pure-phase assemblage
 *---------------------------------------------------------------------- */

EXTERNAL int count_pp_assemblage;
EXTERNAL int max_pp_assemblage;
EXTERNAL struct pp_assemblage *pp_assemblage;
EXTERNAL struct pp_assemblage *dbg_pp_assemblage;
/*----------------------------------------------------------------------
 *   Species_list
 *---------------------------------------------------------------------- */

EXTERNAL int count_species_list;
EXTERNAL int max_species_list;
EXTERNAL struct species_list *species_list;
/*----------------------------------------------------------------------
 *   Jacobian and Mass balance lists
 *---------------------------------------------------------------------- */

EXTERNAL int count_sum_jacob0;	/* number of elements in sum_jacob0 */
EXTERNAL int max_sum_jacob0;	/* calculated maximum number of elements in sum_jacob0 */
EXTERNAL struct list0 *sum_jacob0;	/* array of pointers to targets and coefficients for array */

EXTERNAL int count_sum_mb1;		/* number of elements in sum_mb1 */
EXTERNAL int max_sum_mb1;		/* calculated maximum number of elements in sum_mb1 */
EXTERNAL struct list1 *sum_mb1;	/* array of pointers to sources and targets for mass
								   balance summations with coef = 1.0 */
EXTERNAL int count_sum_jacob1;	/* number of elements in sum_jacob1 */
EXTERNAL int max_sum_jacob1;	/* calculated maximum number of elements in sum_jacob1 */
EXTERNAL struct list1 *sum_jacob1;	/* array of pointers to sources and targets for array
									   equations with coef = 1.0 */
EXTERNAL int count_sum_mb2;		/* number of elements in sum_mb2 */
EXTERNAL int max_sum_mb2;		/* calculated maximum number of elements in sum_mb2 */
EXTERNAL struct list2 *sum_mb2;	/* array of coefficients and pointers to sources and
								   targets for mass balance summations with coef != 1.0 */
EXTERNAL int count_sum_jacob2;	/* number of elements in sum_jacob2 */
EXTERNAL int max_sum_jacob2;	/* calculated maximum number of elements in sum_jacob2 */
EXTERNAL struct list2 *sum_jacob2;	/* array of coefficients and pointers to sources and
									   targets, coef != 1.0 */
EXTERNAL int count_sum_delta;	/* number of elements in sum_delta */
EXTERNAL int max_sum_delta;		/* calculated maximum number of elements in sum_delta */
EXTERNAL struct list2 *sum_delta;	/* array of pointers to sources, targets and coefficients for
									   summing deltas for mass balance equations */
/*----------------------------------------------------------------------
 *   Solution
 *---------------------------------------------------------------------- */

EXTERNAL struct solution **solution;
EXTERNAL struct solution **dbg_solution;
EXTERNAL int count_solution;
EXTERNAL int max_solution;

#if !defined(PHREEQC_CLASS)
#ifdef MAINSUBS
struct iso iso_defaults[] = {
	{"13C", -10, 1},
	{"13C(4)", -10, 1},
	{"13C(-4)", -50, 5},
	{"34S", 10, 1},
	{"34S(6)", 10, 1},
	{"34S(-2)", -30, 5},
	{"2H", -28, 1},
	{"18O", -5, .1},
	{"87Sr", .71, .01},
	{"11B", 20, 5}
};
int count_iso_defaults = (sizeof(iso_defaults) / sizeof(struct iso));
#else   /* MAINSUBS */
extern struct iso iso_defaults[];
extern int count_iso_defaults;
#endif  /* MAINSUBS */
#else   /* PHREEQC_CLASS */
struct iso *iso_defaults;
int count_iso_defaults;
#endif  /* PHREEQC_CLASS */

/*----------------------------------------------------------------------
 *   Global solution
 *---------------------------------------------------------------------- */
EXTERNAL char *title_x;
EXTERNAL int new_x;
EXTERNAL char *description_x;
EXTERNAL LDBLE tc_x;
EXTERNAL LDBLE tk_x;
EXTERNAL LDBLE ph_x;
EXTERNAL LDBLE solution_pe_x;
EXTERNAL LDBLE mu_x;
EXTERNAL LDBLE ah2o_x;
EXTERNAL LDBLE density_x;
EXTERNAL LDBLE total_h_x;
EXTERNAL LDBLE total_o_x;
EXTERNAL LDBLE cb_x;
EXTERNAL LDBLE total_ions_x;
EXTERNAL LDBLE mass_water_aq_x;
EXTERNAL LDBLE mass_water_surfaces_x;
EXTERNAL LDBLE mass_water_bulk_x;
EXTERNAL char *units_x;
EXTERNAL struct pe_data *pe_x;
EXTERNAL int count_isotopes_x;
EXTERNAL struct isotope *isotopes_x;
EXTERNAL int default_pe_x;
/*EXTERNAL int diffuse_layer_x;*/
EXTERNAL enum DIFFUSE_LAYER_TYPE dl_type_x;
EXTERNAL LDBLE total_carbon;
EXTERNAL LDBLE total_co2;
EXTERNAL LDBLE total_alkalinity;
EXTERNAL LDBLE gfw_water;
EXTERNAL LDBLE step_x;
EXTERNAL LDBLE kin_time_x;
/*----------------------------------------------------------------------
 *   Transport data
 *---------------------------------------------------------------------- */
EXTERNAL int count_cells;
EXTERNAL int count_shifts;
EXTERNAL int ishift;
EXTERNAL int bcon_first;
EXTERNAL int bcon_last;
EXTERNAL int correct_disp;
EXTERNAL LDBLE tempr;
EXTERNAL LDBLE timest;
EXTERNAL int simul_tr;
EXTERNAL LDBLE diffc;
EXTERNAL LDBLE heat_diffc;
EXTERNAL int cell;
EXTERNAL LDBLE mcd_substeps;
EXTERNAL struct stag_data *stag_data;
EXTERNAL int print_modulus;
EXTERNAL int punch_modulus;
EXTERNAL int dump_in;
EXTERNAL int dump_modulus;
EXTERNAL int transport_warnings;
EXTERNAL struct cell_data *cell_data;
EXTERNAL int multi_Dflag;		/* signals calc'n of multicomponent diffusion */
EXTERNAL int interlayer_Dflag;	/* multicomponent diffusion and diffusion through interlayer porosity */
EXTERNAL LDBLE default_Dw;		/* default species diffusion coefficient in water at 25oC, m2/s */
EXTERNAL LDBLE multi_Dpor;		/* uniform porosity of free porewater in solid medium */
EXTERNAL LDBLE interlayer_Dpor;	/* uniform porosity of interlayer space of montmorillonite in solid medium */
EXTERNAL LDBLE multi_Dpor_lim;	/* limiting free porewater porosity where transport stops */
EXTERNAL LDBLE interlayer_Dpor_lim;	/* limiting interlayer porosity where transport stops */
EXTERNAL LDBLE multi_Dn;		/* exponent to calculate pore water diffusion coefficient,
								   Dp = Dw * (multi_Dpor)^multi_Dn */
EXTERNAL LDBLE interlayer_tortf;	/* tortuosity_factor in interlayer porosity,
									   Dpil = Dw / interlayer_tortf */

EXTERNAL int cell_no;
/*----------------------------------------------------------------------
 *   Advection data
 *---------------------------------------------------------------------- */
EXTERNAL int count_ad_cells;
EXTERNAL int count_ad_shifts;
EXTERNAL int print_ad_modulus;
EXTERNAL int punch_ad_modulus;
EXTERNAL int *advection_punch, *advection_print;
EXTERNAL LDBLE advection_kin_time;
EXTERNAL LDBLE advection_kin_time_defined;
EXTERNAL int advection_warnings;

/*----------------------------------------------------------------------
 *   Keywords
 *---------------------------------------------------------------------- */

#if !defined(PHREEQC_CLASS)
#ifdef MAINSUBS
/* list of valid keywords */
struct const_key keyword[] = {
	{"eof", 0},
	{"end", 0},
	{"solution_species", 0},
	{"solution_master_species", 0},
	{"solution", 0},
	{"phases", 0},
	{"pure_phases", 0},
	{"reaction", 0},
	{"mix", 0},
	{"use", 0},
	{"save", 0},
	{"exchange_species", 0},
	{"exchange_master_species", 0},
	{"exchange", 0},
	{"surface_species", 0},
	{"surface_master_species", 0},
	{"surface", 0},
	{"reaction_temperature", 0},
	{"inverse_modeling", 0},
	{"gas_phase", 0},
	{"transport", 0},
	{"debug", 0},
	{"selected_output", 0},
	{"select_output", 0},
	{"knobs", 0},
	{"print", 0},
	{"equilibrium_phases", 0},
	{"equilibria", 0},
	{"equilibrium", 0},
	{"pure", 0},
	{"title", 0},
	{"comment", 0},
	{"advection", 0},
	{"kinetics", 0},
	{"incremental_reactions", 0},
	{"incremental", 0},
	{"rates", 0},
	{"solution_s", 0},
	{"user_print", 0},
	{"user_punch", 0},
	{"solid_solutions", 0},
	{"solid_solution", 0},
	{"solution_spread", 0},
	{"spread_solution", 0},
	{"selected_out", 0},
	{"select_out", 0},
	{"user_graph", 0},
	{"llnl_aqueous_model_parameters", 0},
	{"llnl_aqueous_model", 0},
	{"database", 0},
	{"named_analytical_expression", 0},
	{"named_analytical_expressions", 0},
	{"named_expressions", 0},
	{"named_log_k", 0},
	{"isotopes", 0},
	{"calculate_values", 0},
	{"isotope_ratios", 0},
	{"isotope_alphas", 0},
	{"copy", 0},
	{"pitzer", 0},
	{"sit", 0}
#ifdef PHREEQC_CPP
	,
	{"solution_raw", 0},
	{"exchange_raw", 0},
	{"surface_raw", 0},
	{"equilibrium_phases_raw", 0},
	{"kinetics_raw", 0},
	{"solid_solutions_raw", 0},
	{"gas_phase_raw", 0},
	{"reaction_raw", 0},
	{"mix_raw", 0},
	{"reaction_temperature_raw", 0},
	{"dump", 0},
	{"solution_modify", 0},
	{"equilibrium_phases_modify", 0},
	{"exchange_modify", 0},
	{"surface_modify", 0},
	{"solid_solutions_modify", 0},
	{"gas_phase_modify", 0},
	{"kinetics_modify", 0},
	{"delete", 0},
	{"run_cells", 0}
#endif /* PHREEQC_CPP */
	
};
int NKEYS = (sizeof(keyword) / sizeof(struct const_key));	/* Number of valid keywords */
#else   /* MAINSUBS */
extern struct const_key keyword[];
extern int NKEYS;
#endif  /* MAINSUBS */
#else   /* PHREEQC_CLASS */
struct const_key *keyword;
int NKEYS;
#endif  /* PHREEQC_CLASS */

EXTERNAL struct key *keyword_hash;
EXTERNAL int new_model, new_exchange, new_pp_assemblage, new_surface,
	new_reaction, new_temperature, new_mix, new_solution, new_gas_phase,
	new_inverse, new_punch, new_s_s_assemblage, new_kinetics, new_copy,
	new_pitzer;
/*----------------------------------------------------------------------
 *   Elements
 *---------------------------------------------------------------------- */

EXTERNAL struct element **elements;
EXTERNAL int count_elements;
EXTERNAL int max_elements;
EXTERNAL struct element *element_h_one;

/*----------------------------------------------------------------------
 *   Element List
 *---------------------------------------------------------------------- */

EXTERNAL struct elt_list *elt_list;	/* structure array of working space while reading equations
									   names are in "strings", initially in input order */
EXTERNAL int count_elts;		/* number of elements in elt_list = position of next */
EXTERNAL int max_elts;
/*----------------------------------------------------------------------
 *   Reaction
 *---------------------------------------------------------------------- */

/*----------------------------------------------------------------------
 *   Species
 *---------------------------------------------------------------------- */

EXTERNAL struct logk **logk;
EXTERNAL int count_logk;
EXTERNAL int max_logk;

EXTERNAL char *moles_per_kilogram_string;
EXTERNAL char *pe_string;

EXTERNAL struct species **s;
EXTERNAL int count_s;
EXTERNAL int max_s;

EXTERNAL struct species **s_x;
EXTERNAL int count_s_x;
EXTERNAL int max_s_x;

EXTERNAL struct species *s_h2o;
EXTERNAL struct species *s_hplus;
EXTERNAL struct species *s_h3oplus;
EXTERNAL struct species *s_eminus;
EXTERNAL struct species *s_co3;
EXTERNAL struct species *s_h2;
EXTERNAL struct species *s_o2;
/*----------------------------------------------------------------------
 *   Phases
 *---------------------------------------------------------------------- */

EXTERNAL struct phase **phases;
EXTERNAL int count_phases;
EXTERNAL int max_phases;
/*----------------------------------------------------------------------
 *   Master species
 *---------------------------------------------------------------------- */

EXTERNAL struct master **master;	/* structure array of master species */
EXTERNAL struct master **dbg_master;
EXTERNAL int count_master;
EXTERNAL int max_master;
/*----------------------------------------------------------------------
 *   Unknowns
 *---------------------------------------------------------------------- */

EXTERNAL struct unknown **x;
EXTERNAL int count_unknowns;
EXTERNAL int max_unknowns;

EXTERNAL struct unknown *ah2o_unknown;
EXTERNAL struct unknown *alkalinity_unknown;
EXTERNAL struct unknown *carbon_unknown;
EXTERNAL struct unknown *charge_balance_unknown;
EXTERNAL struct unknown *exchange_unknown;
EXTERNAL struct unknown *mass_hydrogen_unknown;
EXTERNAL struct unknown *mass_oxygen_unknown;
EXTERNAL struct unknown *mb_unknown;
EXTERNAL struct unknown *mu_unknown;
EXTERNAL struct unknown *pe_unknown;
EXTERNAL struct unknown *ph_unknown;
EXTERNAL struct unknown *pure_phase_unknown;
EXTERNAL struct unknown *solution_phase_boundary_unknown;
EXTERNAL struct unknown *surface_unknown;
EXTERNAL struct unknown *gas_unknown;
EXTERNAL struct unknown *s_s_unknown;
/*----------------------------------------------------------------------
 *   Reaction work space
 *---------------------------------------------------------------------- */

EXTERNAL struct reaction_temp trxn;	/* structure array of working space while reading equations
									   species names are in "temp_strings" */
EXTERNAL int count_trxn;		/* number of reactants in trxn = position of next */
EXTERNAL int max_trxn;

EXTERNAL struct unknown_list *mb_unknowns;
EXTERNAL int count_mb_unknowns;
EXTERNAL int max_mb_unknowns;
/* ----------------------------------------------------------------------
 *   Print
 * ---------------------------------------------------------------------- */

EXTERNAL struct prints pr;
EXTERNAL int status_on;
EXTERNAL int count_warnings;

/* ----------------------------------------------------------------------
 *   RATES
 * ---------------------------------------------------------------------- */

EXTERNAL struct rate *rates;
EXTERNAL int count_rates;
EXTERNAL LDBLE rate_m, rate_m0, *rate_p, rate_time, rate_sim_time_start,
	rate_sim_time_end, rate_sim_time, rate_moles, initial_total_time;
EXTERNAL int count_rate_p;
/* ----------------------------------------------------------------------
 *   USER PRINT COMMANDS
 * ---------------------------------------------------------------------- */
EXTERNAL struct rate *user_print;
EXTERNAL struct rate *user_punch;
EXTERNAL char **user_punch_headings;
EXTERNAL int user_punch_count_headings;
#ifdef PHREEQ98
EXTERNAL struct rate *user_graph;
EXTERNAL char **user_graph_headings;
EXTERNAL int user_graph_count_headings;
#endif

/* ----------------------------------------------------------------------
 *   GLOBAL DECLARATIONS
 * ---------------------------------------------------------------------- */
EXTERNAL char error_string[10 * MAX_LENGTH];
EXTERNAL int simulation;
EXTERNAL int state;
EXTERNAL int reaction_step;
EXTERNAL int transport_step;
EXTERNAL int transport_start;
EXTERNAL int advection_step;
EXTERNAL int stop_program;
EXTERNAL int incremental_reactions;

EXTERNAL int count_strings;
EXTERNAL int max_strings;

EXTERNAL LDBLE *array;
EXTERNAL LDBLE *delta;
EXTERNAL LDBLE *residual;

EXTERNAL int input_error;
EXTERNAL int next_keyword;
EXTERNAL int parse_error;
EXTERNAL int paren_count;
EXTERNAL int iterations;
EXTERNAL int gamma_iterations;
EXTERNAL int run_reactions_iterations;

EXTERNAL int max_line;
EXTERNAL char *line;
EXTERNAL char *line_save;

EXTERNAL LDBLE LOG_10;

EXTERNAL int debug_model;
EXTERNAL int debug_prep;
EXTERNAL int debug_set;
EXTERNAL int debug_diffuse_layer;
EXTERNAL int debug_inverse;

EXTERNAL LDBLE inv_tol_default;
EXTERNAL int itmax;
EXTERNAL LDBLE ineq_tol;
EXTERNAL LDBLE convergence_tolerance;
EXTERNAL LDBLE step_size;
EXTERNAL LDBLE pe_step_size;
EXTERNAL LDBLE step_size_now;
EXTERNAL LDBLE pe_step_size_now;
EXTERNAL LDBLE pp_scale;
EXTERNAL LDBLE pp_column_scale;
EXTERNAL int diagonal_scale;	/* 0 not used, 1 used */
EXTERNAL int mass_water_switch;
EXTERNAL int delay_mass_water;
EXTERNAL LDBLE censor;
EXTERNAL int aqueous_only;
EXTERNAL int negative_concentrations;
EXTERNAL int calculating_deriv;
EXTERNAL int numerical_deriv;

EXTERNAL int count_total_steps;
EXTERNAL int phast;
EXTERNAL LDBLE *llnl_temp, *llnl_adh, *llnl_bdh, *llnl_bdot, *llnl_co2_coefs;
EXTERNAL int llnl_count_temp, llnl_count_adh, llnl_count_bdh, llnl_count_bdot,
	llnl_count_co2_coefs;

EXTERNAL char *selected_output_file_name;
EXTERNAL char *dump_file_name;

#ifdef PHREEQCI_GUI
EXTERNAL struct spread_sheet g_spread_sheet;
#endif

/* ---------------------------------------------------------------------- */
/*
 *   Hash definitions
 */
/*
** Constants
*/

# define SegmentSize		    256
# define SegmentSizeShift	  8	/* log2(SegmentSize) */
# define DirectorySize	    256
# define DirectorySizeShift      8	/* log2(DirectorySize)  */
# define Prime1			  37
# define Prime2			  1048583
# define DefaultMaxLoadFactor   5


typedef struct Element
{
	/*
	 ** The user only sees the first two fields,
	 ** as we pretend to pass back only a pointer to ENTRY.
	 ** {S}he doesn't know what else is in here.
	 */
	char *Key;
	char *Data;
	struct Element *Next;		/* secret from user    */
} Element, *Segment;

typedef struct
{
	short p;					/* Next bucket to be split      */
	short maxp;					/* upper bound on p during expansion */
	long KeyCount;				/* current # keys       */
	short SegmentCount;			/* current # segments   */
	short MinLoadFactor;
	short MaxLoadFactor;
	Segment *Directory[DirectorySize];
} HashTable;

typedef unsigned long Address;

EXTERNAL HashTable *strings_hash_table;
EXTERNAL HashTable *elements_hash_table;
EXTERNAL HashTable *species_hash_table;
EXTERNAL HashTable *phases_hash_table;
EXTERNAL HashTable *keyword_hash_table;
EXTERNAL HashTable *logk_hash_table;
EXTERNAL HashTable *master_isotope_hash_table;

#if defined(PHREEQCI_GUI)
#include "../../phreeqci_gui.h"
#endif /* defined(PHREEQCI_GUI) */
/* ----------------------------------------------------------------------
 *   ISOTOPES
 * ---------------------------------------------------------------------- */
EXTERNAL struct name_coef match_tokens[50];
EXTERNAL int count_match_tokens;

EXTERNAL int count_master_isotope;
EXTERNAL struct master_isotope **master_isotope;
EXTERNAL int max_master_isotope;
EXTERNAL int initial_solution_isotopes;

#define OPTION_EOF -1
#define OPTION_KEYWORD -2
#define OPTION_ERROR -3
#define OPTION_DEFAULT -4
#define OPT_1 -5


EXTERNAL int count_calculate_value;
EXTERNAL struct calculate_value **calculate_value;
EXTERNAL int max_calculate_value;
EXTERNAL HashTable *calculate_value_hash_table;


EXTERNAL int count_isotope_ratio;
EXTERNAL struct isotope_ratio **isotope_ratio;
EXTERNAL int max_isotope_ratio;
EXTERNAL HashTable *isotope_ratio_hash_table;

EXTERNAL int count_isotope_alpha;
EXTERNAL struct isotope_alpha **isotope_alpha;
EXTERNAL int max_isotope_alpha;
EXTERNAL HashTable *isotope_alpha_hash_table;

EXTERNAL int phreeqc_mpi_myself;

enum entity_type
{ Solution, Reaction, Exchange, Surface, Gas_phase, Pure_phase, Ss_phase,
	Kinetics, Mix, Temperature, UnKnown
};

EXTERNAL int first_read_input;
EXTERNAL char *user_database;
EXTERNAL int pitzer_model, sit_model, pitzer_pe;
EXTERNAL int full_pitzer, always_full_pitzer, ICON, IC;
EXTERNAL LDBLE COSMOT;
EXTERNAL LDBLE AW;
EXTERNAL int have_punch_name;
/* VP: Density Start */
EXTERNAL int print_density;
/* VP: Density End */

EXTERNAL jmp_buf mark;
EXTERNAL LDBLE *zeros;
EXTERNAL int zeros_max;
#if defined(WIN32)
#include <windows.h>
#endif

#if defined(WIN32_MEMORY_DEBUG)
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#endif

EXTERNAL LDBLE pore_volume;

EXTERNAL struct system_species *sys;
EXTERNAL int count_sys, max_sys;

EXTERNAL LDBLE sys_tot;
EXTERNAL LDBLE AA_basic, BB_basic, CC, I_m, rho_0;
EXTERNAL LDBLE solution_mass, solution_volume;
EXTERNAL LDBLE f_rho(LDBLE rho_old);

/* Collect all statics for PHREEQC_CLASS */
#if defined(PHREEQC_CLASS)
/* basic.c ------------------------------- */
#ifdef PHREEQ98
int colnr, rownr;
#endif
int n_user_punch_index;
#define checking	true
#define varnamelen  20
#define maxdims	    4

typedef Char varnamestring[varnamelen + 1];
typedef Char string255[256];

typedef LDBLE numarray[];
typedef Char *strarray[];

#define forloop			0
#define whileloop       1
#define gosubloop       2
#define tokvar			0
#define toknum			1
#define tokstr			2
#define toksnerr		3
#define tokplus			4
#define tokminus		5
#define toktimes		6
#define tokdiv			7
#define tokup			8
#define toklp			9
#define tokrp			10
#define tokcomma		11
#define toksemi			12
#define tokcolon		13
#define tokeq		14
#define toklt	   15
#define tokgt	   16
#define tokle	   17
#define tokge	   18
#define tokne	   19
#define tokand		20
#define tokor	   21
#define tokxor	  22
#define tokmod	  23
#define toknot	  24
#define toksqr	  25
#define toksqrt	 26
#define toksin	  27
#define tokcos	  28
#define toktan	  29
#define tokarctan       30
#define toklog	  31
#define tokexp	  32
#define tokabs	  33
#define toksgn	  34
#define tokstr_	 35
#define tokval	  36
#define tokchr_	 37
#define tokasc	  38
#define toklen	  39
#define tokmid_	 40
#define tokpeek	 41
#define tokrem	  42
#define toklet	  43
#define tokprint	44
#define tokinput	45
#define tokgoto	 46
#define tokif	   47
#define tokend	  48
#define tokstop	 49
#define tokfor	  50
#define toknext	 51
#define tokwhile	52
#define tokwend	 53
#define tokgosub	54
#define tokreturn       55
#define tokread	 56
#define tokdata	 57
#define tokrestore      58
#define tokgotoxy       59
#define tokon	   60
#define tokdim	  61
#define tokpoke	 62
#define toklist	 63
#define tokrun	  64
#define toknew	  65
#define tokload	 66
#define tokmerge	67
#define toksave	 68
#define tokbye	  69
#define tokdel	  70
#define tokrenum	71
#define tokthen	 72
#define tokelse	 73
#define tokto	   74
#define tokstep	 75
#define toktc	   76
#define tokm0	   77
#define tokm	    78
#define tokparm	 79
#define tokact	  80
#define tokmol	  81
#define tokla	   82
#define toklm	   83
#define toksr	   84
#define toksi	   85
#define toktot	  86
#define toktk	   87
#define toktime	 88
#define toklog10	89
#define toksim_time     90
#define tokequi	 91
#define tokgas	  92
#define tokpunch	93
#define tokkin	  94
#define toks_s	  95
#define tokmu	   96
#define tokalk	  97
#define tokrxn	  98
#define tokdist	 99
#define tokmisc1	100
#define tokmisc2	101
#define tokedl	  102
#define tokstep_no      103
#define toksim_no       104
#define toktotal_time   105
#define tokput	  106
#define tokget	  107
#define tokcharge_balance  109
#define tokpercent_error   110
#ifdef PHREEQ98
#define tokgraph_x	111
#define tokgraph_y	112
#define tokgraph_sy       113
#endif
#define tokcell_no      114
#define tokexists       115
#define toksurf	 116
#define toklk_species   117
#define toklk_named     118
#define toklk_phase     119
#define toksum_species  120
#define toksum_gas      121
#define toksum_s_s      122
#define tokcalc_value   123
#define tokdescription  124
#define toksys	  125
#define tokinstr	126
#define tokltrim	127
#define tokrtrim	128
#define toktrim	 129
#define tokpad	  130
#define tokchange_por   131
#define tokget_por	    132
#define tokosmotic	    133
#define tokchange_surf  134
#define tokporevolume   135
#define toksc	136
#define tokgamma	137
#define toklg	   138
/* VP: Density Start */
#define tokrho	   139
/* VP: Density End */
typedef struct tokenrec
{
	struct tokenrec *next;
	int kind;
	union
	{
		struct varrec *vp;
		LDBLE num;
		Char *sp;
		Char snch;
	} UU;
} tokenrec;

typedef struct linerec
{
	long num, num2;
	tokenrec *txt;
	Char inbuf[MAX_LENGTH];
	struct linerec *next;
} linerec;

typedef struct varrec
{
	varnamestring name;
	struct varrec *next;
	long dims[maxdims];
	char numdims;
	boolean stringvar;
	union
	{
		struct
		{
			LDBLE *arr;
			LDBLE *val, rv;
		} U0;
		struct
		{
			Char **sarr;
			Char **sval, *sv;
		} U1;
	} UU;
} varrec;

typedef struct valrec
{
	boolean stringval;
	union
	{
		LDBLE val;
		Char *sval;
	} UU;
} valrec;

typedef struct looprec
{
	struct looprec *next;
	linerec *homeline;
	tokenrec *hometok;
	int kind;
	union
	{
		struct
		{
			varrec *vp;
			LDBLE max, step;
		} U0;
	} UU;
} looprec;

/* Local variables for exec: */
struct LOC_exec
{
	boolean gotoflag, elseflag;
	tokenrec *t;
};

Char *inbuf;
linerec *linebase;
varrec *varbase;
looprec *loopbase;
long curline;
linerec *stmtline, *dataline;
tokenrec *stmttok, *datatok, *buf;
boolean exitflag;
long EXCP_LINE;
HashTable *command_hash_table;
struct const_key *command;
int NCMDS;

/* cl1.c ------------------------------- */

LDBLE *x_arg, *res_arg, *scratch;
int x_arg_max, res_arg_max, scratch_max;

/* dw.c ------------------------------- */

/* COMMON /QQQQ/ */
LDBLE Q0, Q5;
/* COMMON /ACONST/ */
/*LDBLE GASCON = 0.461522e0, TZ = 647.073e0, AA = 1.e0;*/
LDBLE GASCON, TZ, AA;
LDBLE Z, DZ, Y;
/* COMMON /ELLCON/ */
/*LDBLE G1 = 11.e0, G2 = 44.333333333333e0, GF = 3.5e0;*/
LDBLE G1, G2, GF;
LDBLE B1, B2, B1T, B2T, B1TT, B2TT;

/* integrate.c ------------------------------- */

LDBLE z, xd, alpha;
struct surface_charge *surface_charge_ptr;
int max_row_count, max_column_count;
int carbon;
char **col_name, **row_name;
int count_rows, count_optimize;
int col_phases, col_redox, col_epsilon, col_ph, col_water,
	col_isotopes, col_phase_isotopes;
int row_mb, row_fract, row_charge, row_carbon, row_isotopes,
	row_epsilon, row_isotope_epsilon, row_water;
LDBLE *inv_zero, *array1, *res, *inv_delta1, *delta2, *delta3, *inv_cu,
	*delta_save;
LDBLE *min_delta, *max_delta;
int *iu, *is;
int klmd, nklmd, n2d, kode, iter;
LDBLE toler, error, max_pct, scaled_error;
struct master *master_alk;
int *row_back, *col_back;
unsigned long *good, *bad, *minimal;
int max_good, max_bad, max_minimal;
int count_good, count_bad, count_minimal, count_calls;
unsigned long soln_bits, phase_bits, current_bits, temp_bits;

/* inverse.c ------------------------------- */

FILE *netpath_file;
int count_inverse_models, count_pat_solutions;

/* kinetics.c ------------------------------- */

public:
	void *cvode_kinetics_ptr;
	int cvode_test;
	int cvode_error;
	int cvode_n_user;
	int cvode_n_reactions;
	realtype cvode_step_fraction;
	realtype cvode_rate_sim_time;
	realtype cvode_rate_sim_time_start;
	realtype cvode_last_good_time;
	realtype cvode_prev_good_time;
	N_Vector cvode_last_good_y;
	N_Vector cvode_prev_good_y;
	M_Env kinetics_machEnv;
	N_Vector kinetics_y, kinetics_abstol;
	void *kinetics_cvode_mem;
	struct pp_assemblage *cvode_pp_assemblage_save;
	struct s_s_assemblage *cvode_s_s_assemblage_save;
protected:
	LDBLE *m_original;
	LDBLE *m_temp;

/* model.c ------------------------------- */

LDBLE min_value;
LDBLE model_min_value;
LDBLE *normal, *ineq_array, *inv_res, *cu, *zero, *delta1;
int *inv_iu, *inv_is, *back_eq;
int normal_max, ineq_array_max, res_max, cu_max, zero_max, 
	delta1_max, iu_max, is_max, back_eq_max;

/* output.c ------------------------------- */

#define MAX_CALLBACKS 10
struct output_callback *output_callbacks;
size_t count_output_callback;
int forward_output_to_log;

/* phqalloc.c ------------------------------- */

typedef struct PHRQMemHeader
{
	struct PHRQMemHeader *pNext;	/* memory allocated just after this one */
	struct PHRQMemHeader *pPrev;	/* memory allocated just prior to this one */
	size_t size;				/* memory request + sizeof(PHRQMemHeader) */
#if !defined(NDEBUG)
	char *szFileName;			/* file name */
	int nLine;					/* line number */
	int dummy;					/* alignment */
#endif
} PHRQMemHeader;
PHRQMemHeader *s_pTail;

/* phreeqc_files.c ------------------------------- */

char *default_data_base;
FILE *input_file;
FILE *database_file;
FILE *output;		/* OUTPUT_MESSAGE */
FILE *log_file;	    /* OUTPUT_LOG */
FILE *punch_file;	/* OUTPUT_PUNCH */
FILE *error_file;	/* OUTPUT_ERROR */
FILE *dump_file;	/* OUTPUT_DUMP */
#ifdef PHREEQ98
int outputlinenr;
char *LogFileNameC;
char progress_str[512];
#endif
Anyptr __MallocTemp__;
int P_argc;
char **P_argv;
int P_escapecode;
int P_ioresult;
__p2c_jmp_buf *__top_jb;

/* pitzer.c ------------------------------- */

LDBLE A0;
struct species **spec, **cations, **anions, **neutrals;
int count_cations, count_anions, count_neutrals;
int MAXCATIONS, FIRSTANION, MAXNEUTRAL;
struct pitz_param *mcb0, *mcb1, *mcc0;
int *IPRSNT;
LDBLE *M, *LGAMMA;
LDBLE BK[23], DK[23];
#ifdef PHREEQ98
int connect_simulations, graph_initial_solutions;
int shifts_as_points;
int chart_type;
int ShowChart;
int RowOffset, ColumnOffset;
#endif
LDBLE dummy;

/* print.c ------------------------------- */

#ifdef PHREEQ98
int colnr, rownr;
int graph_initial_solutions;
int prev_advection_step, prev_transport_step;	/*, prev_reaction_step */
/* int shifts_as_points; */
int chart_type;
int AddSeries;
int FirstCallToUSER_GRAPH;
#endif

// sit.c -------------------------------

LDBLE sit_A0;
int sit_count_cations, sit_count_anions, sit_count_neutrals;
int sit_MAXCATIONS, sit_FIRSTANION, sit_MAXNEUTRAL;
int *sit_IPRSNT;
LDBLE *sit_M, *sit_LGAMMA;

/* tidy.c ------------------------------- */

LDBLE a0, a1, kc, kb;

/* tally.c ------------------------------- */


struct tally_buffer *t_buffer;
int tally_count_component;

struct tally *tally_table;
int count_tally_table_columns;
int count_tally_table_rows;

/* transport.c ------------------------------- */


struct sol_D *sol_D;
struct sol_D *sol_D_dbg;
struct J_ij *J_ij, *J_ij_il;
int J_ij_count_spec;

struct M_S *m_s;
int count_m_s;
LDBLE tot1_h, tot1_o, tot2_h, tot2_o;
LDBLE diffc_max, diffc_tr, J_ij_sum;
int transp_surf;
LDBLE *heat_mix_array;
LDBLE *temp1, *temp2;
int nmix, heat_nmix;
LDBLE heat_mix_f_imm, heat_mix_f_m;
int warn_MCD_X, warn_fixed_Surf;

/* utilities.c ------------------------------- */

#ifdef PHREEQ98
int AutoLoadOutputFile, CreateToC;
int ProcessMessages, ShowProgress, ShowProgressWindow, ShowChart;
int outputlinenr;
int stop_calculations;
char err_str98[80];
#endif


#endif /* PHREEQC_CLASS) */
#endif /* _INC_GLOBAL_H  */

