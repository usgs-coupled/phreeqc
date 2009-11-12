#ifndef _INC_PHREEQC_H
#define _INC_PHREEQC_H
#if defined(WIN32)
#include <windows.h>
#endif

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
#include <assert.h>
#include <setjmp.h>
#include "phrqtype.h"

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
/* must be defined here and in cl.c */
/* #include <nan.h> */
#ifndef NAN
#   define NAN -99999999
#endif
#define MISSING -9999.999
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
typedef enum
{ kcal, cal, kjoules, joules } DELTA_H_UNIT;

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
//// Remove static
#define STATIC
#define EXTERNAL

class Phreeqc
{
public:
	Phreeqc(void);
	~Phreeqc(void);


private:

struct _generic_N_Vector;
struct calculate_value;
struct conc;
struct element;
struct exchange;
struct exch_comp;
struct elt_list;
struct gas_phase;
struct gas_comp;
struct inverse;
struct inv_elts;
struct inv_phases;
struct inv_isotope;
struct irrev;
struct isotope;
struct kinetics;
struct kinetics_comp;
struct LOC_exec;
struct master;
struct master_activity;
struct master_isotope;
struct mix;
struct mix_comp;
struct name_coef;
struct output_callback;
struct pe_data;
struct phase;
struct PHRQMemHeader;
struct pitz_param;
struct pp_assemblage;
struct pure_phase;
struct reaction;
struct reaction_temp;
struct rxn_token;
struct rxn_token_temp;
struct solution;
struct species;
struct s_s;
struct s_s_assemblage;
struct s_s_comp;
struct species_diff_layer;
struct surface;
struct surface_comp;
struct surface_charge;
struct surface_diff_layer;
struct theta_param;
struct tokenrec;
struct varrec;
struct unknown;

#define PITZER_EXTERNAL 
#include "pitzer.h"
#include "input.h"
#define KINETICS_EXTERNAL 
#include "sundialstypes.h"		/* definitions of types realtype and                        */
							 /* integertype, and the constant FALSE            */
#include "cvode.h"				/* prototypes for CVodeMalloc, CVode, and            */
							 /* CVodeFree, constants OPT_SIZE, BDF, NEWTON,   */
							 /* SV, SUCCESS, NST,NFE,NSETUPS, NNI, NCFN, NETF */
#include "cvdense.h"			/* prototype for CVDense, constant DENSE_NJE    */
#include "nvector_serial.h"		/* definitions of type N_Vector and macro          */
							 /* NV_Ith_S, prototypes for N_VNew, N_VFree      */
#include "dense.h"				/* definitions of type DenseMat, macro DENSE_ELEM */
#include "nvector.h"
#define extern
#include "p2c.h"
#undef extern 
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
 *   STRUCTURES
 * ---------------------------------------------------------------------- */
enum SURFACE_TYPE
{ UNKNOWN_DL, NO_EDL, DDL, CD_MUSIC, CCM };
enum DIFFUSE_LAYER_TYPE
{ NO_DL, BORKOVEK_DL, DONNAN_DL };
enum SITES_UNITS
{ SITES_ABSOLUTE, SITES_DENSITY };
struct model
{
	int force_prep;
	LDBLE temperature;
	int count_exchange;
	struct master **exchange;

	int count_kinetics;
	struct kinetics *kinetics;

	int count_gas_phase;
	struct phase **gas_phase;

	int count_s_s_assemblage;
	char **s_s_assemblage;

	int count_pp_assemblage;
	struct phase **pp_assemblage;
	char **add_formula;
	LDBLE *si;

	/*int diffuse_layer; */
	/*int edl; */
	enum DIFFUSE_LAYER_TYPE dl_type;
	enum SURFACE_TYPE surface_type;
	int only_counter_ions;
	/*int donnan; */
	LDBLE thickness;
	int count_surface_comp;
	char **surface_comp;
	int count_surface_charge;
	char **surface_charge;
};
struct model last_model;
int same_model;
int same_temperature;
struct name_master
{
	char *name;
	struct master *master;
};
struct name_species
{
	char *name;
	struct species *s;
};
struct name_phase
{
	char *name;
	struct phase *phase;
};
struct punch
{
	int in;
	int new_def;
	struct name_master *totals;
	int count_totals;
	struct name_species *molalities;
	int count_molalities;
	struct name_species *activities;
	int count_activities;
	struct name_phase *pure_phases;
	int count_pure_phases;
	struct name_phase *si;
	int count_si;
	struct name_phase *gases;
	int count_gases;
	struct name_phase *s_s;
	int count_s_s;
	struct name_phase *kinetics;
	int count_kinetics;
	struct name_master *isotopes;
	int count_isotopes;
	struct name_master *calculate_values;
	int count_calculate_values;
	int inverse;
	int sim;
	int state;
	int soln;
	int dist;
	int time;
	int step;
	int rxn;
	int temp;
	int ph;
	int pe;
	int alk;
	int mu;
	int water;
	int high_precision;
	int user_punch;
	int charge_balance;
	int percent_error;
};
struct punch punch;
/* ----------------------------------------------------------------------
 *   Temperatures
 * ---------------------------------------------------------------------- */
struct temperature
{
	int n_user;
	int n_user_end;
	char *description;
	LDBLE *t;
	int count_t;
};
struct temperature *temperature;
int count_temperature;
/* ----------------------------------------------------------------------
 *   Surface
 * --------------------------------------------------------------------- */
struct surface
{
	int n_user;
	int n_user_end;
	int new_def;
	/*int diffuse_layer; */
	/*int edl; */
	int only_counter_ions;
	/*int donnan; */
	enum DIFFUSE_LAYER_TYPE dl_type;
	enum SURFACE_TYPE type;
	enum SITES_UNITS sites_units;
	LDBLE thickness;
	LDBLE debye_lengths;
	LDBLE DDL_viscosity;		/* viscosity relative to pure water */
	LDBLE DDL_limit;			/* limits DDL water to this fraction of bulk water */
	char *description;
	int solution_equilibria;
	int n_solution;
	int count_comps;
	struct surface_comp *comps;
	int count_charge;
	struct surface_charge *charge;
	int related_phases;
	int related_rate;
	int transport;				/* transports comp's and charges if true */
};
struct surface_comp
{
	char *formula;
	struct elt_list *formula_totals;
	LDBLE formula_z;
	LDBLE moles;
	struct master *master;
	struct elt_list *totals;
	LDBLE la;
	int charge;
	LDBLE cb;
	char *phase_name;
	LDBLE phase_proportion;
	char *rate_name;
	LDBLE Dw;					/* diffusion coefficient in water, used in MCD. No transport if 0 */
};
struct surface_charge
{
	char *name;
	LDBLE specific_area;
	LDBLE grams;
	LDBLE charge_balance;
	LDBLE mass_water;
	struct elt_list *diffuse_layer_totals;
	int count_g;
	struct surface_diff_layer *g;	/* stores g and dg/dXd for each ionic charge */
	LDBLE la_psi, la_psi1, la_psi2;
	LDBLE psi, psi1, psi2;
	LDBLE capacitance[2];
	LDBLE sigma0, sigma1, sigma2, sigmaddl;
};
struct surface_diff_layer
{
	LDBLE charge;
	LDBLE g;
	LDBLE dg;
	LDBLE psi_to_z;
};
int g_iterations;
LDBLE G_TOL;
struct surface *surface;
struct surface *dbg_surface;
int count_surface;
int max_surface;
struct Charge_Group
{
	LDBLE z;
	LDBLE eq;
} *charge_group;
int change_surf_count;
struct Change_Surf
{
	char *comp_name;
	LDBLE fraction;
	char *new_comp_name;
	LDBLE new_Dw;
	int cell_no;
	int next;
} *change_surf;
/* ----------------------------------------------------------------------
 *   Exchange
 * ---------------------------------------------------------------------- */
struct exchange
{
	int n_user;
	int n_user_end;
	int new_def;
	char *description;
	int solution_equilibria;
	int n_solution;
	int count_comps;
	struct exch_comp *comps;
	int related_phases;
	int related_rate;
	int pitzer_exchange_gammas;
};
struct exch_comp
{
	char *formula;
	LDBLE formula_z;
	struct elt_list *formula_totals;
	LDBLE moles;
	struct master *master;
	struct elt_list *totals;
	LDBLE la;
	LDBLE charge_balance;
	char *phase_name;
	LDBLE phase_proportion;
	char *rate_name;
};
struct exchange *exchange;
struct exchange *dbg_exchange;
int count_exchange;
int max_exchange;
/* ----------------------------------------------------------------------
 *   Kinetics
 * ---------------------------------------------------------------------- */
struct kinetics
{
	int n_user;
	int n_user_end;
	char *description;
	int count_comps;
	struct kinetics_comp *comps;
	int count_steps;
	LDBLE *steps;
	LDBLE step_divide;
	/*char *units; */
	struct elt_list *totals;
	int rk;
	int bad_step_max;
	int use_cvode;
	int cvode_order;
	int cvode_steps;
};
struct kinetics_comp
{
	char *rate_name;
#ifdef SKIP
	char *formula;
#endif
	struct name_coef *list;
	int count_list;
	/*    struct phase *phase; */
	LDBLE tol;
	LDBLE m;
	LDBLE initial_moles;
	LDBLE m0;
	LDBLE moles;
	int count_c_params;
	char **c_params;
	int count_d_params;
	LDBLE *d_params;
};
struct kinetics *kinetics;
struct kinetics *dbg_kinetics;
int count_kinetics;
int max_kinetics;

struct save_values
{
	LDBLE value;
	int count_subscripts;
	int *subscripts;
};
int count_save_values;
struct save_values *save_values;

#ifdef SKIP
struct kin_exch
{
	char *exch_name;
	char *phase_name;
	LDBLE phase_proportion;
};
struct kin_exch *kin_exch;
int count_kin_exch;
struct kin_surf
{
	char *surf_name;
	char *phase_name;
	LDBLE phase_proportion;
};
struct kin_surf *kin_surf;
int count_kin_surf;
#endif
/*----------------------------------------------------------------------
 *   Save
 *---------------------------------------------------------------------- */
struct save
{
	int solution;
	int n_solution_user;
	int n_solution_user_end;
	int mix;
	int n_mix_user;
	int n_mix_user_end;
	int irrev;
	int n_irrev_user;
	int n_irrev_user_end;
	int pp_assemblage;
	int n_pp_assemblage_user;
	int n_pp_assemblage_user_end;
	int exchange;
	int n_exchange_user;
	int n_exchange_user_end;
	int kinetics;
	int n_kinetics_user;
	int n_kinetics_user_end;
	int surface;
	int n_surface_user;
	int n_surface_user_end;
	int gas_phase;
	int n_gas_phase_user;
	int n_gas_phase_user_end;
	int s_s_assemblage;
	int n_s_s_assemblage_user;
	int n_s_s_assemblage_user_end;
};
struct save save;
/*----------------------------------------------------------------------
 *   Use
 *---------------------------------------------------------------------- */
struct Use
{
	int solution_in;
	int n_solution_user;
	int n_solution;
	struct solution *solution_ptr;

	int pp_assemblage_in;
	int n_pp_assemblage_user;
	int n_pp_assemblage;
	struct pp_assemblage *pp_assemblage_ptr;

	int mix_in;
	int n_mix_user;
	int n_mix;
	struct mix *mix_ptr;
	int n_mix_user_orig;

	int irrev_in;
	int n_irrev_user;
	int n_irrev;
	struct irrev *irrev_ptr;

	int exchange_in;
	int n_exchange_user;
	int n_exchange;
	struct exchange *exchange_ptr;

	int kinetics_in;
	int n_kinetics_user;
	int n_kinetics;
	struct kinetics *kinetics_ptr;

	int surface_in;
	int n_surface_user;
	int n_surface;
	struct surface *surface_ptr;

	int temperature_in;
	int n_temperature_user;
	int n_temperature;
	struct temperature *temperature_ptr;

	int inverse_in;
	int n_inverse_user;
	int n_inverse;
	struct inverse *inverse_ptr;

	int gas_phase_in;
	int n_gas_phase_user;
	int n_gas_phase;
	struct gas_phase *gas_phase_ptr;

	int s_s_assemblage_in;
	int n_s_s_assemblage_user;
	int n_s_s_assemblage;
	struct s_s_assemblage *s_s_assemblage_ptr;

	int trans_in;
	int advect_in;
};
struct Use use;
struct Use *dbg_use;
/*----------------------------------------------------------------------
 *   Copy
 *---------------------------------------------------------------------- */
struct copier
{
	int count;
	int max;
	int *n_user;
	int *start;
	int *end;
};
struct copier copy_solution;
struct copier copy_pp_assemblage;
struct copier copy_exchange;
struct copier copy_surface;
struct copier copy_s_s_assemblage;
struct copier copy_gas_phase;
struct copier copy_kinetics;
struct copier copy_mix;
struct copier copy_irrev;
struct copier copy_temperature;


/*----------------------------------------------------------------------
 *   Inverse
 *---------------------------------------------------------------------- */
struct inverse
{
	int n_user;
	char *description;
	int new_def;
	int minimal;
	int range;
	int mp;
	LDBLE mp_censor;
	LDBLE range_max;
	LDBLE tolerance;
	LDBLE mp_tolerance;
	int count_uncertainties;
	LDBLE *uncertainties;
	int count_ph_uncertainties;
	LDBLE *ph_uncertainties;
#ifdef SKIP
	LDBLE *alk_uncertainties;
#endif
	LDBLE water_uncertainty;
	int mineral_water;
	int carbon;
	LDBLE *dalk_dph;
	LDBLE *dalk_dc;
	int count_solns;
	int *solns;
	int count_force_solns;
	int *force_solns;
	int count_elts;
	struct inv_elts *elts;
	int count_phases;
	struct inv_phases *phases;
	int count_master_list;
	struct master **master_list;
	int count_redox_rxns;
	int count_isotopes;
	struct inv_isotope *isotopes;
	int count_i_u;
	struct inv_isotope *i_u;
	int count_isotope_unknowns;
	struct isotope *isotope_unknowns;
	char *netpath;
	char *pat;
};
struct inv_elts
{
	char *name;
	struct master *master;
	int row;
	int count_uncertainties;
	LDBLE *uncertainties;
};
struct inv_isotope
{
	char *isotope_name;
	LDBLE isotope_number;
	char *elt_name;
	int count_uncertainties;
	LDBLE *uncertainties;
};
struct inv_phases
{
	char *name;
	struct phase *phase;
	int column;
	int constraint;
	int force;
	int count_isotopes;
	struct isotope *isotopes;
};
struct inverse *inverse;
int count_inverse;

/*----------------------------------------------------------------------
 *   Mix
 *---------------------------------------------------------------------- */
struct mix
{
	int n_user;
	int n_user_end;
	char *description;
	int count_comps;
	struct mix_comp *comps;
};
struct mix_comp
{
	int n_solution;
	LDBLE fraction;
};
struct mix *mix;
struct mix *dbg_mix;
int count_mix;
/*----------------------------------------------------------------------
 *   Irreversible reaction
 *---------------------------------------------------------------------- */
struct irrev
{
	int n_user;
	int n_user_end;
	char *description;
	struct name_coef *list;
	struct elt_list *elts;
	LDBLE *steps;
	char *units;
	int count_steps;
	int count_list;
};
struct name_coef
{
	char *name;
	LDBLE coef;
};
struct irrev *irrev;
struct irrev *dbg_irrev;
int count_irrev;
/*----------------------------------------------------------------------
 *   Gas phase
 *---------------------------------------------------------------------- */
struct gas_phase
{
	int n_user;
	int n_user_end;
	char *description;
	int new_def;
	int solution_equilibria;
	int n_solution;
	int type;
	LDBLE total_p;
	LDBLE total_moles;
	LDBLE volume;
	LDBLE temperature;
	int count_comps;
	struct gas_comp *comps;
};
struct gas_comp
{
	struct phase *phase;
	char *name;
	LDBLE p_read;
	LDBLE moles;
	LDBLE initial_moles;
};
int count_gas_phase;
int max_gas_phase;
struct gas_phase *gas_phase;
/*----------------------------------------------------------------------
 *   Solid solution
 *---------------------------------------------------------------------- */
struct s_s_assemblage
{
	int n_user;
	int n_user_end;
	char *description;
	int new_def;
	/*    int type; */
	/*    int solution_equilibria; */
	/*    int n_solution; */
	int count_s_s;
	struct s_s *s_s;
};
struct s_s
{
	char *name;
	struct s_s_comp *comps;
	int count_comps;
	LDBLE total_moles;
	LDBLE dn;
	LDBLE a0, a1;
	LDBLE ag0, ag1;
	int s_s_in;
	int miscibility;
	int spinodal;
	LDBLE tk, xb1, xb2;
	int input_case;
	LDBLE p[4];
};
struct s_s_comp
{
	char *name;
	struct phase *phase;
	LDBLE initial_moles;
	LDBLE moles;
	LDBLE init_moles;
	LDBLE delta;
	LDBLE fraction_x;
	LDBLE log10_lambda;
	LDBLE log10_fraction_x;
	LDBLE dn, dnc, dnb;
};
int count_s_s_assemblage;
int max_s_s_assemblage;
struct s_s_assemblage *s_s_assemblage;
/*----------------------------------------------------------------------
 *   Pure-phase assemblage
 *---------------------------------------------------------------------- */
struct pp_assemblage
{
	int n_user;
	int n_user_end;
	char *description;
	int new_def;
	struct elt_list *next_elt;
	int count_comps;
	struct pure_phase *pure_phases;
};
struct pure_phase
{
	struct phase *phase;
	char *name;
	char *add_formula;
	LDBLE si;
	LDBLE moles;
	LDBLE delta;
	LDBLE initial_moles;
	int force_equality;
	int dissolve_only;
	int precipitate_only;
};
int count_pp_assemblage;
int max_pp_assemblage;
struct pp_assemblage *pp_assemblage;
struct pp_assemblage *dbg_pp_assemblage;
/*----------------------------------------------------------------------
 *   Species_list
 *---------------------------------------------------------------------- */
struct species_list
{
	struct species *master_s;
	struct species *s;
	LDBLE coef;
};
int count_species_list;
int max_species_list;
struct species_list *species_list;
/*----------------------------------------------------------------------
 *   Jacobian and Mass balance lists
 *---------------------------------------------------------------------- */
struct list0
{
	LDBLE *target;
	LDBLE coef;
};
int count_sum_jacob0;	/* number of elements in sum_jacob0 */
int max_sum_jacob0;	/* calculated maximum number of elements in sum_jacob0 */
struct list0 *sum_jacob0;	/* array of pointers to targets and coefficients for array */

struct list1
{
	LDBLE *source;
	LDBLE *target;
};
int count_sum_mb1;		/* number of elements in sum_mb1 */
int max_sum_mb1;		/* calculated maximum number of elements in sum_mb1 */
struct list1 *sum_mb1;	/* array of pointers to sources and targets for mass
								   balance summations with coef = 1.0 */
int count_sum_jacob1;	/* number of elements in sum_jacob1 */
int max_sum_jacob1;	/* calculated maximum number of elements in sum_jacob1 */
struct list1 *sum_jacob1;	/* array of pointers to sources and targets for array
									   equations with coef = 1.0 */
struct list2
{
	LDBLE *source;
	LDBLE *target;
	LDBLE coef;
};
int count_sum_mb2;		/* number of elements in sum_mb2 */
int max_sum_mb2;		/* calculated maximum number of elements in sum_mb2 */
struct list2 *sum_mb2;	/* array of coefficients and pointers to sources and
								   targets for mass balance summations with coef != 1.0 */
int count_sum_jacob2;	/* number of elements in sum_jacob2 */
int max_sum_jacob2;	/* calculated maximum number of elements in sum_jacob2 */
struct list2 *sum_jacob2;	/* array of coefficients and pointers to sources and
									   targets, coef != 1.0 */
int count_sum_delta;	/* number of elements in sum_delta */
int max_sum_delta;		/* calculated maximum number of elements in sum_delta */
struct list2 *sum_delta;	/* array of pointers to sources, targets and coefficients for
									   summing deltas for mass balance equations */
/*----------------------------------------------------------------------
 *   Solution
 *---------------------------------------------------------------------- */
struct solution
{
	int new_def;
	int n_user;
	int n_user_end;
	char *description;
	LDBLE tc;
	LDBLE ph;
	LDBLE solution_pe;
	LDBLE mu;
	LDBLE ah2o;
	LDBLE density;
	LDBLE total_h;
	LDBLE total_o;
	LDBLE cb;
	LDBLE mass_water;
	LDBLE total_alkalinity;
	char *units;
	struct pe_data *pe;
	int default_pe;
	struct conc *totals;
	struct master_activity *master_activity;
	int count_master_activity;
	int count_isotopes;
	struct isotope *isotopes;
	struct master_activity *species_gamma;
	int count_species_gamma;
};
struct master_activity
{
	char *description;
	LDBLE la;
};
struct conc
{
	char *description;
	/*int skip; */
	LDBLE moles;
	LDBLE input_conc;
	char *units;
	char *equation_name;
	struct phase *phase;
	LDBLE phase_si;
	int n_pe;
	char *as;
	LDBLE gfw;
};
struct pe_data
{
	char *name;
	struct reaction *rxn;
};
struct isotope
{
	LDBLE isotope_number;
	char *elt_name;
	char *isotope_name;
	LDBLE total;
	LDBLE ratio;
	LDBLE ratio_uncertainty;
	LDBLE x_ratio_uncertainty;
	struct master *master;
	struct master *primary;
	LDBLE coef;					/* coefficient of element in phase */
};
struct solution **solution;
struct solution **dbg_solution;
int count_solution;
int max_solution;
struct iso
{
	char *name;
	LDBLE value;
	LDBLE uncertainty;
};
#ifdef SKIP
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
#else
struct iso iso_defaults[];
int count_iso_defaults;
#endif
#endif
struct iso *iso_defaults;
int count_iso_defaults;

/*----------------------------------------------------------------------
 *   Global solution
 *---------------------------------------------------------------------- */
char *title_x;
int new_x;
char *description_x;
LDBLE tc_x;
LDBLE tk_x;
LDBLE ph_x;
LDBLE solution_pe_x;
LDBLE mu_x;
LDBLE ah2o_x;
LDBLE density_x;
LDBLE total_h_x;
LDBLE total_o_x;
LDBLE cb_x;
LDBLE total_ions_x;
LDBLE mass_water_aq_x;
LDBLE mass_water_surfaces_x;
LDBLE mass_water_bulk_x;
char *units_x;
struct pe_data *pe_x;
int count_isotopes_x;
struct isotope *isotopes_x;
int default_pe_x;
/*int diffuse_layer_x;*/
enum DIFFUSE_LAYER_TYPE dl_type_x;
LDBLE total_carbon;
LDBLE total_co2;
LDBLE total_alkalinity;
LDBLE gfw_water;
LDBLE step_x;
LDBLE kin_time_x;
/*----------------------------------------------------------------------
 *   Transport data
 *---------------------------------------------------------------------- */
int count_cells;
int count_shifts;
int ishift;
int bcon_first;
int bcon_last;
int correct_disp;
LDBLE tempr;
LDBLE timest;
int simul_tr;
LDBLE diffc;
LDBLE heat_diffc;
int cell;
LDBLE mcd_substeps;
struct stag_data
{
	int count_stag;
	LDBLE exch_f;
	LDBLE th_m;
	LDBLE th_im;
} *stag_data;
int print_modulus;
int punch_modulus;
int dump_in;
int dump_modulus;
int transport_warnings;
struct cell_data
{
	LDBLE length;
	LDBLE mid_cell_x;
	LDBLE disp;
	LDBLE temp;
	LDBLE por;					/* free (uncharged) porewater porosities */
	LDBLE por_il;				/* interlayer water porosities */
	int punch;
	int print;
} *cell_data;
int multi_Dflag;		/* signals calc'n of multicomponent diffusion */
int interlayer_Dflag;	/* multicomponent diffusion and diffusion through interlayer porosity */
LDBLE default_Dw;		/* default species diffusion coefficient in water at 25oC, m2/s */
LDBLE multi_Dpor;		/* uniform porosity of free porewater in solid medium */
LDBLE interlayer_Dpor;	/* uniform porosity of interlayer space of montmorillonite in solid medium */
LDBLE multi_Dpor_lim;	/* limiting free porewater porosity where transport stops */
LDBLE interlayer_Dpor_lim;	/* limiting interlayer porosity where transport stops */
LDBLE multi_Dn;		/* exponent to calculate pore water diffusion coefficient,
								   Dp = Dw * (multi_Dpor)^multi_Dn */
LDBLE interlayer_tortf;	/* tortuosity_factor in interlayer porosity,
									   Dpil = Dw / interlayer_tortf */

int cell_no;
/*----------------------------------------------------------------------
 *   Advection data
 *---------------------------------------------------------------------- */
int count_ad_cells;
int count_ad_shifts;
int print_ad_modulus;
int punch_ad_modulus;
int *advection_punch, *advection_print;
LDBLE advection_kin_time;
LDBLE advection_kin_time_defined;
int advection_warnings;
/*----------------------------------------------------------------------
 *   Keywords
 *---------------------------------------------------------------------- */
struct key
{
	char *name;
	int keycount;
};
struct const_key
{
	const char *name;
	int keycount;
};

	struct const_key *keyword;
	int NKEYS;

struct key *keyword_hash;
int new_model, new_exchange, new_pp_assemblage, new_surface,
	new_reaction, new_temperature, new_mix, new_solution, new_gas_phase,
	new_inverse, new_punch, new_s_s_assemblage, new_kinetics, new_copy,
	new_pitzer;

/*----------------------------------------------------------------------
 *   Elements
 *---------------------------------------------------------------------- */
struct element
{
	char *name;					/* element name */
	/*    int in; */
	struct master *master;
	struct master *primary;
	LDBLE gfw;
};
struct element **elements;
int count_elements;
int max_elements;
struct element *element_h_one;

/*----------------------------------------------------------------------
 *   Element List
 *---------------------------------------------------------------------- */
struct elt_list
{								/* list of name and number of elements in an equation */
	struct element *elt;		/* pointer to element structure */
	LDBLE coef;					/* number of element e's in eqn */
};
struct elt_list *elt_list;	/* structure array of working space while reading equations
									   names are in "strings", initially in input order */
int count_elts;		/* number of elements in elt_list = position of next */
int max_elts;
/*----------------------------------------------------------------------
 *   Reaction
 *---------------------------------------------------------------------- */
struct reaction
{
	LDBLE logk[8];
	LDBLE dz[3];
	struct rxn_token *token;
};
struct rxn_token
{
	struct species *s;
	LDBLE coef;
	char *name;
};
/*----------------------------------------------------------------------
 *   Species
 *---------------------------------------------------------------------- */
struct species
{								/* all data pertinent to an aqueous species */
	char *name;					/* name of species */
	char *mole_balance;			/* formula for mole balance */
	int in;						/* species used in model if TRUE */
	int number;
	struct master *primary;		/* points to master species list, NULL if not primary master */
	struct master *secondary;	/* points to master species list, NULL if not secondary master */
	LDBLE gfw;					/* gram formula wt of species */
	LDBLE z;					/* charge of species */
	LDBLE dw;					/* tracer diffusion coefficient in water at 25oC, m2/s */
	LDBLE erm_ddl;				/* enrichment factor in DDL */
	LDBLE equiv;				/* equivalents in exchange species */
	LDBLE alk;					/* alkalinity of species, used for cec in exchange */
	LDBLE carbon;				/* stoichiometric coefficient of carbon in species */
	LDBLE co2;					/* stoichiometric coefficient of C(4) in species */
	LDBLE h;					/* stoichiometric coefficient of H in species */
	LDBLE o;					/* stoichiometric coefficient of O in species */
	LDBLE dha, dhb, a_f;		/* WATEQ Debye Huckel a and b-dot; active_fraction coef for exchange species */
	LDBLE lk;					/* log10 k at working temperature */
	LDBLE logk[8];				/* log kt0, delh, 6 coefficients analalytical expression */
/* VP: Density Start */
	LDBLE millero[6];		    /* regression coefficients to calculate temperature dependent phi_0 and b_v of Millero density model */
/* VP: Density End */
	DELTA_H_UNIT original_units;	/* enum with original delta H units */
	int count_add_logk;
	struct name_coef *add_logk;
	LDBLE lg;					/* log10 activity coefficient, gamma */
	LDBLE lg_pitzer;			/* log10 activity coefficient, from pitzer calculation */
	LDBLE lm;					/* log10 molality */
	LDBLE la;					/* log10 activity */
	LDBLE dg;					/* gamma term for jacobian */
	LDBLE dg_total_g;
	LDBLE moles;				/* moles in solution; moles/mass_water = molality */
	int type;					/* flag indicating presence in model and types of equations */
	int gflag;					/* flag for preferred activity coef eqn */
	int exch_gflag;				/* flag for preferred activity coef eqn */
	struct elt_list *next_elt;	/* pointer to next element */
	struct elt_list *next_secondary;
	struct elt_list *next_sys_total;
	int check_equation;			/* switch to check equation for charge and element balance */
	struct reaction *rxn;		/* pointer to data base reaction */
	struct reaction *rxn_s;		/* pointer to reaction converted to secondary and primary
								   master species */
	struct reaction *rxn_x;		/* reaction to be used in model */
	LDBLE tot_g_moles;			/* (1 + sum(g)) * moles */
	LDBLE tot_dh2o_moles;		/* sum(moles*g*Ws/Waq) */
	struct species_diff_layer *diff_layer;	/* information related to diffuse layer factors for each
											   surface */
	LDBLE cd_music[5];
	LDBLE dz[3];
};
struct logk
{								/* Named log K's */
	char *name;					/* name of species */
	LDBLE lk;					/* log10 k at working temperature */
	LDBLE log_k[8];				/* log kt0, delh, 6 coefficients analalytical expression */
	DELTA_H_UNIT original_units;	/* enum with original delta H units */
	int count_add_logk;
	int done;
	struct name_coef *add_logk;
	LDBLE log_k_original[8];	/* log kt0, delh, 5 coefficients analalytical expression */
};
struct logk **logk;
int count_logk;
int max_logk;
struct species_diff_layer
{
	struct surface_charge *charge;
	int count_g;
	LDBLE g_moles;
	LDBLE dg_g_moles;			/* g_moles*dgterm */
	LDBLE dx_moles;
	LDBLE dh2o_moles;			/* moles*g*Ws/Waq */
	LDBLE drelated_moles;		/* for related phase */
};
char *moles_per_kilogram_string;
char *pe_string;

struct species **s;
int count_s;
int max_s;

struct species **s_x;
int count_s_x;
int max_s_x;

struct species *s_h2o;
struct species *s_hplus;
struct species *s_h3oplus;
struct species *s_eminus;
struct species *s_co3;
struct species *s_h2;
struct species *s_o2;
/*----------------------------------------------------------------------
 *   Phases
 *---------------------------------------------------------------------- */
struct phase
{								/* all data pertinent to a pure solid phase */
	char *name;					/* name of species */
	char *formula;				/* chemical formula */
	int in;						/* species used in model if TRUE */
	LDBLE lk;					/* log10 k at working temperature */
	LDBLE logk[8];				/* log kt0, delh, 6 coefficients analalytical expression */
	DELTA_H_UNIT original_units;	/* enum with original delta H units */
	int count_add_logk;
	struct name_coef *add_logk;
	LDBLE moles_x;
	LDBLE p_soln_x;
	LDBLE fraction_x;
	LDBLE log10_lambda, log10_fraction_x;
	LDBLE dn, dnb, dnc;
	LDBLE gn, gntot;
	LDBLE gn_n, gntot_n;

	int type;					/* flag indicating presence in model and types of equations */
	struct elt_list *next_elt;	/* pointer to list of elements in phase */
	struct elt_list *next_sys_total;
	int check_equation;			/* switch to check equation for charge and element balance */
	struct reaction *rxn;		/* pointer to data base reaction */
	struct reaction *rxn_s;		/* pointer to reaction converted to secondary and primary
								   master species */
	struct reaction *rxn_x;		/* reaction to be used in model */
	int in_system;
};
struct phase **phases;
int count_phases;
int max_phases;
/*----------------------------------------------------------------------
 *   Master species
 *---------------------------------------------------------------------- */
struct master
{								/* list of name and number of elements in an equation */
	int in;						/* TRUE if in model, FALSE if out, REWRITE if other mb eq */
	int number;					/* sequence number in list of masters */
	int last_model;				/* saved to determine if model has changed */
	int type;					/* AQ or EX */
	int primary;				/* TRUE if master species is primary */
	LDBLE coef;					/* coefficient of element in master species */
	LDBLE total;				/* total concentration for element or valence state */
	LDBLE isotope_ratio;
	LDBLE isotope_ratio_uncertainty;
	int isotope;
	LDBLE total_primary;
	/*    LDBLE la;  *//* initial guess of master species log activity */
	struct element *elt;		/* element structure */
	LDBLE alk;					/* alkalinity of species */
	LDBLE gfw;					/* default gfw for species */
	char *gfw_formula;			/* formula from which to calcuate gfw */
	struct unknown *unknown;	/* pointer to unknown structure */
	struct species *s;			/* pointer to species structure */
	struct reaction *rxn_primary;	/* reaction writes master species in terms of primary
									   master species */
	struct reaction *rxn_secondary;	/* reaction writes master species in terms of secondary
									   master species */
	struct reaction **pe_rxn;	/* e- written in terms of redox couple (or e-), points
								   to location */
	int minor_isotope;
};
struct master **master;	/* structure array of master species */
struct master **dbg_master;
int count_master;
int max_master;
/*----------------------------------------------------------------------
 *   Unknowns
 *---------------------------------------------------------------------- */
struct unknown
{
	int type;
	LDBLE moles;
	LDBLE ln_moles;
	LDBLE f;
	LDBLE sum;
	LDBLE delta;
	LDBLE la;
	int number;
	char *description;
	struct master **master;
	struct phase *phase;
	LDBLE si;
	struct gas_phase *gas_phase;
	struct conc *total;
	struct species *s;
	struct exch_comp *exch_comp;
	struct pure_phase *pure_phase;
	struct s_s *s_s;
	struct s_s_comp *s_s_comp;
	int s_s_comp_number;
	int s_s_in;
	struct surface_comp *surface_comp;
	LDBLE related_moles;
	struct unknown *potential_unknown, *potential_unknown1,
		*potential_unknown2;
	int count_comp_unknowns;
	struct unknown **comp_unknowns;	/* list for CD_MUSIC of comps that contribute to 0 plane mass-balance term */
	struct unknown *phase_unknown;
	struct surface_charge *surface_charge;
	LDBLE mass_water;
	int dissolve_only;
	LDBLE inert_moles;
};
struct unknown **x;
int count_unknowns;
int max_unknowns;

struct unknown *ah2o_unknown;
struct unknown *alkalinity_unknown;
struct unknown *carbon_unknown;
struct unknown *charge_balance_unknown;
struct unknown *exchange_unknown;
struct unknown *mass_hydrogen_unknown;
struct unknown *mass_oxygen_unknown;
struct unknown *mb_unknown;
struct unknown *mu_unknown;
struct unknown *pe_unknown;
struct unknown *ph_unknown;
struct unknown *pure_phase_unknown;
struct unknown *solution_phase_boundary_unknown;
struct unknown *surface_unknown;
struct unknown *gas_unknown;
struct unknown *s_s_unknown;
/*----------------------------------------------------------------------
 *   Reaction work space
 *---------------------------------------------------------------------- */
struct reaction_temp
{
	LDBLE logk[8];
	LDBLE dz[3];
	struct rxn_token_temp *token;
};
struct rxn_token_temp
{								/* data for equations, aq. species or minerals */
	char *name;					/* pointer to a species name (formula) */
	LDBLE z;					/* charge on species */
	struct species *s;
	struct unknown *unknown;
	LDBLE coef;					/* coefficient of species name */
};
struct reaction_temp trxn;	/* structure array of working space while reading equations
									   species names are in "temp_strings" */
int count_trxn;		/* number of reactants in trxn = position of next */
int max_trxn;
struct unknown_list
{
	struct unknown *unknown;
	LDBLE *source;
	LDBLE *gamma_source;
	/*    int row; */
	/*    int col; */
	LDBLE coef;
};
struct unknown_list *mb_unknowns;
int count_mb_unknowns;
int max_mb_unknowns;
/* ----------------------------------------------------------------------
 *   Print
 * ---------------------------------------------------------------------- */
struct prints
{
	int all;
	int initial_solutions;
	int initial_exchangers;
	int reactions;
	int gas_phase;
	int s_s_assemblage;
	int pp_assemblage;
	int surface;
	int exchange;
	int kinetics;
	int totals;
	int eh;
	int species;
	int saturation_indices;
	int irrev;
	int mix;
	int reaction;
	int use;
	int logfile;
	int punch;
	int status;
	int inverse;
	int dump;
	int user_print;
	int headings;
	int user_graph;
	int echo_input;
	int warnings;
	int initial_isotopes;
	int isotope_ratios;
	int isotope_alphas;
	int hdf;
	int alkalinity;
};
struct prints pr;
int status_on;
int count_warnings;

/* ----------------------------------------------------------------------
 *   RATES
 * ---------------------------------------------------------------------- */
struct rate
{
	char *name;
	char *commands;
	int new_def;
	void *linebase;
	void *varbase;
	void *loopbase;
};
struct rate *rates;
int count_rates;
LDBLE rate_m, rate_m0, *rate_p, rate_time, rate_sim_time_start,
	rate_sim_time_end, rate_sim_time, rate_moles, initial_total_time;
int count_rate_p;
/* ----------------------------------------------------------------------
 *   USER PRINT COMMANDS
 * ---------------------------------------------------------------------- */
struct rate *user_print;
struct rate *user_punch;
char **user_punch_headings;
int user_punch_count_headings;
#ifdef PHREEQ98
struct rate *user_graph;
char **user_graph_headings;
int user_graph_count_headings;
#endif

/* ----------------------------------------------------------------------
 *   GLOBAL DECLARATIONS
 * ---------------------------------------------------------------------- */
char error_string[10 * MAX_LENGTH];
int simulation;
int state;
int reaction_step;
int transport_step;
int transport_start;
int advection_step;
int stop_program;
int incremental_reactions;

int count_strings;
int max_strings;

LDBLE *array;
LDBLE *delta;
LDBLE *residual;

int input_error;
int next_keyword;
int parse_error;
int paren_count;
int iterations;
int gamma_iterations;
int run_reactions_iterations;

int max_line;
char *line;
char *line_save;

LDBLE LOG_10;

int debug_model;
int debug_prep;
int debug_set;
int debug_diffuse_layer;
int debug_inverse;

LDBLE inv_tol_default;
int itmax;
LDBLE ineq_tol;
LDBLE convergence_tolerance;
LDBLE step_size;
LDBLE pe_step_size;
LDBLE step_size_now;
LDBLE pe_step_size_now;
LDBLE pp_scale;
LDBLE pp_column_scale;
int diagonal_scale;	/* 0 not used, 1 used */
int mass_water_switch;
int delay_mass_water;
LDBLE censor;
int aqueous_only;
int negative_concentrations;
int calculating_deriv;
int numerical_deriv;
int count_total_steps;

int phast;
LDBLE *llnl_temp, *llnl_adh, *llnl_bdh, *llnl_bdot, *llnl_co2_coefs;
int llnl_count_temp, llnl_count_adh, llnl_count_bdh, llnl_count_bdot,
	llnl_count_co2_coefs;

char *selected_output_file_name;
char *dump_file_name;
struct spread_row
{
	int count;
	int empty, string, number;
	char **char_vector;
	LDBLE *d_vector;
	int *type_vector;
};
struct defaults
{
	LDBLE temp;
	LDBLE density;
	char *units;
	char *redox;
	LDBLE ph;
	LDBLE pe;
	LDBLE water;
	int count_iso;
	struct iso *iso;
};
struct spread_sheet
{
	struct spread_row *heading;
	struct spread_row *units;
	int count_rows;
	struct spread_row **rows;
	struct defaults defaults;
};
#ifdef PHREEQCI_GUI
struct spread_sheet g_spread_sheet;
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

HashTable *strings_hash_table;
HashTable *elements_hash_table;
HashTable *species_hash_table;
HashTable *phases_hash_table;
HashTable *keyword_hash_table;
HashTable *logk_hash_table;
HashTable *master_isotope_hash_table;

#if defined(PHREEQCI_GUI)
#include "../../phreeqci_gui.h"
#endif /* defined(PHREEQCI_GUI) */

struct name_coef match_tokens[50];
int count_match_tokens;
struct master_isotope
{
	char *name;
	struct master *master;
	struct element *elt;
	char *units;
	LDBLE standard;
	LDBLE ratio;
	LDBLE moles;
	int total_is_major;
	int minor_isotope;
};
int count_master_isotope;
struct master_isotope **master_isotope;
int max_master_isotope;
int initial_solution_isotopes;

#define OPTION_EOF -1
#define OPTION_KEYWORD -2
#define OPTION_ERROR -3
#define OPTION_DEFAULT -4
#define OPT_1 -5

struct calculate_value
{
	char *name;
	LDBLE value;
	char *commands;
	int new_def;
	int calculated;
	void *linebase;
	void *varbase;
	void *loopbase;
};
int count_calculate_value;
struct calculate_value **calculate_value;
int max_calculate_value;
HashTable *calculate_value_hash_table;

struct isotope_ratio
{
	char *name;
	char *isotope_name;
	LDBLE ratio;
	LDBLE converted_ratio;
};
int count_isotope_ratio;
struct isotope_ratio **isotope_ratio;
int max_isotope_ratio;
HashTable *isotope_ratio_hash_table;
struct isotope_alpha
{
	char *name;
	char *named_logk;
	LDBLE value;
};
int count_isotope_alpha;
struct isotope_alpha **isotope_alpha;
int max_isotope_alpha;
HashTable *isotope_alpha_hash_table;

int phreeqc_mpi_myself;

enum entity_type
{ Solution, Reaction, Exchange, Surface, Gas_phase, Pure_phase, Ss_phase,
	Kinetics, Mix, Temperature, UnKnown
};

int first_read_input;
char *user_database;
int pitzer_model, sit_model, pitzer_pe;
int full_pitzer, always_full_pitzer, ICON, IC;
LDBLE COSMOT;
LDBLE AW;
int have_punch_name;
/* VP: Density Start */
int print_density;
/* VP: Density End */

jmp_buf mark;
LDBLE *zeros;
int zeros_max;

#if defined(WIN32_MEMORY_DEBUG)
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#endif

struct system
{
	struct solution *solution;
	struct exchange *exchange;
	struct pp_assemblage *pp_assemblage;
	struct gas_phase *gas_phase;
	struct s_s_assemblage *s_s_assemblage;
	struct kinetics *kinetics;
	struct surface *surface;
};

LDBLE pore_volume;

//***************************** end of global.h *****************************

// basic.c -------------------------------
//static char const svnid[] = "$Id: basic.c 3490 2009-05-13 17:00:08Z dlpark $";


#ifdef PHREEQ98
void GridChar(char *s, char *a);
extern int colnr, rownr;
#endif

int n_user_punch_index;
int sget_logical_line(char **ptr, int *l, char *return_line);

#define checking	true
#define varnamelen      20
#define maxdims	 4

typedef Char varnamestring[varnamelen + 1];
typedef Char string255[256];

typedef LDBLE numarray[];
typedef Char *strarray[];

#define forloop	 0
#define whileloop       1
#define gosubloop       2
#define tokvar	  0
#define toknum	  1
#define tokstr	  2
#define toksnerr	3
#define tokplus	 4
#define tokminus	5
#define toktimes	6
#define tokdiv	  7
#define tokup	   8
#define toklp	   9
#define tokrp	   10
#define tokcomma	11
#define toksemi	 12
#define tokcolon	13
#define tokeq	   14
#define toklt	   15
#define tokgt	   16
#define tokle	   17
#define tokge	   18
#define tokne	   19
#define tokand	  20
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

int free_dim_stringvar(struct varrec *varbase);
long EXCP_LINE;
void exec(void);

/*$if not checking$
   $range off$
$end$*/
HashTable *command_hash_table;

struct const_key *command;
int NCMDS;



int basic_renumber(char *commands, void **lnbase, void **vbase, void **lpbase);
void restoredata(void);
void clearloops(void);
void clearvar(varrec * v);
void clearvars(void);
Char * numtostr(Char * Result, LDBLE n);
void parse(Char * inbuf, tokenrec ** buf);

void listtokens(FILE * f, tokenrec * buf);
void disposetokens(tokenrec ** tok);
void parseinput(tokenrec ** buf);
void errormsg(const Char * s);
void snerr(void);
void tmerr(void);
void badsubscr(void);
LDBLE realfactor(struct LOC_exec *LINK);
Char * strfactor(struct LOC_exec * LINK);
Char *stringfactor(Char * Result, struct LOC_exec *LINK);
long intfactor(struct LOC_exec *LINK);
LDBLE realexpr(struct LOC_exec *LINK);
Char * strexpr(struct LOC_exec * LINK);
Char * stringexpr(Char * Result, struct LOC_exec * LINK);
long intexpr(struct LOC_exec *LINK);
void require(int k, struct LOC_exec *LINK);
void skipparen(struct LOC_exec *LINK);
varrec * findvar(struct LOC_exec *LINK);
valrec factor(struct LOC_exec *LINK);
valrec upexpr(struct LOC_exec * LINK);
valrec term(struct LOC_exec * LINK);
valrec sexpr(struct LOC_exec * LINK);
valrec relexpr(struct LOC_exec * LINK);
valrec andexpr(struct LOC_exec * LINK);
valrec expr(struct LOC_exec *LINK);

void checkextra(struct LOC_exec *LINK);
boolean iseos(struct LOC_exec *LINK);
void skiptoeos(struct LOC_exec *LINK);
linerec * findline(long n);
linerec * mustfindline(long n);
void cmdend(struct LOC_exec *LINK);
void cmdnew(struct LOC_exec *LINK);
void cmdlist(struct LOC_exec *LINK);
void cmdload(boolean merging, Char * name, struct LOC_exec *LINK);
void cmdrun(struct LOC_exec *LINK);
void cmdsave(struct LOC_exec *LINK);
void cmdput(struct LOC_exec *LINK);
void cmdchange_por(struct LOC_exec *LINK);
void cmdchange_surf(struct LOC_exec *LINK);
void cmdbye(void);
void cmddel(struct LOC_exec *LINK);
void cmdrenum(struct LOC_exec *LINK);
void cmdprint(struct LOC_exec *LINK);
void cmdpunch(struct LOC_exec *LINK);
#ifdef PHREEQ98
void cmdgraph_x(struct LOC_exec *LINK);
void cmdgraph_y(struct LOC_exec *LINK);
void cmdgraph_sy(struct LOC_exec *LINK);
#endif
void cmdlet(boolean implied, struct LOC_exec *LINK);
void cmdgoto(struct LOC_exec *LINK);
void cmdif(struct LOC_exec *LINK);
void cmdelse(struct LOC_exec *LINK);
boolean skiploop(int up, int dn, struct LOC_exec *LINK);
void cmdfor(struct LOC_exec *LINK);
void cmdnext(struct LOC_exec *LINK);
void cmdwhile(struct LOC_exec *LINK);
void cmdwend(struct LOC_exec *LINK);
void cmdgosub(struct LOC_exec *LINK);
void cmdreturn(struct LOC_exec *LINK);
void cmdread(struct LOC_exec *LINK);
void cmddata(struct LOC_exec *LINK);
void cmdrestore(struct LOC_exec *LINK);
void cmdgotoxy(struct LOC_exec *LINK);
void cmdon(struct LOC_exec *LINK);
void cmddim(struct LOC_exec *LINK);
void cmdpoke(struct LOC_exec *LINK);


// basicsubs.c -------------------------------
struct system_species
{
	char *name;
	char *type;
	LDBLE moles;
};
struct system_species *sys;
int count_sys, max_sys;

int system_total_solids(struct exchange *exchange_ptr,
					struct pp_assemblage *pp_assemblage_ptr,
					struct gas_phase *gas_phase_ptr,
					struct s_s_assemblage *s_s_assemblage_ptr,
					struct surface *surface_ptr);

LDBLE sys_tot;
LDBLE AA_basic, BB_basic, CC, I_m, rho_0;
LDBLE solution_mass, solution_volume;
static LDBLE f_rho(LDBLE rho_old, void *cookie);

// cl1.c -------------------------------
void cl1_space(int check, int n2d, int klm, int nklmd);
//void zero_double(LDBLE * target, int n);
//LDBLE *x_arg = NULL, *res_arg = NULL, *scratch = NULL;
//int x_arg_max = 0, res_arg_max = 0, scratch_max = 0;
LDBLE *x_arg, *res_arg, *scratch;
int x_arg_max, res_arg_max, scratch_max;

// cl1mp.c -------------------------------
int cl1mp(int k, int l, int m, int n,
		  int nklmd, int n2d,
		  LDBLE * q_arg,
		  int *kode, LDBLE toler,
		  int *iter, LDBLE * x_arg, LDBLE * res_arg, LDBLE * error,
		  LDBLE * cu_arg, int *iu, int *s, int check, LDBLE censor_arg);

// class_main.c -------------------------------
#ifdef DOS
int write_banner(void);
#endif

// dw.c -------------------------------
int BB(LDBLE T);
LDBLE PS(LDBLE T);
LDBLE VLEST(LDBLE T);
int DFIND(LDBLE * DOUT, LDBLE P, LDBLE D, LDBLE T);
int QQ(LDBLE T, LDBLE D);
LDBLE BASE(LDBLE D);

/* COMMON /QQQQ/ */
LDBLE Q0, Q5;
/* COMMON /ACONST/ */
//LDBLE GASCON = 0.461522e0, TZ = 647.073e0, AA = 1.e0;
LDBLE GASCON, TZ, AA;
LDBLE Z, DZ, Y;
/* COMMON /ELLCON/ */
//LDBLE G1 = 11.e0, G2 = 44.333333333333e0, GF = 3.5e0;
LDBLE G1, G2, GF;
LDBLE B1, B2, B1T, B2T, B1TT, B2TT;

// input.c -------------------------------

int reading_database(void);
struct read_callback s_read_callback;
int check_line(const char *string, int allow_empty, int allow_eof,
		   int allow_keyword, int print);
// integrate.c -------------------------------
LDBLE g_function(LDBLE x_value);
LDBLE midpnt(LDBLE x1, LDBLE x2, int n);
void polint(LDBLE * xa, LDBLE * ya, int n, LDBLE xv, LDBLE * yv,
				   LDBLE * dy);
LDBLE qromb_midpnt(LDBLE x1, LDBLE x2);
LDBLE z, xd, alpha;
struct surface_charge *surface_charge_ptr;

LDBLE calc_psi_avg(LDBLE surf_chrg_eq);
int calc_all_donnan_music(void);
int calc_init_donnan_music(void);
//integrate.c
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

// inverse.c -------------------------------

int add_to_file(const char *filename, char *string);
int bit_print(unsigned long bits, int l);
int carbon_derivs(struct inverse *inv_ptr);
int check_isotopes(struct inverse *inv_ptr);
int check_solns(struct inverse *inv_ptr);
int count_isotope_unknowns(struct inverse *inv_ptr,
								  struct isotope **isotope_unknowns);
struct isotope *get_isotope(struct solution *solution_ptr, const char *elt);
struct conc *get_inv_total(struct solution *solution_ptr, const char *elt);
int isotope_balance_equation(struct inverse *inv_ptr, int row, int n);
int post_mortem(void);
unsigned long get_bits(unsigned long bits, int position, int number);
unsigned long minimal_solve(struct inverse *inv_ptr,
								   unsigned long minimal_bits);
void dump_netpath(struct inverse *inv_ptr);
int dump_netpath_pat(struct inverse *inv_ptr);
int next_set_phases(struct inverse *inv_ptr, int first_of_model_size,
						   int model_size);
int phase_isotope_inequalities(struct inverse *inv_ptr);
int print_model(struct inverse *inv_ptr);
int punch_model_heading(struct inverse *inv_ptr);
int punch_model(struct inverse *inv_ptr);
void print_isotope(FILE * netpath_file, struct solution *solution_ptr,
				   const char *elt, const char *string);
void print_total(FILE * netpath_file, struct solution *solution_ptr,
				 const char *elt, const char *string);
void print_total_multi(FILE * netpath_file, struct solution *solution_ptr,
					   const char *string, const char *elt0,
					   const char *elt1, const char *elt2, const char *elt3,
					   const char *elt4);
void print_total_pat(FILE * netpath_file, const char *elt,
					 const char *string);
int range(struct inverse *inv_ptr, unsigned long cur_bits);
int save_bad(unsigned long bits);
int save_good(unsigned long bits);
int save_minimal(unsigned long bits);
unsigned long set_bit(unsigned long bits, int position, int value);
int setup_inverse(struct inverse *inv_ptr);
int set_initial_solution(int n_user_old, int n_user_new);
int set_ph_c(struct inverse *inv_ptr,
					int i,
					struct solution *soln_ptr_orig,
					int n_user_new,
					LDBLE d_alk, LDBLE ph_factor, LDBLE alk_factor);
int shrink(struct inverse *inv_ptr, LDBLE * array_in,
				  LDBLE * array_out, int *k, int *l, int *m, int *n,
				  unsigned long cur_bits, LDBLE * delta_l, int *col_back_l,
				  int *row_back_l);
int solve_inverse(struct inverse *inv_ptr);
int solve_with_mask(struct inverse *inv_ptr, unsigned long cur_bits);
int subset_bad(unsigned long bits);
int subset_minimal(unsigned long bits);
int superset_minimal(unsigned long bits);
int write_optimize_names(struct inverse *inv_ptr);
FILE *netpath_file;
int count_inverse_models, count_pat_solutions;

// isotopes.c -------------------------------
int calculate_value_init(struct calculate_value *calculate_value_ptr);
int isotope_alpha_init(struct isotope_alpha *isotope_alpha_ptr);
int isotope_ratio_init(struct isotope_ratio *isotope_ratio_ptr);
int master_isotope_init(struct master_isotope *master_isotope_ptr);

//#ifdef SKIP_SOME
// kinetics.c -------------------------------
//#define KINETICS_EXTERNAL 
//#include "sundialstypes.h"		/* definitions of types realtype and                        */
//							 /* integertype, and the constant FALSE            */
//#include "cvode.h"				/* prototypes for CVodeMalloc, CVode, and            */
//							 /* CVodeFree, constants OPT_SIZE, BDF, NEWTON,   */
//							 /* SV, SUCCESS, NST,NFE,NSETUPS, NNI, NCFN, NETF */
//#include "cvdense.h"			/* prototype for CVDense, constant DENSE_NJE    */
//#include "nvector_serial.h"		/* definitions of type N_Vector and macro          */
//							 /* NV_Ith_S, prototypes for N_VNew, N_VFree      */
//#include "dense.h"				/* definitions of type DenseMat, macro DENSE_ELEM */
//#include "nvector.h"
//#include "kinetics.h"
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
static void f(integertype N, realtype t, N_Vector y, N_Vector ydot,
			  void *f_data);

static void Jac(integertype N, DenseMat J, RhsFn f, void *f_data, realtype t,
				N_Vector y, N_Vector fy, N_Vector ewt, realtype h,
				realtype uround, void *jac_data, long int *nfePtr,
				N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

int calc_final_kinetic_reaction(struct kinetics *kinetics_ptr);
int calc_kinetic_reaction(struct kinetics *kinetics_ptr,
								 LDBLE time_step);
int rk_kinetics(int i, LDBLE kin_time, int use_mix, int nsaver,
					   LDBLE step_fraction);
int set_reaction(int i, int use_mix, int use_kinetics);
int set_transport(int i, int use_mix, int use_kinetics, int nsaver);
int store_get_equi_reactants(int k, int kin_end);

LDBLE *m_original;
LDBLE *m_temp;

//#endif /* SKIP_SOME */

// mainsubs  -------------------------------
int copy_use(int i);
int set_use(void);

// model.c -------------------------------
LDBLE min_value;
LDBLE s_s_root(LDBLE a0, LDBLE a1, LDBLE kc, LDBLE kb, LDBLE xcaq,
					  LDBLE xbaq);
LDBLE s_s_halve(LDBLE a0, LDBLE a1, LDBLE x0, LDBLE x1, LDBLE kc,
					   LDBLE kb, LDBLE xcaq, LDBLE xbaq);
LDBLE s_s_f(LDBLE xb, LDBLE a0, LDBLE a1, LDBLE kc, LDBLE kb,
				   LDBLE xcaq, LDBLE xbaq);
int numerical_jacobian(void);
void set_inert_moles(void);
void unset_inert_moles(void);

#ifdef SLNQ
int add_trivial_eqns(int rows, int cols, LDBLE * matrix);
int slnq(int n, LDBLE * a, LDBLE * delta, int ncols, int print);
#endif

int calc_gas_pressures(void);
int calc_s_s_fractions(void);
int gammas(LDBLE mu);
int initial_guesses(void);
int revise_guesses(void);
int s_s_binary(struct s_s *s_s_ptr);
int s_s_ideal(struct s_s *s_s_ptr);

int remove_unstable_phases;
int gas_in;
void ineq_init(int max_row_count, int max_column_count);

//LDBLE min_value = 1e-10;
//
//LDBLE *normal = NULL, *ineq_array = NULL, *res = NULL, *cu = NULL, *zero =
//	NULL, *delta1 = NULL;
//int *iu = NULL, *is = NULL, *back_eq = NULL;
//int normal_max = 0, ineq_array_max = 0, res_max = 0, cu_max = 0, zero_max =
//	0, delta1_max = 0, iu_max = 0, is_max = 0, back_eq_max = 0;

LDBLE model_min_value;
LDBLE *normal, *ineq_array, *inv_res, *cu, *zero, *delta1;
int *inv_iu, *inv_is, *back_eq;
int normal_max, ineq_array_max, res_max, cu_max, zero_max, 
	delta1_max, iu_max, is_max, back_eq_max;

// output.c -------------------------------
private:
  #include "output.h"
#define MAX_CALLBACKS 10
//struct output_callback output_callbacks[MAX_CALLBACKS];
//size_t count_output_callback = 0;
//int forward_output_to_log = 0;
struct output_callback *output_callbacks;
size_t count_output_callback;
int forward_output_to_log;
int output_message(const int type, const char *err_str, const int stop,
			   const char *format, va_list args);

// parse.c -------------------------------
int get_coef(LDBLE * coef, char **eqnaddr);
int get_secondary(char **t_ptr, char *element, int *i);
int get_species(char **ptr);

// phqalloc.c -------------------------------
#if !defined(NDEBUG)
void *PHRQ_malloc(size_t, const char *, int);
void *PHRQ_calloc(size_t, size_t, const char *, int);
void *PHRQ_realloc(void *, size_t, const char *, int);
#else
void *PHRQ_malloc(size_t);
void *PHRQ_calloc(size_t, size_t);
void *PHRQ_realloc(void *, size_t);
#endif
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
//PHRQMemHeader *s_pTail = NULL;
PHRQMemHeader *s_pTail;
void PHRQ_free(void *ptr);
void PHRQ_free_all(void);


// phreeqc_files.c -------------------------------

//const char *default_data_base = "phreeqc.dat";
char *default_data_base;

//FILE *input_file = NULL;
//FILE *database_file = NULL;
//FILE *output = NULL;		/* OUTPUT_MESSAGE */
//FILE *log_file = NULL;	/* OUTPUT_LOG */
//FILE *punch_file = NULL;	/* OUTPUT_PUNCH */
//FILE *error_file = NULL;	/* OUTPUT_ERROR */
//FILE *dump_file = NULL;	/* OUTPUT_DUMP */
FILE *input_file;
FILE *database_file;
FILE *output;		/* OUTPUT_MESSAGE */
FILE *log_file;	    /* OUTPUT_LOG */
FILE *punch_file;	/* OUTPUT_PUNCH */
FILE *error_file;	/* OUTPUT_ERROR */
FILE *dump_file;	/* OUTPUT_DUMP */

#ifdef PHREEQ98
int outputlinenr;
void check_line_breaks(char *s);
char *prefix_database_dir(char *s);
char *LogFileNameC;
void show_progress(const int type, char *s);
char progress_str[512];
#endif

int fileop_handler(const int type, int (*PFN) (FILE *));
int open_handler(const int type, const char *file_name);
int output_handler(const int type, const char *err_str,
						  const int stop, void *cookie, const char *format,
						  va_list args);
static int rewind_wrapper(FILE * file_ptr);

//#ifdef SKIP_SOME
// p2clib.c
//#define extern
//#include "p2c.h"
//#undef extern 
//int P_argc;
//char **P_argv;
char *_ShowEscape(char *buf, int code, int ior, char *prefix);

//int P_escapecode;
//int P_ioresult;

//long EXCP_LINE;					/* Used by Pascal workstation system */

//Anyptr __MallocTemp__;

//__p2c_jmp_buf *__top_jb;
Void PASCAL_MAIN(int argc, char **argv);
int _Escape(int code);
int _EscIO(int code);
int  _OutMem(void);
int  _CaseCheck(void);
int  _NilCheck(void);
int P_escapecode;
int P_ioresult;
__p2c_jmp_buf *__top_jb;
//#endif /* SKIP_SOME */

// pitzer.c -------------------------------
#define PITZER_EXTERNAL
#include "pitzer.h"

/* variables */
LDBLE A0;
struct species **spec, **cations, **anions, **neutrals;
int count_cations, count_anions, count_neutrals;
int MAXCATIONS, FIRSTANION, MAXNEUTRAL;
struct pitz_param *mcb0, *mcb1, *mcc0;
int *IPRSNT;
LDBLE *M, *LGAMMA;
LDBLE BK[23], DK[23];

/* routines */
int calc_pitz_param(struct pitz_param *pz_ptr, LDBLE TK, LDBLE TR);
int check_gammas_pz(void);
int ISPEC(char *name);
/*int DH_AB (LDBLE TK, LDBLE *A, LDBLE *B);*/
LDBLE G(LDBLE Y);
LDBLE GP(LDBLE Y);
#ifdef SKIP
LDBLE ETHETAP(LDBLE ZJ, LDBLE ZK, LDBLE I);
LDBLE ETHETA(LDBLE ZJ, LDBLE ZK, LDBLE I);
#endif
int ETHETAS(LDBLE ZJ, LDBLE ZK, LDBLE I, LDBLE * etheta,
				   LDBLE * ethetap);
int BDK(LDBLE X);
int pitzer_initial_guesses(void);
int pitzer_revise_guesses(void);
int pitzer_remove_unstable_phases;
int PTEMP(LDBLE TK);
LDBLE JAY(LDBLE X);
LDBLE JPRIME(LDBLE Y);
int jacobian_pz(void);

// pitzer_structures.c -------------------------------

struct pitz_param *pitz_param_alloc(void);
int pitz_param_init(struct pitz_param *pitz_param_ptr);
struct pitz_param *pitz_param_duplicate(struct pitz_param *old_ptr);
int pitz_param_copy(struct pitz_param *old_ptr,
						   struct pitz_param *new_ptr);

// pitzer_structures.c -------------------------------
int add_potential_factor(void);
int add_cd_music_factors(int n);
int add_surface_charge_balance(void);
int add_cd_music_charge_balances(int i);
int build_gas_phase(void);
int build_jacobian_sums(int k);
int build_mb_sums(void);
int build_min_exch(void);
int build_model(void);
int build_pure_phases(void);
int build_s_s_assemblage(void);
int build_solution_phase_boundaries(void);
int build_species_list(int n);
int build_min_surface(void);
int change_hydrogen_in_elt_list(LDBLE charge);
int clear(void);
int convert_units(struct solution *solution_ptr);
struct unknown *find_surface_charge_unknown(char *str_ptr, int plane);
struct master **get_list_master_ptrs(char *ptr,
											struct master *master_ptr);
int inout(void);
int is_special(struct species *spec);
int mb_for_species_aq(int n);
int mb_for_species_ex(int n);
int mb_for_species_surf(int n);
int quick_setup(void);
int resetup_master(void);
int save_model(void);
int setup_exchange(void);
int setup_gas_phase(void);
int setup_master_rxn(struct master **master_ptr_list,
							struct reaction **pe_rxn);
int setup_pure_phases(void);
int setup_related_surface(void);
int setup_s_s_assemblage(void);
int setup_solution(void);
int setup_surface(void);
int setup_unknowns(void);
int store_dn(int k, LDBLE * source, int row, LDBLE coef_in,
					LDBLE * gamma_source);
int store_jacob(LDBLE * source, LDBLE * target, LDBLE coef);
int store_jacob0(int row, int column, LDBLE coef);
int store_mb(LDBLE * source, LDBLE * target, LDBLE coef);
int store_mb_unknowns(struct unknown *unknown_ptr, LDBLE * LDBLE_ptr,
							 LDBLE coef, LDBLE * gamma_ptr);
int store_sum_deltas(LDBLE * source, LDBLE * target, LDBLE coef);
int tidy_redox(void);
struct master **unknown_alloc_master(void);
int write_mb_eqn_x(void);
int write_mb_for_species_list(int n);
int write_mass_action_eqn_x(int stop);

// print.c -------------------------------
int print_alkalinity(void);
int print_diffuse_layer(struct surface_charge *surface_charge_ptr);
int print_eh(void);
int print_irrev(void);
int print_kinetics(void);
int print_mix(void);
int print_pp_assemblage(void);
int print_s_s_assemblage(void);
int print_saturation_indices(void);
int print_surface_cd_music(void);
int print_totals(void);
int print_using(void);
/*int print_user_print(void);*/
int punch_gas_phase(void);
int punch_identifiers(void);
int punch_kinetics(void);
int punch_molalities(void);
int punch_activities(void);
int punch_pp_assemblage(void);
int punch_s_s_assemblage(void);
int punch_saturation_indices(void);
int punch_totals(void);
int punch_user_punch(void);


// read.c -------------------------------


int add_psi_master_species(char *token);
int read_advection(void);
int read_analytical_expression_only(char *ptr, LDBLE * log_k);
/* VP: Density Start */
int read_millero_abcdef (char *ptr, LDBLE * abcdef);
/* VP: Density End */
int read_copy(void);
int read_debug(void);
int read_delta_h_only(char *ptr, LDBLE * delta_h,
							 DELTA_H_UNIT * units);
int read_llnl_aqueous_model_parameters(void);
int read_exchange(void);
int read_exchange_master_species(void);
int read_exchange_species(void);
int read_gas_phase(void);
int read_incremental_reactions(void);
int read_inverse(void);
int read_inv_balances(struct inverse *inverse_ptr, char *next_char);
int read_inv_isotopes(struct inverse *inverse_ptr, char *ptr);
int read_inv_phases(struct inverse *inverse_ptr, char *next_char);
int read_kinetics(void);
int read_line_doubles(char *next_char, LDBLE ** d, int *count_d,
							 int *count_alloc);
int read_lines_doubles(char *next_char, LDBLE ** d, int *count_d,
							  int *count_alloc, const char **opt_list,
							  int count_opt_list, int *opt);
LDBLE *read_list_doubles(char **ptr, int *count_doubles);
int *read_list_ints(char **ptr, int *count_ints, int positive);
int *read_list_t_f(char **ptr, int *count_ints);
int read_master_species(void);
int read_mix(void);
int read_named_logk(void);
int read_phases(void);
int read_print(void);
int read_pure_phases(void);
int read_rates(void);
int read_reaction(void);
int read_reaction_reactants(struct irrev *irrev_ptr);
int read_reaction_steps(struct irrev *irrev_ptr);
int read_solid_solutions(void);
int read_temperature(void);
int read_reaction_temps(struct temperature *temperature_ptr);
int read_save(void);
int read_selected_output(void);
int read_solution(void);
int read_species(void);
int read_surf(void);
int read_surface_master_species(void);
int read_surface_species(void);
int read_use(void);
int read_title(void);
int read_user_print(void);
int read_user_punch(void);
#ifdef PHREEQ98
int read_user_graph(void);
 int connect_simulations, graph_initial_solutions;
/*extern*/ int shifts_as_points;
 int chart_type;
 int ShowChart;
 int RowOffset, ColumnOffset;
#endif

//extern int reading_database(void);
//extern int check_line(const char *string, int allow_empty, int allow_eof,
//					  int allow_keyword, int print);

#ifdef PHREEQ98
int copy_title(char *token_ptr, char **ptr, int *length);
int OpenCSVFile(char file_name[MAX_LENGTH]);
void GridHeadings(char *s, int i);
void SetAxisTitles(char *s, int i);
void SetAxisScale(char *a, int c, char *v, int l);
void SetChartTitle(char *s);
#endif

LDBLE dummy;
int next_keyword_or_option(const char **opt_list, int count_opt_list);

// readtr.c -------------------------------

int read_line_LDBLEs(char *next_char, LDBLE ** d, int *count_d,
							int *count_alloc);

// sit.c -------------------------------

LDBLE sit_A0;
//extern struct species **spec, **cations, **anions, **neutrals;
int sit_count_cations, sit_count_anions, sit_count_neutrals;
int sit_MAXCATIONS, sit_FIRSTANION, sit_MAXNEUTRAL;
//extern struct pitz_param *mcb0, *mcb1, *mcc0;
int *sit_IPRSNT;
LDBLE *sit_M, *sit_LGAMMA;

/* routines */
int calc_sit_param(struct pitz_param *pz_ptr, LDBLE TK, LDBLE TR);
int check_gammas_sit(void);
int sit_ISPEC(char *name);
/*int DH_AB (LDBLE TK, LDBLE *A, LDBLE *B);*/
int sit_initial_guesses(void);
int sit_revise_guesses(void);
int sit_remove_unstable_phases;
int PTEMP_SIT(LDBLE tk);
int jacobian_sit(void);

// spread.c -------------------------------

int copy_token_tab(char *token_ptr, char **ptr, int *length);
int get_option_string(const char **opt_list, int count_opt_list,
							 char **next_char);
int spread_row_free(struct spread_row *spread_row_ptr);
int spread_row_to_solution(struct spread_row *heading,
								  struct spread_row *units,
								  struct spread_row *data,
								  struct defaults defaults);
struct spread_row *string_to_spread_row(char *string);
#ifdef PHREEQCI_GUI
void add_row(struct spread_row *spread_row_ptr);
void copy_defaults(struct defaults *dest_ptr,
						  struct defaults *src_ptr);
void free_spread(void);
struct spread_row *copy_row(struct spread_row *spread_row_ptr);
#endif

// step.c -------------------------------
//static char const svnid[] = "$Id: step.c 3453 2009-04-14 20:42:49Z dlpark $";

int check_pp_assemblage(struct pp_assemblage *pp_assemblage_ptr);
int gas_phase_check(struct gas_phase *gas_phase_ptr);
int pp_assemblage_check(struct pp_assemblage *pp_assemblage_ptr);
int reaction_calc(struct irrev *irrev_ptr);
int solution_check(void);
int s_s_assemblage_check(struct s_s_assemblage *s_s_assemblage_ptr);
// structures.c -------------------------------

//static char const svnid[] =
//	"$Id: structures.c 3716 2009-10-22 19:24:31Z dlpark $";

static int exchange_compare_int(const void *ptr1, const void *ptr2);
static int gas_phase_compare_int(const void *ptr1, const void *ptr2);

static int inverse_compare(const void *ptr1, const void *ptr2);
int inverse_free(struct inverse *inverse_ptr);

static int irrev_compare(const void *ptr1, const void *ptr2);
static int irrev_compare_int(const void *ptr1, const void *ptr2);

static int kinetics_compare_int(const void *ptr1, const void *ptr2);

int logk_init(struct logk *logk_ptr);

static int master_compare_string(const void *ptr1, const void *ptr2);
int master_free(struct master *master_ptr);

#if defined(PHREEQCI_GUI)
static int mix_compare(const void *ptr1, const void *ptr2);
#else
static int mix_compare(const void *ptr1, const void *ptr2);
#endif
static int mix_compare_int(const void *ptr1, const void *ptr2);

struct phase *phase_alloc(void);
static int phase_compare_string(const void *ptr1, const void *ptr2);
int phase_free(struct phase *phase_ptr);
int phase_init(struct phase *phase_ptr);

static int pp_assemblage_compare_int(const void *ptr1, const void *ptr2);

static int rate_compare(const void *ptr1, const void *ptr2);
static int rate_compare_string(const void *ptr1, const void *ptr2);

struct species *s_alloc(void);
int s_free(struct species *s_ptr);
int s_init(struct species *s_ptr);

static int s_s_assemblage_compare_int(const void *ptr1, const void *ptr2);

static int solution_compare(const void *ptr1, const void *ptr2);
static int solution_compare_int(const void *ptr1, const void *ptr2);

static int species_list_compare(const void *ptr1, const void *ptr2);

static int surface_compare_int(const void *ptr1, const void *ptr2);

static int temperature_compare(const void *ptr1, const void *ptr2);
static int temperature_compare_int(const void *ptr1, const void *ptr2);

static int rxn_token_temp_compare(const void *ptr1, const void *ptr2);
int trxn_multiply(LDBLE coef);

#ifdef PHREEQCI_GUI
//extern void free_spread(void);
#endif
#if defined(USE_MPI) && defined(HDF5_CREATE) && defined(MERGE_FILES)
extern void MergeFinalize(void);
#endif
//extern LDBLE *scratch, *x_arg, *res_arg;
//extern LDBLE *normal, *ineq_array, *zero, *res, *delta1, *cu;
//extern int *iu, *is, *back_eq;
//extern int x_arg_max, res_arg_max, scratch_max;

// tidy.c -------------------------------

//static char const svnid[] = "$Id: tidy.c 3684 2009-09-28 15:40:35Z dlpark $";

int check_species_input(void);
LDBLE coef_in_master(struct master *master_ptr);
int phase_rxn_to_trxn(struct phase *phase_ptr,
							 struct reaction *rxn_ptr);
int reset_last_model(void);
int rewrite_eqn_to_primary(void);
int rewrite_eqn_to_secondary(void);
int species_rxn_to_trxn(struct species *s_ptr);
int tidy_logk(void);
int tidy_exchange(void);
int tidy_min_exchange(void);
int tidy_kin_exchange(void);
int tidy_gas_phase(void);
int tidy_inverse(void);
int tidy_isotopes(void);
int tidy_isotope_ratios(void);
int tidy_isotope_alphas(void);
int tidy_kin_surface(void);
int tidy_master_isotope(void);
int tidy_min_surface(void);
int tidy_phases(void);
int tidy_pp_assemblage(void);
int tidy_solutions(void);
int tidy_s_s_assemblage(void);
int tidy_species(void);
int tidy_surface(void);

LDBLE a0, a1, kc, kb;
#if !defined(PHREEQC_CLASS)
int scan(LDBLE f(LDBLE x), LDBLE * xx0, LDBLE * xx1);
static LDBLE f_spinodal(LDBLE x);
#else
int scan(LDBLE f(LDBLE x, void *), LDBLE * xx0, LDBLE * xx1);
static LDBLE f_spinodal(LDBLE x, void *);
#endif
int solve_misc(LDBLE * xxc1, LDBLE * xxc2, LDBLE tol);
int s_s_calc_a0_a1(struct s_s *s_s_ptr);

// tally.c -------------------------------

//static char const svnid[] = "$Id: tally.c 3164 2008-10-29 22:15:35Z dlpark $";
/*
 *  storage
 */
struct tally_buffer
{
	char *name;
	struct master *master;
	LDBLE moles;
	LDBLE gfw;
};
struct tally_buffer *t_buffer;
int tally_count_component;
struct tally
{
	char *name;
	enum entity_type type;
	char *add_formula;
	LDBLE moles;
	struct elt_list *formula;
	/*
	 * first total is initial
	 * second total is final
	 * third total is difference (final - initial)
	 */
	struct tally_buffer *total[3];
};
struct tally *tally_table;
int count_tally_table_columns;
int count_tally_table_rows;



// transport.c -------------------------------

//static char const svnid[] =
//	"$Id: transport.c 3694 2009-10-02 17:54:58Z dlpark $";

struct spec
{
	char *name;					/* name of species */
	char *aq_name;				/* name of aqueous species in EX species */
	int type;					/* type: AQ or EX */
	LDBLE a;					/* activity */
	LDBLE lm;					/* log(concentration) */
	LDBLE lg;					/* log(gamma) */
	LDBLE c;					/* concentration for AQ, equivalent fraction for EX */
	LDBLE z;					/* charge number */
	LDBLE Dwt;					/* temperature corrected free water diffusion coefficient, m2/s */
	LDBLE erm_ddl;				/* enrichment factor in ddl */
};
struct sol_D
{
	int count_spec;				/* number of aqueous + exchange species */
	int count_exch_spec;		/* number of exchange species */
	LDBLE exch_total;			/* total moles of X- */
	struct spec *spec;
} *sol_D;
struct sol_D *sol_D_dbg;
struct J_ij
{
	char *name;
	LDBLE tot1, tot2;
} *J_ij, *J_ij_il;
int J_ij_count_spec;

struct M_S
{
	char *name;
	LDBLE tot1, tot2;
} *m_s;
int count_m_s;
LDBLE tot1_h, tot1_o, tot2_h, tot2_o;

int multi_D(LDBLE DDt, int mobile_cell, int stagnant);
int find_J(int icell, int jcell, LDBLE mixf, LDBLE DDt, int stagnant);
int fill_spec(int cell_no);
int fill_m_s(struct J_ij *J_ij, int J_ij_count_spec);
static int sort_species_name(const void *ptr1, const void *ptr2);

LDBLE diffc_max, diffc_tr, J_ij_sum;
int transp_surf;
int disp_surf(LDBLE stagkin_time);
int diff_stag_surf(int mobile_cell);
int check_surfaces(struct surface *surface_ptr1,
						  struct surface *surface_ptr2);
int mobile_surface_copy(struct surface *surface_old_ptr,
							   struct surface *surf_ptr1, int n_user_new,
							   int move_old);
int init_mix(void);
int init_heat_mix(int nmix);
int heat_mix(int heat_nmix);
int mix_stag(int i, LDBLE stagkin_time, int punch,
					LDBLE step_fraction_kin);
LDBLE *heat_mix_array;
LDBLE *temp1, *temp2;
int nmix, heat_nmix;
LDBLE heat_mix_f_imm, heat_mix_f_m;
int warn_MCD_X, warn_fixed_Surf;

// utilities.c -------------------------------

//static char const svnid[] =
//	"$Id: utilities.c 3164 2008-10-29 22:15:35Z dlpark $";

#ifdef PHREEQ98
extern int AutoLoadOutputFile, CreateToC;
extern int ProcessMessages, ShowProgress, ShowProgressWindow, ShowChart;
extern int outputlinenr;
extern int stop_calculations;
void AddToCEntry(char *a, int l, int i);
void ApplicationProcessMessages(void);
/* void check_line_breaks(char *s); */
char err_str98[80];
int copy_title(char *token_ptr, char **ptr, int *length);
extern int clean_up_null(void);
#endif

int isamong(char c, const char *s_l);

Address Hash_multi(HashTable * Table, char *Key);
void ExpandTable_multi(HashTable * Table);





#ifdef PHREEQ98
int punch_user_graph(void);
int colnr, rownr;
int graph_initial_solutions;
int prev_advection_step, prev_transport_step;	/*, prev_reaction_step */
/* int shifts_as_points; */
int chart_type;
int AddSeries;
int FirstCallToUSER_GRAPH;
#endif

#if defined(SWIG_SHARED_OBJ)
int EndRow(void);
void AddSelectedOutput(const char *name, const char *format,
							  va_list argptr);
#endif


// Remove static
#define STATIC
#define EXTERNAL

private:
  #include "phrqproto.h"


public:
	void set_phast(int);
	int main_method(int argc, char *argv[]);

};
#endif /* _INC_PHREEQC_H */

