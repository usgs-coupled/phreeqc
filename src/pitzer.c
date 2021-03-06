#define EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "output.h"
#include "phrqproto.h"
#define PITZER
#define PITZER_EXTERNAL
#include "pitzer.h"

static char const svnid[] = "$Id: pitzer.c 248 2005-04-14 17:10:53Z dlpark $";
/* variables */
static double A0;
struct species **spec, **cations, **anions, **neutrals;
static int count_cations, count_anions, count_neutrals;
static int MAXCATIONS, FIRSTANION, MAXNEUTRAL;
struct pitz_param *mcb0, *mcb1, *mcc0;
static int *IPRSNT;
static double *M, *LGAMMA;
static double BK[23], DK[23];

/* routines */
static int calc_pitz_param (struct pitz_param *pz_ptr, double TK, double TR);
static int check_gammas_pz(void);
static int ISPEC(char * name);
/*static int DH_AB (double TK, double *A, double *B);*/
static double G (double Y);
static double GP (double Y);
#ifdef SKIP
static double ETHETAP (double ZJ, double ZK, double I);
static double ETHETA (double ZJ, double ZK, double I);
#endif
static int ETHETAS (double ZJ, double ZK, double I, double *etheta, double *ethetap);
static int BDK (double X);
static int initial_guesses(void);
static int revise_guesses(void);
static int remove_unstable_phases;


/* ---------------------------------------------------------------------- */
int pitzer_init (void)
/* ---------------------------------------------------------------------- */
{
	int i;
/*
 *      Initialization for pitzer
 */
	pitzer_model = FALSE;
	max_pitz_param = 100;
	count_pitz_param = 0;
	space ((void **) ((void *) &pitz_params), INIT, &max_pitz_param, sizeof(struct pitz_param *));

	max_theta_param = 100;
	count_theta_param = 0;
	space ((void **) ((void *) &theta_params), INIT, &max_theta_param, sizeof(struct theta_param *));

	ICON = TRUE;
	OTEMP=0.0;
	for (i = 0; i < 23; i++) {
		BK[i] = 0.0;
		DK[i] = 0.0;
	}
	return OK;
}
/* ---------------------------------------------------------------------- */
int pitzer_tidy (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Make lists of species for cations, anions, neutral
 */
	char * string1, *string2;
	int i, j, order;
	double z0, z1;
	struct pitz_param *pzp_ptr;
	struct theta_param *theta_param_ptr;
	/*
	 *  allocate pointers to species structures
	 */
	if (spec != NULL) spec = (struct species **) free_check_null(spec);
	spec = (struct species **) PHRQ_malloc((size_t) (3*count_s*sizeof(struct species *)));
	if (spec == NULL) malloc_error();
	for (i = 0; i < 3*count_s; i++) spec[i] = NULL;
	cations = spec;
	neutrals = &(spec[count_s]);
	anions = &(spec[2*count_s]);
	MAXCATIONS = count_s;
	FIRSTANION = 2*count_s;
	MAXNEUTRAL = count_s;
	count_cations = 0;
	count_anions = 0;
	count_neutrals = 0;
	if (itmax < 200) itmax = 200;
	/*
	 *  allocate other arrays for Pitzer
	 */
	if (IPRSNT != NULL) IPRSNT = (int *) free_check_null(IPRSNT);
	IPRSNT = (int *) PHRQ_malloc((size_t) (3*count_s*sizeof(int)));
	if (IPRSNT == NULL) malloc_error();
	if (M != NULL) M = (double *) free_check_null(M);
	M = (double *) PHRQ_malloc((size_t) (3*count_s*sizeof(double)));
	if (M == NULL) malloc_error();
	if (LGAMMA != NULL) LGAMMA = (double *) free_check_null(LGAMMA);
	LGAMMA = (double *) PHRQ_malloc((size_t) (3*count_s*sizeof(double)));
	if (LGAMMA == NULL) malloc_error();
	

	for (i = 0; i < count_s; i++) {
		if (s[i] == s_eminus) continue;
		if (s[i] == s_h2o) continue;
		if (s[i]->z < -.001) {
			anions[count_anions++] = s[i];
		} else if (s[i]->z > .001) {
			cations[count_cations++] = s[i];
		} else {
			neutrals[count_neutrals++] = s[i];
		}
	}
	/*
	 *  Add etheta to parameter list in case theta not defined for 
         *  cation-cation or anion-anion pair
	 *  Remove old TYPE_ETHETA definitions
	 */
	j = 0;
	for (i = 0; i < count_pitz_param; i++) {
		if (pitz_params[i]->type == TYPE_ETHETA) {
			pitz_params[i] = (struct pitz_param *) free_check_null(pitz_params[i]);
		} else {
			pitz_params[j++] = pitz_params[i];
		}
	}
	count_pitz_param = j;
	for (i = 0; i < count_cations - 1; i++) {
		for (j = i+1; j < count_cations; j++) {
			sprintf(line,"%s %s 1", spec[i]->name, spec[j]->name);
			pzp_ptr = pitz_param_read(line, 2);
			pzp_ptr->type = TYPE_ETHETA;
			if (count_pitz_param >= max_pitz_param) {
				space ((void **) ((void *) &pitz_params), count_pitz_param, &max_pitz_param, sizeof(struct pitz_param *));
			}
			pitz_params[count_pitz_param++] = pzp_ptr;
			
		}
	}
	for (i = 2*count_s; i < 2*count_s + count_anions - 1; i++) {
		for (j = i+1; j < 2*count_s + count_anions; j++) {
			sprintf(line,"%s %s 1", spec[i]->name, spec[j]->name);
			pzp_ptr = pitz_param_read(line, 2);
			pzp_ptr->type = TYPE_ETHETA;
			if (count_pitz_param >= max_pitz_param) {
				space ((void **) ((void *) &pitz_params), count_pitz_param, &max_pitz_param, sizeof(struct pitz_param *));
			}
			pitz_params[count_pitz_param] = pzp_ptr;
			count_pitz_param++;
		}
	}
	/*
	 *  put species numbers in pitz_params
	 */
	for (i = 0; i < count_pitz_param; i++) {
		for (j = 0; j < 3; j++) {
			if (pitz_params[i]->species[j] == NULL) continue;
			pitz_params[i]->ispec[j] = ISPEC(pitz_params[i]->species[j]);
			if ((j < 2 && pitz_params[i]->ispec[j] == -1) ||
			     (j == 3 && (pitz_params[i]->type == TYPE_PSI || pitz_params[i]->type == TYPE_ZETA) && pitz_params[i]->ispec[j] == -1)) {
				input_error++;
				sprintf(error_string, "Species for Pitzer parameter not defined in SOLUTION_SPECIES, %s", pitz_params[i]->species[j]);
				error_msg(error_string, CONTINUE);
				return(ERROR);
			}
		}
	}
	/*
	 * McGinnis data
	 */
	string1 = string_hsave("K+");
	string2 = string_hsave("Cl-");
	IC = ISPEC(string2);
	for (i = 0; i < count_pitz_param; i++) {
		if (pitz_params[i]->species[0] == string1 &&
		    pitz_params[i]->species[1] == string2) {
			switch (pitz_params[i]->type) {
			case TYPE_B0:
				mcb0 = pitz_params[i];
				break;
			case TYPE_B1:
				mcb1 = pitz_params[i];
				break;
			case TYPE_C0:
				mcc0 = pitz_params[i];
				break;
			case TYPE_B2:
			case TYPE_THETA:
			case TYPE_LAMDA:
			case TYPE_ZETA:
			case TYPE_PSI:
			case TYPE_ETHETA:
			case TYPE_Other:
				break;
			}
		}
	}
	/*
	 * Set alpha values
	 */
	for (i = 0; i < count_pitz_param; i++) {
		z0 = fabs(spec[pitz_params[i]->ispec[0]]->z);
		z1 = fabs(spec[pitz_params[i]->ispec[1]]->z);
		if (equal(z0, 1.0, 1e-8) || equal(z1, 1.0, 1e-8)) {
			order = 1;
		} else if (equal(z0,2.0, 1e-8) && equal(z1, 2.0, 1e-8)) {
			order = 2;
		} else {
			order = 3;
		}
		if (pitz_params[i]->type == TYPE_B1) {
			switch (order) {
			case 1:
			case 3:
				pitz_params[i]->alpha = 2.0;
				break;
			case 2:
				pitz_params[i]->alpha = 1.4;
				break;
			}
		} else if (pitz_params[i]->type == TYPE_B2) {
			switch (order) {
			case 1:
				pitz_params[i]->alpha = 12.0;
				break;
			case 2:
				pitz_params[i]->alpha = 12.0;
				break;
			case 3:
				pitz_params[i]->alpha = 50.0;
				break;
			}
		}
	}
	/*
	 *   Add thetas pointer to etheta pitzer parameters
	 */

	if (count_theta_param > 0 ) {
		for (i = 0; i < count_theta_param; i++) {
			theta_params[i] = (struct theta_param *) free_check_null(theta_params[i]);
		}
	}
	count_theta_param = 0;
	for (i = 0; i < count_pitz_param; i++) {
		if (pitz_params[i]->type == TYPE_ETHETA) {
			z0 = spec[pitz_params[i]->ispec[0]]->z;
			z1 = spec[pitz_params[i]->ispec[1]]->z;
			theta_param_ptr = theta_param_search(z0, z1);
			if (theta_param_ptr == NULL) {
				if (count_theta_param >= max_theta_param) {
					space ((void **) ((void *) &theta_params), count_theta_param, &max_theta_param, sizeof(struct theta_param *));
				}
				theta_params[count_theta_param] = theta_param_alloc();
				theta_param_init(theta_params[count_theta_param]);
				theta_params[count_theta_param]->zj = z0;
				theta_params[count_theta_param]->zk = z1;
				theta_param_ptr = theta_params[count_theta_param];
				count_theta_param++;
			}
			pitz_params[i]->thetas = theta_param_ptr;
		}
	}
	return OK;
}
/* ---------------------------------------------------------------------- */
int ISPEC(char * name)
/* ---------------------------------------------------------------------- */
/*
 *      Find species number in spec for character string species name
 */
{
	int i;
	for (i = 0; i < 3*count_s; i++) {
		if (spec[i] == NULL) continue;
		if (name == spec[i]->name) {
			return(i);
		}
	}
	return (-1);
}
/* ---------------------------------------------------------------------- */
int read_pitzer (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads advection information
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
/*
 *   Read advection parameters: 
 *        number of cells;
 *        number of shifts;
 */
	int n, j;
	struct pitz_param *pzp_ptr;
	pitz_param_type pzp_type;

	int return_value, opt, opt_save;
	char *next_char;
	const char *opt_list[] = {
		"b0",                    /* 0 */
		"b1",                    /* 1 */
		"b2",                    /* 2 */
		"c0",                    /* 3 */
		"theta",                 /* 4 */
		"lamda",                 /* 5 */
		"zeta",                  /* 6 */
		"psi",                   /* 7 */
		"macinnes",              /* 8 */
		"macinnis",              /* 9 */
		"mac"                    /* 10 */
	};
	int count_opt_list = 11;
/*
 *   Read lines
 */
	opt_save = OPTION_ERROR;
	return_value = UNKNOWN;
	n = -1;
	pzp_type = TYPE_Other;
	for (;;) {
		opt = get_option(opt_list, count_opt_list, &next_char);
		if (opt == OPTION_DEFAULT) {
			opt = opt_save;
		}
		switch (opt) {
		case OPTION_EOF:               /* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:           /* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_DEFAULT:
			    pzp_ptr = pitz_param_read(line, n);
			    if (pzp_ptr != NULL) {
				    pzp_ptr->type = pzp_type;
				    j = pitz_param_search(pzp_ptr);
				    if (j < 0) {
					    if (count_pitz_param >= max_pitz_param) {
						    space ((void **) ((void *) &pitz_params), count_pitz_param, &max_pitz_param, sizeof(struct pitz_param *));
					    }
					    
					    pitz_params[count_pitz_param] = pzp_ptr;
					    count_pitz_param++;
				    } else {
					    pitz_params[j] = (struct pitz_param *) free_check_null(pitz_params[j]);
					    pitz_params[j] = pzp_ptr;
				    }
			    }
			    break;
		case OPTION_ERROR:
			input_error++;
			error_msg("Unknown input in PITZER keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case 0:                        /* b0 */
			pzp_type = TYPE_B0;
			n = 2;
			opt_save = OPTION_DEFAULT;
			break;
		case 1:                        /* b1 */
			pzp_type = TYPE_B1;
			n = 2;
			opt_save = OPTION_DEFAULT;
			break;
		case 2:                        /* b2 */
			pzp_type = TYPE_B2;
			n = 2;
			opt_save = OPTION_DEFAULT;
			break;
		case 3:                        /* c0 */
			pzp_type = TYPE_C0;
			n = 2;
			opt_save = OPTION_DEFAULT;
			break;
		case 4:                        /* theta */
			pzp_type = TYPE_THETA;
			n = 2;
			opt_save = OPTION_DEFAULT;
			break;
		case 5:                        /* lambda */
			pzp_type = TYPE_LAMDA;
			n = 2;
			opt_save = OPTION_DEFAULT;
			break;
		case 6:                        /* zeta */
			pzp_type = TYPE_ZETA;
			n = 3;
			opt_save = OPTION_DEFAULT;
			break;
		case 7:                        /* psi */
			pzp_type = TYPE_PSI;
			n = 3;
			opt_save = OPTION_DEFAULT;
			break;
		case 8:                        /* macinnes */
		case 9:                        /* macinnis */
		case 10:                       /* mac */
			opt_save = OPTION_ERROR;
			ICON = get_true_false(next_char, TRUE);
			break;
		}
		if (return_value == EOF || return_value == KEYWORD) break;
	}
	if (count_pitz_param > 0) pitzer_model = TRUE;
	return(return_value);
}
/* ---------------------------------------------------------------------- */
int PTEMP (double TK)
/* ---------------------------------------------------------------------- */
{
/*
C
C     SUBROUTINE TO CALUCLATE TEMPERATURE DEPENDENCE OF PITZER PARAMETER
C
*/
	double DC0;
	int i;
	double TR=298.15;

	if (fabs(TK-OTEMP) < 0.01e0) return OK;
	OTEMP=TK;
/*
C     Set DW0
*/
	DW(TK);
	for (i = 0; i < count_pitz_param; i++) {
		calc_pitz_param(pitz_params[i], TK, TR);
	}
	DC0=DC(TK);
	if (fabs(TK-TR) < 0.01e0) {
		A0=0.392e0;
	} else {
		DC0=DC(TK);
		A0=1.400684e6*sqrt(DW0/(pow((DC0*TK),3.0e0)));
		/*A0=1.400684D6*(DW0/(DC0*TK)**3.0D0)**0.5D0*/
	}
	return OK;
}

/* ---------------------------------------------------------------------- */
int calc_pitz_param (struct pitz_param *pz_ptr, double TK, double TR)
/* ---------------------------------------------------------------------- */
{
	double param;
	/*
	*/

	if (fabs(TK-TR) < 0.01) {
		param = pz_ptr->a[0];
	} else {
		param = (pz_ptr->a[0] + 
			 pz_ptr->a[1]*(1.e0/TK-1.e0/TR) +
			 pz_ptr->a[2]*log(TK/TR) +
			 pz_ptr->a[3]*(TK - TR) +
			 pz_ptr->a[4]*(TK*TK - TR*TR)); 
	}
	pz_ptr->p = param;
	switch 	(pz_ptr->type) {
	    case TYPE_B0:
	        pz_ptr->U.b0 = param;
	        break;
	    case TYPE_B1:
		pz_ptr->U.b1 = param;
		break;
	    case TYPE_B2:
		pz_ptr->U.b2 = param;
		break;
	    case TYPE_C0:
		pz_ptr->U.c0 = param;
		break;
	    case TYPE_THETA:
		pz_ptr->U.theta = param;
		break;
	    case TYPE_LAMDA:
		pz_ptr->U.lamda = param;
		break;
	    case TYPE_ZETA:
		pz_ptr->U.zeta = param;
		break;
	    case TYPE_ETHETA:
		    break;
	    case TYPE_PSI:
		pz_ptr->U.psi = param;
		break;
	    case TYPE_Other:
	        error_msg("Should not be TYPE_Other in function calc_pitz_param", STOP);
	        break;
	}
	return OK;
}
/* ---------------------------------------------------------------------- */
int pitzer (void)
/* ---------------------------------------------------------------------- */
{
	int i, i0, i1, i2;
	double param, alpha, z0, z1, z2;
	double etheta, ethetap;
	double dummy;
	/*
	double CONV, XI, XX, OSUM, BIGZ, DI, F, XXX, GAMCLM, 
		CSUM, PHIMAC, OSMOT, BMXP, ETHEAP, CMX, BMX, PHI,
		BMXPHI, PHIPHI, AW, A, B;
	*/
	double CONV, XI, XX, OSUM, BIGZ, DI, F, XXX, GAMCLM, 
		CSUM, PHIMAC, OSMOT, B;
	double I, TK;
	int LNEUT;
	/*
	  C
	  C     INITIALIZE
	  C
	*/
	CONV = 1.0/log(10.0);
	XI=0.0e0;
	XX=0.0e0;
	OSUM=0.0e0;
	LNEUT=FALSE;
	/*n
	I = *I_X;
	TK = *TK_X;
	*/
	I = mu_x;
	TK = tk_x;
	/*	DH_AB(TK, &A, &B); */
	/*
	  C
	  C     TRANSFER DATA FROM TO M
	  C
	*/
	for (i = 0; i < 3*count_s; i++) {
		IPRSNT[i] = FALSE;
		M[i] = 0.0;
		if (spec[i] != NULL && spec[i]->in == TRUE) {
			if (spec[i]->type == EX ||
			    spec[i]->type == SURF ||
			    spec[i]->type == SURF_PSI) continue;
			M[i] = under(spec[i]->lm);
			if (M[i] > MIN_TOTAL) IPRSNT[i] = TRUE;
		}
	}
	if (ICON == TRUE) {
		IPRSNT[IC] = TRUE;
	}
#ifdef SKIP
	for (i = count_s; i < count_s + count_neutrals; i++) {
		if (M[i] > MIN_TOTAL) LNEUT = TRUE;
	}
#endif
	/*
	ICON = 0;
	M[1] = 1.40070736;
	M[4] = 2.52131086E-05;
	M[140] = 4.59985435E-09;
	*/

/*
C
C     COMPUTE PITZER COEFFICIENTS' TEMPERATURE DEPENDENCE
C
*/
	PTEMP(TK);
	for (i = 0; i < 2*count_s + count_anions; i++) {
		LGAMMA[i] = 0.0;
		if (IPRSNT[i] == TRUE) {
			XX=XX+M[i]*fabs(spec[i]->z);
			XI=XI+M[i]*spec[i]->z*spec[i]->z;
			OSUM=OSUM+M[i];
		}
	}
	I=XI/2.0e0;
/*
C
C     EQUATION (8)
C
*/
	BIGZ=XX;
	DI=sqrt(I);
/*
C
C     CALCULATE F & GAMCLM
C
*/
	B = 1.2;
	F=-A0*(DI/(1.0e0+B*DI)+2.0e0*log(1.0e0+B*DI)/B);
	XXX=2.0e0*DI;
	XXX=(1.0e0-(1.0e0+XXX-XXX*XXX*0.5e0)*exp(-XXX))/(XXX*XXX);
	/*GAMCLM=F+I*2.0e0*(BCX(1,IK,IC)+BCX(2,IK,IC)*XXX)+1.5e0*BCX(4,IK,IC)*I*I;*/
    	/*GAMCLM=F+I*2.0e0*(mcb0->U.b0 + mcb1->U.b1*XXX) + 1.5e0*mcc0->U.c0*I*I;*/
	GAMCLM=F+I*2.0e0*(mcb0->p + mcb1->p*XXX) + 1.5e0*mcc0->p*I*I;
	CSUM=0.0e0;
	OSMOT=-(A0)*pow(I,1.5e0)/(1.0e0+B*DI);
/*
 *  Calculate ethetas
 */
	for (i = 0; i < count_theta_param; i++) {
		z0 = theta_params[i]->zj;
		z1 = theta_params[i]->zk;
		ETHETAS(z0, z1, I, &etheta, &ethetap);
		theta_params[i]->etheta = etheta;
		theta_params[i]->ethetap = ethetap;
	}
/*
 *  Sums for F, LGAMMA, and OSMOT
 */
	dummy = LGAMMA[1];
	for (i = 0; i < count_pitz_param; i++) {
		i0 = pitz_params[i]->ispec[0];
		i1 = pitz_params[i]->ispec[1];
		if (IPRSNT[i0] == FALSE || IPRSNT[i1] == FALSE) continue;
		z0 = spec[i0]->z;
		z1 = spec[i1]->z;
		param = pitz_params[i]->p;
		alpha = pitz_params[i]->alpha;
		switch (pitz_params[i]->type) {
		case TYPE_B0:
			LGAMMA[i0] += M[i1]*2.0*param;
			LGAMMA[i1] += M[i0]*2.0*param;
			OSMOT += M[i0]*M[i1]*param;
			break;
		case TYPE_B1:
			F += M[i0]*M[i1]*param*GP(alpha*DI)/I; 
			LGAMMA[i0] += M[i1]*2.0*param*G(alpha*DI);
			LGAMMA[i1] += M[i0]*2.0*param*G(alpha*DI);
			OSMOT += M[i0]*M[i1]*param*exp(-alpha*DI);
			break;
		case TYPE_B2:
			F += M[i0]*M[i1]*param*GP(alpha*DI)/I; 
			LGAMMA[i0] += M[i1]*2.0*param*G(alpha*DI);
			LGAMMA[i1] += M[i0]*2.0*param*G(alpha*DI);
			OSMOT += M[i0]*M[i1]*param*exp(-alpha*DI);
			break;
		case TYPE_C0:
			CSUM += M[i0]*M[i1]*pitz_params[i]->p/(2.0e0*sqrt(fabs(z0*z1)));
			LGAMMA[i0] += M[i1]*BIGZ*param/(2.0*sqrt(fabs(z0*z1))); 
			LGAMMA[i1] += M[i0]*BIGZ*param/(2.0*sqrt(fabs(z0*z1))); 
			OSMOT += M[i0]*M[i1]*BIGZ*param/(2.0*sqrt(fabs(z0*z1))); 
			break;
		case TYPE_THETA:
			LGAMMA[i0] += 2.0*M[i1]*(param /*+ ETHETA(z0, z1, I) */ ); 
			LGAMMA[i1] += 2.0*M[i0]*(param /*+ ETHETA(z0, z1, I) */ ); 
			OSMOT += M[i0]*M[i1]*param;
			break;
		case TYPE_ETHETA:
			/*
			ETHETAS(z0, z1, I, &etheta, &ethetap);
			*/
			etheta = pitz_params[i]->thetas->etheta;
			ethetap = pitz_params[i]->thetas->ethetap;
			F += M[i0]*M[i1]*ethetap;
			LGAMMA[i0] += 2.0*M[i1]*etheta; 
			LGAMMA[i1] += 2.0*M[i0]*etheta; 
			OSMOT += M[i0]*M[i1]*(etheta + I*ethetap); 
			/*
			F += M[i0]*M[i1]*ETHETAP(z0, z1, I);
			LGAMMA[i0] += 2.0*M[i1]*(ETHETA(z0, z1, I) ); 
			LGAMMA[i1] += 2.0*M[i0]*(ETHETA(z0, z1, I) ); 
			OSMOT += M[i0]*M[i1]*(ETHETA(z0, z1, I) + I*ETHETAP(z0, z1, I) ); 
			*/
			break;
		case TYPE_PSI:
			i2 = pitz_params[i]->ispec[2];
			if (IPRSNT[i2] == FALSE) continue;
			z2 = spec[i2]->z;
			LGAMMA[i0] += M[i1]*M[i2]*param;
			LGAMMA[i1] += M[i0]*M[i2]*param;
			LGAMMA[i2] += M[i0]*M[i1]*param;
			OSMOT += M[i0]*M[i1]*M[i2]*param;
			break;
		case TYPE_LAMDA:
			LGAMMA[i0] += 2.0*M[i1]*param;
			LGAMMA[i1] += 2.0*M[i0]*param;
			OSMOT += M[i0]*M[i1]*param;
			break;
		case TYPE_ZETA:
			i2 = pitz_params[i]->ispec[2];
			if (IPRSNT[i2] == FALSE) continue;
			LGAMMA[i0] += M[i1]*M[i2]*param;
			LGAMMA[i1] += M[i0]*M[i2]*param;
			LGAMMA[i2] += M[i0]*M[i1]*param;
			OSMOT += M[i0]*M[i1]*M[i2]*param;
			break;
		case TYPE_Other:
			error_msg("TYPE_Other in pitz_param list.", STOP);
			break;
		}
	}
	
	/*
	 *  Add F and CSUM terms to LGAMMA
	 */

	for (i = 0; i < count_cations; i++) {
		z0 = spec[i]->z;
		LGAMMA[i] += z0*z0*F+fabs(z0)*CSUM;
	}
	for (i = 2*count_s; i < 2*count_s + count_anions; i++) {
		z0 = spec[i]->z;
		LGAMMA[i] += z0*z0*F+fabs(z0)*CSUM;
	}
/*
C
C     CONVERT TO MACINNES CONVENTION
C
*/
      if (ICON == TRUE) {
	      PHIMAC=LGAMMA[IC]-GAMCLM;
/*
C
C     CORRECTED ERROR IN PHIMAC, NOVEMBER, 1989
C
*/
	      for (i = 0; i < 2*count_s + count_anions; i++) {
		      if (IPRSNT[i] == TRUE) {
			      LGAMMA[i]=LGAMMA[i]+spec[i]->z*PHIMAC;
		      }
	      }
      }

      COSMOT = 1.0e0 + 2.0e0*OSMOT/OSUM;
/*
C
C     CALCULATE THE ACTIVITY OF WATER
C
*/
      AW=exp(-OSUM*COSMOT/55.50837e0);
      if (AW > 1.0) AW = 1.0;
      /*s_h2o->la=log10(AW);*/
      mu_x = I;
      for (i = 0; i < 2*count_s + count_anions; i++) {
	      if (IPRSNT[i] == FALSE) continue;
	      /*spec[i]->lg=LGAMMA[i]*CONV;*/
	      spec[i]->lg_pitzer=LGAMMA[i]*CONV;
	      /*
	      output_msg(OUTPUT_MESSAGE, "%d %s:\t%e\t%e\t%e\t%e \n", i, spec[i]->name, M[i], spec[i]->la, spec[i]->lg_pitzer, spec[i]->lg);
	      */
      }
      /*
      output_msg(OUTPUT_MESSAGE, "OSUM: %e\n", OSUM);
      output_msg(OUTPUT_MESSAGE, "OSMOT: %e\n", OSMOT);
      output_msg(OUTPUT_MESSAGE, "COSMOT: %e\n", COSMOT);
      output_msg(OUTPUT_MESSAGE, "F: %e\n", F);
      output_msg(OUTPUT_MESSAGE, "AW: %e\n", AW);
      */
      /*
      *I_X = I;
      *COSMOT_X = COSMOT;
      */
      return(OK);
}
/* ---------------------------------------------------------------------- */
double JAY (double X)
/* ---------------------------------------------------------------------- */
/*
C
C     FUNCTION TO CALCULATE JAY AND JPRIME
C
C     J0 AND J1, USED IN CALCULATION OF ETHETA AND ETHEAP
C
*/
{
	double JAY;
	BDK (X);
	JAY=X/4.0e0-1.0e0+0.5e0*(BK[0]-BK[2]);
	return JAY;
}
/* ---------------------------------------------------------------------- */
double JPRIME (double Y)
/* ---------------------------------------------------------------------- */
{
	double DZ;
	BDK (Y);
	if (Y > 1.0e0) {
		DZ=-4.0e0*pow(Y,-1.1e0)/9.0e0;
	} else {
		DZ=0.8e0*pow(Y,-0.8e0);
	}			   
	return (Y*(.25e0+DZ*(DK[0]-DK[2])/2.0e0));
}


/* ---------------------------------------------------------------------- */
int BDK (double X)
/* ---------------------------------------------------------------------- */
/*
C
C     NUMERICAL APPROXIMATION TO THE INTEGRALS IN THE EXPRESSIONS FOR J0
C     AND J1.  CHEBYSHEV APPROXIMATION IS USED.  THE CONSTANTS 'AK' ARE
C     DEFINED IN BLOCK COMMON.
C
*/
/*
C
C     AK IS USED TO CALCULATE HIGHER ORDER ELECTROSTATIC TERMS IN
C     SUBROUTINE PITZER
C
*/
{
     double AKX[42] = { 
         1.925154014814667e0, -.060076477753119e0, -.029779077456514e0,
         -.007299499690937e0, 0.000388260636404e0, 0.000636874599598e0,
         0.000036583601823e0, -.000045036975204e0, -.000004537895710e0,
         0.000002937706971e0, 0.000000396566462e0, -.000000202099617e0,
         -.000000025267769e0, 0.000000013522610e0, 0.000000001229405e0,
         -.000000000821969e0, -.000000000050847e0, 0.000000000046333e0,
         0.000000000001943e0, -.000000000002563e0, -.000000000010991e0,
         0.628023320520852e0, 0.462762985338493e0, 0.150044637187895e0,
         -.028796057604906e0, -.036552745910311e0, -.001668087945272e0,
         0.006519840398744e0, 0.001130378079086e0, -.000887171310131e0,
         -.000242107641309e0, 0.000087294451594e0, 0.000034682122751e0,
         -.000004583768938e0, -.000003548684306e0, -.000000250453880e0,
         0.000000216991779e0, 0.000000080779570e0, 0.000000004558555e0,
         -.000000006944757e0, -.000000002849257e0, 0.000000000237816e0};
/*
      DOUBLE PRECISION AK, BK, DK
      COMMON / MX8 / AK(0:20,2),BK(0:22),DK(0:22)
*/
        double *AK;
	double Z;
	int II;
	int i;

	if (X <= 1.0e0) {
		II=1;
		Z=4.0e0*pow(X,0.2e0)-2.0e0;
		AK = &AKX[0];
	} else {
		II=2;
		Z=40.0e0*pow(X,-1.0e-1)/9.0e0-22.0e0/9.0e0;
		AK = &AKX[21];
	}
	for (i = 20; i >= 0; i--) {
		BK[i]=Z*BK[i+1]-BK[i+2]+AK[i];
		DK[i]=BK[i+1]+Z*DK[i+1]-DK[i+2];
	}
	return OK;
}
/* ---------------------------------------------------------------------- */
double G (double Y)
/* ---------------------------------------------------------------------- */
{
	return (2.0e0*(1.0e0-(1.0e0+Y)*exp(-Y))/(Y*Y));
}
/* ---------------------------------------------------------------------- */
double GP (double Y)
/* ---------------------------------------------------------------------- */
{
	return (-2.0e0*(1.0e0-(1.0e0+Y+Y*Y/2.0e0)*exp(-Y))/(Y*Y));
}
#ifdef SKIP
/* ---------------------------------------------------------------------- */
double ETHETA (double ZJ, double ZK, double I)
/* ---------------------------------------------------------------------- */
{
	double XCON, ZZ;
	double XJK, XJJ, XKK;

	if (ZJ == ZK) return(0.0);
	XCON=6.0e0*A0*sqrt(I);
	ZZ=ZJ*ZK;
/*
C
C     NEXT 3 ARE EQUATION (A1)
C
*/
	XJK=XCON*ZZ;
	XJJ=XCON*ZJ*ZJ;
	XKK=XCON*ZK*ZK;
/*
C
C     EQUATION (A2)
C
*/
	
	return (ZZ*(JAY(XJK)-JAY(XJJ)/2.0e0-JAY(XKK)/2.0e0)/(4.0e0*I));
}
/* ---------------------------------------------------------------------- */
double ETHETAP (double ZJ, double ZK, double I)
/* ---------------------------------------------------------------------- */
{
	double XCON, ZZ, ETHETA, ETHETAP;
	double XJK, XJJ, XKK;

	if (ZJ == ZK) return(0.0);
	XCON=6.0e0*A0*sqrt(I);
	ZZ=ZJ*ZK;
/*
C
C     NEXT 3 ARE EQUATION (A1)
C
*/
	XJK=XCON*ZZ;
	XJJ=XCON*ZJ*ZJ;
	XKK=XCON*ZK*ZK;
/*
C
C     EQUATION (A3)
C
*/
	ETHETA=ZZ*(JAY(XJK)-JAY(XJJ)/2.0e0-JAY(XKK)/2.0e0)/(4.0e0*I);
	ETHETAP=ZZ*(JPRIME(XJK)-JPRIME(XJJ)/2.0e0-JPRIME(XKK)/2.0e0)/(8.0e0*I*I) - ETHETA/I;
	return (ETHETAP);
}
#endif
/* ---------------------------------------------------------------------- */
int ETHETAS (double ZJ, double ZK, double I, double *etheta, double *ethetap)
/* ---------------------------------------------------------------------- */
{
	double XCON, ZZ;
	double XJK, XJJ, XKK;

	*etheta = 0.0;
	*ethetap = 0.0;
	if (ZJ == ZK) return(OK);
	XCON=6.0e0*A0*sqrt(I);
	ZZ=ZJ*ZK;
/*
C
C     NEXT 3 ARE EQUATION (A1)
C
*/
	XJK=XCON*ZZ;
	XJJ=XCON*ZJ*ZJ;
	XKK=XCON*ZK*ZK;
/*
C
C     EQUATION (A3)
C
*/
	*etheta=ZZ*(JAY(XJK)-JAY(XJJ)/2.0e0-JAY(XKK)/2.0e0)/(4.0e0*I);
	*ethetap=ZZ*(JPRIME(XJK)-JPRIME(XJJ)/2.0e0-JPRIME(XKK)/2.0e0)/(8.0e0*I*I) - *etheta/I;
	return (OK);
}
/* ---------------------------------------------------------------------- */
int pitzer_clean_up(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Free all allocated memory, except strings
 */
	int i;

	if (svnid == NULL) fprintf(stderr," ");
	for (i = 0; i < count_pitz_param; i++) {
		pitz_params[i] = (struct pitz_param *) free_check_null(pitz_params[i]);
	}
	pitz_params = (struct pitz_param **) free_check_null(pitz_params);
	for (i = 0; i < count_theta_param; i++) {
		theta_params[i] = (struct theta_param *) free_check_null(theta_params[i]);
	}
	theta_params = (struct theta_param **) free_check_null(theta_params);
	LGAMMA = (double *) free_check_null(LGAMMA);
	IPRSNT = (int *) free_check_null(IPRSNT);
	spec = (struct species **) free_check_null(spec);
	M = (double *) free_check_null(M);

	return OK;
}

/* ---------------------------------------------------------------------- */
int set_pz(int initial)
/* ---------------------------------------------------------------------- */
{
/*
 *   Sets initial guesses for unknowns if initial == TRUE
 *   Revises guesses whether initial is true or not
 */
	int i;
	struct solution *solution_ptr;
/*
 *   Set initial log concentrations to zero
 */
	iterations = -1;
	solution_ptr = use.solution_ptr;
	for (i=0; i < count_s_x; i++) {
		s_x[i]->lm = LOG_ZERO_MOLALITY;
		/*s_x[i]->lg = 0.0;*/
		s_x[i]->lg_pitzer = 0.0;
	}
/*
 *   Set master species activities
 */

	tc_x=solution_ptr->tc;
	tk_x=tc_x+273.15;
/*
 *   H+, e-, H2O
 */
	mass_water_aq_x = solution_ptr->mass_water;
	mu_x = solution_ptr->mu;
	s_h2o->moles = mass_water_aq_x/gfw_water;
	s_h2o->la = log10(solution_ptr->ah2o);
	AW = pow(10.0, s_h2o->la);
	s_hplus->la = - solution_ptr->ph;
	s_hplus->lm = s_hplus->la;
	s_hplus->moles = exp(s_hplus->lm * LOG_10)*mass_water_aq_x;
	s_eminus->la= - solution_ptr->solution_pe;
	if (initial == TRUE) initial_guesses(); 
	if (diffuse_layer_x == TRUE) initial_surface_water();
	revise_guesses();
	return(OK);
}
/* ---------------------------------------------------------------------- */
int initial_guesses(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Make initial guesses for activities of master species and
 *   ionic strength
 */
	int i;
	struct solution *solution_ptr;

	solution_ptr = use.solution_ptr;
	mu_x = s_hplus->moles + exp((solution_ptr->ph - 14.) * LOG_10) * mass_water_aq_x;
	mu_x /= mass_water_aq_x;
	s_h2o->la=0.0;
	for ( i=0; i < count_unknowns; i++ ) {
		if (x[i] == ph_unknown || x[i] == pe_unknown ) continue;
		if (x[i]->type < CB) {
			mu_x += x[i]->moles / mass_water_aq_x * 0.5 * x[i]->master[0]->s->z *
				x[i]->master[0]->s->z;
			x[i]->master[0]->s->la = log10(x[i]->moles/mass_water_aq_x);
		} else if (x[i]->type == CB) {
			x[i]->master[0]->s->la = log10(0.001 * x[i]->moles/mass_water_aq_x);
		} else if (x[i]->type == SOLUTION_PHASE_BOUNDARY) {
			x[i]->master[0]->s->la = log10(0.001 * x[i]->moles/mass_water_aq_x);
		} else if (x[i]->type == EXCH) {
			if (x[i]->moles <= 0) {
				x[i]->master[0]->s->la = MIN_RELATED_LOG_ACTIVITY;
			} else {
				x[i]->master[0]->s->la = log10(x[i]->moles);
			}
		} else if (x[i]->type == SURFACE) {
			if (x[i]->moles <= 0) {
				x[i]->master[0]->s->la = MIN_RELATED_LOG_ACTIVITY;
			} else {
				x[i]->master[0]->s->la = log10(0.1 * x[i]->moles);
			}
		} else if (x[i]->type == SURFACE_CB) {
			x[i]->master[0]->s->la = 0.0;
		}
	}
	return(OK);
}
/* ---------------------------------------------------------------------- */
int revise_guesses(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Revise molalities species
 */
	int i;
	int iter, max_iter, repeat, fail;
	LDBLE weight, f;

	max_iter = 10;
	/* gammas(mu_x);*/
	iter = 0;
	repeat = TRUE;
 	fail = FALSE;;
	while ( repeat == TRUE ) {
		iter++;
		if (debug_set == TRUE) {
			output_msg(OUTPUT_MESSAGE,"\nBeginning set iteration %d.\n", iter);
		}			
 		if (iter == max_iter + 1) {
			output_msg(OUTPUT_LOG, "Did not converge in set, iteration %d.\n", iterations);
 			fail = TRUE;
 		}
 		if (iter > 2*max_iter) {
			output_msg(OUTPUT_LOG, "Did not converge with relaxed criteria in set.\n");
 			return(OK);
  		}
		molalities(TRUE);
		/*pitzer();*/
		/*s_h2o->la = 0.0;*/
		/*molalities(TRUE);*/
		mb_sums();
		if (state < REACTION) {
			sum_species();
		} else {
			for (i = 0; i < count_unknowns; i++)  {
				x[i]->sum = x[i]->f;
			}
		}
		/*n
		if (debug_set == TRUE) {
			pr.species = TRUE;
			pr.all = TRUE;
			print_species();
		}
		*/
		repeat=FALSE;
		for ( i=0; i < count_unknowns; i++ ) {
			if (x[i] == ph_unknown || x[i] == pe_unknown) continue;
			if (x[i]->type == MB || 
/*			    x[i]->type == ALK || */
			    x[i]->type == CB || 
			    x[i]->type == SOLUTION_PHASE_BOUNDARY || 
			    x[i]->type == EXCH || 
			    x[i]->type == SURFACE ) {
				
				if ( debug_set == TRUE ) {
					output_msg(OUTPUT_MESSAGE,"\n\t%5s  at beginning of set %d: %e\t%e\t%e\n", x[i]->description, iter, (double) x[i]->sum, (double) x[i]->moles, (double) x[i]->master[0]->s->la);
				}
				if (fabs(x[i]->moles) < 1e-30) x[i]->moles = 0;
				f = fabs(x[i]->sum);
				if (f == 0 && x[i]->moles == 0) {
					x[i]->master[0]->s->la = MIN_RELATED_LOG_ACTIVITY;
					continue;
				} else if (f == 0) {
					repeat = TRUE;
					x[i]->master[0]->s->la += 5;
/*!!!!*/				if (x[i]->master[0]->s->la < -999.) x[i]->master[0]->s->la = MIN_RELATED_LOG_ACTIVITY;
 				} else if (fail == TRUE && f < 1.5 * fabs(x[i]->moles)) {
 					continue;
				} else if (f > 1.5 * fabs(x[i]->moles) || f < 1e-5 * fabs(x[i]->moles) ) {
					weight = (f < 1e-5 * fabs(x[i]->moles)) ? 0.3 : 1.0;
					if (x[i]->moles <= 0) {
						x[i]->master[0]->s->la = MIN_RELATED_LOG_ACTIVITY;
					} else {
						repeat = TRUE;
						x[i]->master[0]->s->la += weight * log10(fabs(x[i]->moles / x[i]->sum));
					}
					if ( debug_set == TRUE ) {
						output_msg(OUTPUT_MESSAGE,"\t%5s not converged in set %d: %e\t%e\t%e\n", x[i]->description, iter, (double) x[i]->sum, (double) x[i]->moles, (double) x[i]->master[0]->s->la);
					}
				}
			} else if (x[i]->type == ALK) {
				f = total_co2;
 				if (fail == TRUE && f < 1.5 * fabs(x[i]->moles)) {
 					continue;
 				}
				if (f > 1.5 * fabs(x[i]->moles) || f < 1e-5 * fabs(x[i]->moles) ) {
					repeat = TRUE;
					weight = (f < 1e-5 * fabs(x[i]->moles)) ? 0.3 : 1.0;
					x[i]->master[0]->s->la += weight * 
						log10(fabs(x[i]->moles / x[i]->sum));
					if ( debug_set == TRUE ) {
						output_msg(OUTPUT_MESSAGE,"%s not converged in set. %e\t%e\t%e\n", x[i]->description, (double) x[i]->sum, (double) x[i]->moles, (double) x[i]->master[0]->s->la);
					}
				}
			}
		}
	}
	output_msg(OUTPUT_LOG,"Iterations in revise_guesses: %d\n", iter);
	/*mu_x = mu_unknown->f * 0.5 / mass_water_aq_x;*/
	if (mu_x <= 1e-8) {
		mu_x = 1e-8;
	}
	/*gammas(mu_x);*/
	return(OK);
}
/* ---------------------------------------------------------------------- */
int jacobian_pz(void)
/* ---------------------------------------------------------------------- */
{
	double *base;
	double d, d1, d2;
	int i, j;

	if (full_pitzer == TRUE) {
		molalities(TRUE);
		pitzer();
		residuals();
	}
	base = (LDBLE *) PHRQ_malloc((size_t) count_unknowns * sizeof(LDBLE));
	if (base == NULL) malloc_error();
	for (i = 0; i < count_unknowns; i++) {
		base[i] = residual[i];
	}
	d = 0.0001;
	d1 = d*log(10.0);
	d2 = 0;
	for (i = 0; i < count_unknowns; i++) {
		switch (x[i]->type) {
		case MB:
		case ALK:
		case CB:
		case SOLUTION_PHASE_BOUNDARY:
		case EXCH:
		case SURFACE:
		case SURFACE_CB:
			x[i]->master[0]->s->la += d;
			d2 = d1;
			break;
		case AH2O:
			x[i]->master[0]->s->la += d;
			d2 = d1;
			break;
		case PITZER_GAMMA:
			x[i]->s->lg += d;
			d2 = d;
			break;
		case MH2O:
			mass_water_aq_x *= (1.0 + d);
			x[i]->master[0]->s->moles = mass_water_aq_x/gfw_water;
			d2 = log(1.0 + d);
			break;
		case MU:
		case MH:
		case PP:
		case S_S_MOLES:
			continue;
			break;
		}
		molalities(TRUE);
		if (full_pitzer == TRUE) pitzer();
		mb_sums();
		residuals();
		for (j = 0; j < count_unknowns; j++) {
			array[j*(count_unknowns + 1) + i] = -(residual[j] - base[j])/d2;
		}
		switch (x[i]->type) {
		case MB:
		case ALK:
		case CB:
		case SOLUTION_PHASE_BOUNDARY:
		case EXCH:
		case SURFACE:
		case SURFACE_CB:
		case AH2O:
			x[i]->master[0]->s->la -= d;
			break;
		case PITZER_GAMMA:
			x[i]->s->lg -= d;
			break;
		case MH2O:
			mass_water_aq_x /= (1 + d);
			x[i]->master[0]->s->moles = mass_water_aq_x/gfw_water;
			break;
		}
	}
	molalities(TRUE);
	if (full_pitzer == TRUE) pitzer();
	mb_sums();
	residuals();
	free_check_null(base);
	return OK;
}
/* ---------------------------------------------------------------------- */
int model_pz(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   model is called after the equations have been set up by prep
 *   and initial guesses have been made in set.
 * 
 *   Here is the outline of the calculation sequence:
 *      residuals--residuals are calculated, if small we are done
 *      sum_jacobian--jacobian is calculated 
 *      ineq--inequality solver is called
 *      reset--estimates of unknowns revised, if changes are small solution
 *         has been found, usually convergence is found in residuals.
 *      gammas--new activity coefficients
 *      molalities--calculate molalities
 *      mb_sums--calculate mass-balance sums
 *      mb_gases--decide if gas_phase exists
 *      mb_s_s--decide if solid_solutions exists
 *      switch_bases--check to see if new basis species is needed
 *         reprep--rewrite equations with new basis species if needed
 *         revise_guesses--revise unknowns to get initial mole balance
 *      check_residuals--check convergence one last time
 *         sum_species--calculate sums of elements from species concentrations
 *
 *      An additional pass through may be needed if unstable phases still exist
 *         in the phase assemblage. 
 */
	int kode, return_kode;
	int r;
	int count_infeasible, count_basis_change;
	int debug_model_save;
	int mass_water_switch_save;
	if (svnid == NULL) fprintf(stderr," ");

/*	debug_model = TRUE; */
/*	debug_prep = TRUE; */
/*	debug_set = TRUE; */
	/* mass_water_switch == TRUE, mass of water is constant */
	mass_water_switch_save = mass_water_switch;
	if (mass_water_switch_save == FALSE && delay_mass_water == TRUE) {
		mass_water_switch = TRUE;
	}
	debug_model_save = debug_model;
	pe_step_size_now = pe_step_size;
	step_size_now = step_size;
	status(0, NULL);
	iterations=0;
	gamma_iterations = 0;
	count_basis_change = count_infeasible = 0;
	stop_program = FALSE;
	remove_unstable_phases = FALSE;
	if (always_full_pitzer == TRUE) {
		full_pitzer = TRUE;
	} else {
		full_pitzer = FALSE;
	}
	for (; ; ) {
		mb_gases();
		mb_s_s();
		kode = 1;
		while ( ( r = residuals() ) != CONVERGED || remove_unstable_phases == TRUE) {
#if defined(PHREEQCI_GUI)
			if (WaitForSingleObject(g_hKill /*g_eventKill*/, 0) == WAIT_OBJECT_0)
				{
					error_msg("Execution canceled by user.", CONTINUE);
					RaiseException(USER_CANCELED_RUN, 0, 0, NULL);
				}
#endif
			iterations++;
			if (iterations > itmax - 1 && debug_model == FALSE && pr.logfile == TRUE) {
				set_forward_output_to_log(TRUE);
				debug_model = TRUE;
			}
			if (debug_model == TRUE) {
				output_msg(OUTPUT_MESSAGE,"\nIteration %d\tStep_size = %f\n", 
					   iterations, (double) step_size_now);
				output_msg(OUTPUT_MESSAGE,"\t\tPe_step_size = %f\n\n", (double) pe_step_size_now);
			}
			/*
			 *   Iterations exceeded
			 */
			if (iterations > itmax ) {
				sprintf(error_string,"Maximum iterations exceeded, %d\n",itmax);
				warning_msg(error_string);
				stop_program = TRUE;
				break;
			}
			/*
			 *   Calculate jacobian
			 */
			gammas_pz();
			jacobian_sums();
			jacobian_pz();
			/*
			 *   Full matrix with pure phases
			 */
			if ( r == OK || remove_unstable_phases == TRUE) {
				return_kode = ineq(kode);
				if ( return_kode != OK ) {
					if (debug_model == TRUE) {
						output_msg(OUTPUT_MESSAGE, "Ineq had infeasible solution, "
							   "kode %d, iteration %d\n", 
							   return_kode, iterations);
					}
					output_msg(OUTPUT_LOG, "Ineq had infeasible solution, "
						   "kode %d, iteration %d\n", return_kode, iterations);
					count_infeasible++;
				}
				if ( return_kode == 2 ) { 
					ineq(0);
				}
				reset();
			}
			gammas_pz();
			if (full_pitzer == TRUE) pitzer();
			if (always_full_pitzer == TRUE) {
				full_pitzer = TRUE;
			} else {
				full_pitzer = FALSE;
			}
			molalities(TRUE);
			if(use.surface_ptr != NULL && 
			   use.surface_ptr->diffuse_layer == TRUE &&
			   use.surface_ptr->related_phases == TRUE)
				initial_surface_water();
			mb_sums();
			mb_gases();
			mb_s_s();
			/* debug
			   species_list_sort();
			   sum_species();
			   print_species();
			   print_exchange();
			   print_surface();
			*/
			if (stop_program == TRUE) {
				break;
			}
		}
/*
 *   Check for stop_program
 */

		if (stop_program == TRUE) {
			break;
		}
		if (check_residuals() == ERROR) {
			stop_program = TRUE;
			break;
		}
		if (remove_unstable_phases == FALSE && mass_water_switch_save == FALSE &&
		    mass_water_switch == TRUE) {
			output_msg(OUTPUT_LOG,"\nChanging water switch to FALSE. Iteration %d.\n", iterations);
			mass_water_switch = FALSE;
			continue;
		}
		gamma_iterations++;
		if (gamma_iterations > itmax ) {
			sprintf(error_string,"Maximum gamma iterations exceeded, %d\n", itmax);
			warning_msg(error_string);
			stop_program = TRUE;
			break;
		}
		if (check_gammas_pz() != TRUE) {
			full_pitzer = TRUE;
			continue;
		}
		if (remove_unstable_phases == FALSE) break;
		if (debug_model == TRUE) {
			output_msg(OUTPUT_MESSAGE,"\nRemoving unstable phases. Iteration %d.\n", iterations);
		}
		output_msg(OUTPUT_LOG,"\nRemoving unstable phases. Iteration %d.\n", iterations);
        }
	output_msg(OUTPUT_LOG,"\nNumber of infeasible solutions: %d\n",count_infeasible);
	output_msg(OUTPUT_LOG,"Number of basis changes: %d\n\n",count_basis_change);
	output_msg(OUTPUT_LOG,"Number of iterations: %d\n\n", iterations);
	debug_model = debug_model_save;
	set_forward_output_to_log(FALSE);
	if (stop_program == TRUE) {
		return(ERROR);
	}
	return(OK);
}
/* ---------------------------------------------------------------------- */
int check_gammas_pz(void)
/* ---------------------------------------------------------------------- */
{
	double old_aw, old_mu, tol;
	int converge, i;

	old_mu = mu_x;
	old_aw = s_h2o->la;
	pitzer();
	molalities(TRUE);
	mb_sums();
	converge = TRUE;
	tol = convergence_tolerance*10.;
	for (i = 0; i < count_unknowns; i++) {
		if (x[i]->type != PITZER_GAMMA) continue;
		if (fabs(x[i]->s->lg - x[i]->s->lg_pitzer) > tol) {
			converge = FALSE;
		}
	}
	if (fabs(old_mu - mu_x) > tol) converge = FALSE;
	if ((pow(10.0,s_h2o->la) - AW) > tol) converge = FALSE;
	return converge;
}
/* ---------------------------------------------------------------------- */
int gammas_pz ()
/* ---------------------------------------------------------------------- */
{
/*
 *   Need exchange gammas for pitzer
 */
	int i, j;
	double coef;
	/* Initialize */
/*
 *   Calculate activity coefficients
 */
	for (i=0; i < count_s_x; i++) {
		switch (s_x[i]->gflag) {
		    case 0:                   /* uncharged */
		    case 1:                   /* Davies */
		    case 2:                   /* Extended D-H, WATEQ D-H */
		    case 3:                   /* Always 1.0 */
			    break;
		    case 4:		      /* Exchange */
			    /* Now calculated in next loop */
			    break;
		    case 5:                   /* Always 1.0 */
			    break;
		    case 6:		      /* Surface */
/*
 *   Find moles of sites. 
 *   s_x[i]->equiv is stoichiometric coefficient of sites in species
 */
			for (j=1; s_x[i]->rxn_x->token[j].s != NULL; j++) {
				if (s_x[i]->rxn_x->token[j].s->type == SURF) {
					s_x[i]->alk = s_x[i]->rxn_x->token[j].s->primary->unknown->moles;
					break;
				}
			}
			if (s_x[i]->alk > 0) {
				s_x[i]->lg = log10(s_x[i]->equiv / s_x[i]->alk);
				s_x[i]->dg = 0.0;
			} else {
				s_x[i]->lg = 0.0;
				s_x[i]->dg = 0.0;
			}
			break;
		    case 7:		   /* LLNL */
			    break;
		    case 8:		   /* LLNL CO2*/
			    break;
		    case 9:		   /* activity water */
			s_x[i]->lg = log10(exp( s_h2o->la * LOG_10) * gfw_water);
			s_x[i]->dg = 0.0;
			break;
		}
/*
		if (mu_unknown != NULL) {
			if (fabs(residual[mu_unknown->number]) > 0.1 &&
			    fabs(residual[mu_unknown->number])/mu_x > 0.5) {
				s_x[i]->dg = 0.0;
			}
		}
 */
	}
	/*
	 *  calculate exchange gammas 
	 */

	if (use.exchange_ptr != NULL) {
		for (i=0; i < count_s_x; i++) {
			switch (s_x[i]->gflag) {
			case 0:               /* uncharged */
			case 1:               /* Davies */
			case 2:               /* Extended D-H, WATEQ D-H */
			case 3:               /* Always 1.0 */
			case 5:               /* Always 1.0 */
			case 6:		      /* Surface */
			case 7:		      /* LLNL */
			case 8:		      /* LLNL CO2*/
			case 9:		      /* activity water */
				break;
			case 4:		      /* Exchange */

				/*
				 *   Find CEC
				 *   z contains valence of cation for exchange species, alk contains cec
				 */
				/* !!!!! */
				for (j=1; s_x[i]->rxn_x->token[j].s != NULL; j++) {
					if (s_x[i]->rxn_x->token[j].s->type == EX) {
						s_x[i]->alk = s_x[i]->rxn_x->token[j].s->primary->unknown->moles;
						break;
					}
				}
				/*
				 *   Master species is a dummy variable with meaningless activity and mass
				 */
				s_x[i]->lg = 0.0;
				s_x[i]->dg = 0.0;
				if (s_x[i]->primary != NULL) {
					break;
				} 
				/*
				 *   All other species
				 */

				/* modific 29 july 2005... */
				if (s_x[i]->equiv != 0 && s_x[i]->alk > 0) {
					s_x[i]->lg = log10(fabs(s_x[i]->equiv) / s_x[i]->alk);
				}
				if (use.exchange_ptr->pitzer_exchange_gammas == TRUE) {
					/* Assume equal gamma's of solute and exchangeable species...  */
					for (j=1; s_x[i]->rxn_x->token[j].s != NULL; j++) {
						if (s_x[i]->rxn_x->token[j].s->type == EX) continue;
						coef = s_x[i]->rxn_x->token[j].coef;
						s_x[i]->lg += coef * s_x[i]->rxn_x->token[j].s->lg;
						s_x[i]->dg += coef * s_x[i]->rxn_x->token[j].s->dg;
					}
				}
			}
		}
	}
/* ...end modific 29 july 2005 */

	return(OK);
}
