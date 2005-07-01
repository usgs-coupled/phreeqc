#define EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "output.h"
#include "phrqproto.h"
#define PITZER
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
double COSMOT;

/* routines */
static int calc_pitz_param (struct pitz_param *pz_ptr, double TK, double TR);
static int ISPEC(char * name);
/*static int DH_AB (double TK, double *A, double *B);*/
static double G (double Y);
static double GP (double Y);
static double ETHETAP (double ZJ, double ZK, double I);
static double ETHETA (double ZJ, double ZK, double I);
static int BDK (double X);

/* ---------------------------------------------------------------------- */
int pitzer_init (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Initialization for pitzer
 */
	pitzer_model = FALSE;
	max_pitz_param = 100;
	count_pitz_param = 0;
	space ((void *) &pitz_params, INIT, &max_pitz_param, sizeof(struct pitz_param *));
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
	/*
	 *  allocate pointers to species structures
	 */
	if (spec != NULL) spec = free_check_null(spec);
	spec = PHRQ_malloc((size_t) (3*count_s*sizeof(struct species *)));
	if (spec == NULL) malloc_error();
	for (i = 1; i < 3*count_s; i++) spec[i] = NULL;
	cations = spec;
	neutrals = &(spec[count_s]);
	anions = &(spec[2*count_s]);
	MAXCATIONS = count_s;
	FIRSTANION = 2*count_s;
	MAXNEUTRAL = count_s;
	count_cations = 0;
	count_anions = 0;
	count_neutrals = 0;
	/*
	 *  allocate other arrays for Pitzer
	 */
	if (IPRSNT != NULL) IPRSNT = free_check_null(IPRSNT);
	IPRSNT = PHRQ_malloc((size_t) (3*count_s*sizeof(int)));
	if (IPRSNT == NULL) malloc_error();
	if (M != NULL) M = free_check_null(M);
	M = PHRQ_malloc((size_t) (3*count_s*sizeof(double)));
	if (M == NULL) malloc_error();
	if (LGAMMA != NULL) LGAMMA = free_check_null(LGAMMA);
	LGAMMA = PHRQ_malloc((size_t) (3*count_s*sizeof(double)));
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
	 */
	for (i = 0; i < count_cations - 1; i++) {
		for (j = i+1; j < count_cations; j++) {
			sprintf(line,"%s %s 1", spec[i]->name, spec[j]->name);
			pzp_ptr = pitz_param_read(line, 2);
			pzp_ptr->type = TYPE_ETHETA;
			if (count_pitz_param >= max_pitz_param) {
				space ((void *) &pitz_params, count_pitz_param, &max_pitz_param, sizeof(struct pitz_param *));
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
				space ((void *) &pitz_params, count_pitz_param, &max_pitz_param, sizeof(struct pitz_param *));
			}
			pitz_params[count_pitz_param++] = pzp_ptr;
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
				pitz_params[i]->alpha = 1.0;
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
		"psi"                    /* 7 */
	};
	int count_opt_list = 8;
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
						    space ((void *) &pitz_params, count_pitz_param, &max_pitz_param, sizeof(struct pitz_param *));
					    }
					    
					    pitz_params[count_pitz_param] = pzp_ptr;
					    count_pitz_param++;
				    } else {
					    pitz_params[j] = free_check_null(pitz_params[j]);
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
	static double OTEMP=0.0;

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
		A0=1.400684e6*(DW0/pow(pow((DC0*TK),3.0e0),0.5e0));
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
int pitzer ()
/* ---------------------------------------------------------------------- */
{
	int i, i0, i1, i2;
	double param, alpha, z0, z1, z2;
	/*
	double CONV, XI, XX, OSUM, BIGZ, DI, F, XXX, GAMCLM, 
		CSUM, PHIMAC, OSMOT, BMXP, ETHEAP, CMX, BMX, PHI,
		BMXPHI, PHIPHI, AW, A, B;
	*/
	double CONV, XI, XX, OSUM, BIGZ, DI, F, XXX, GAMCLM, 
		CSUM, PHIMAC, OSMOT, AW, B;
	double COSMOT;
	double I, TK;
	int IC, ICON;
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
		if (spec[i] != NULL && spec[i]->in == TRUE) {
			M[i] = under(spec[i]->lm);
			if (M[i] > MIN_TOTAL) IPRSNT[i] = TRUE;
		}
	}
	for (i = count_s; i < count_s + count_neutrals; i++) {
		if (M[i] > MIN_TOTAL) LNEUT = TRUE;
	}
#ifdef SKIP
	/*  TESTING !!!!!!!!!!! */
	for (i = 0; i < 3*count_s; i++) {
		M[i] = 0.0;
		IPRSNT[i] = FALSE;
	}
	/*
	 * 67 Cl-
	 * 68 HCO3-
	 * 69 HSO4-
	 * 70 OH-
	 * 71 SO4-2
	 */
	M[1] = 1.0; /* Ca+2 */
	IPRSNT[1] = TRUE;
	M[11] = 0.0; /* Na+ */
	IPRSNT[11] = FALSE;
	M[5] = 1.0; /* K+ */
	IPRSNT[5] = TRUE;
	M[2*count_s + 5] = 1.0; /* Cl */
	IPRSNT[2*count_s + 5] = TRUE;
	IC = 2*count_s + 5;
	M[71] = 1.0; /* SO4-2 */
	IPRSNT[71] = TRUE;
#endif
	ICON = 0;
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
	GAMCLM=F+I*2.0e0*(mcb0->U.b0 + mcb1->U.b1*XXX) + 1.5e0*mcc0->U.c0*I*I;
	CSUM=0.0e0;
	OSMOT=-(A0)*pow(I,1.5e0)/(1.0e0+B*DI);
/*
 *  Sums for F, LGAMMA, and OSMOT
 */
	for (i = 0; i < count_pitz_param; i++) {
		if (i == 107) {
			fprintf(stderr,"\n");
		}
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
			F += M[i0]*M[i1]*param*GP(alpha*sqrt(I))/I; 
			LGAMMA[i0] += M[i1]*2.0*param*G(alpha*sqrt(I));
			LGAMMA[i1] += M[i0]*2.0*param*G(alpha*sqrt(I));
			OSMOT += M[i0]*M[i1]*param*exp(-alpha*DI);
			break;
		case TYPE_B2:
			F += M[i0]*M[i1]*param*GP(alpha*sqrt(I))/I; 
			LGAMMA[i0] += M[i1]*2.0*param*G(alpha*sqrt(I));
			LGAMMA[i1] += M[i0]*2.0*param*G(alpha*sqrt(I));
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
			F += M[i0]*M[i1]*ETHETAP(z0, z1, I);
			LGAMMA[i0] += 2.0*M[i1]*(ETHETA(z0, z1, I) ); 
			LGAMMA[i1] += 2.0*M[i0]*(ETHETA(z0, z1, I) ); 
			OSMOT += M[i0]*M[i1]*(ETHETA(z0, z1, I) + I*ETHETAP(z0, z1, I) ); 
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
      if (ICON != 0) {
	      PHIMAC=LGAMMA[IC]-GAMCLM;
/*
C
C     CORRECTED ERROR IN PHIMAC, NOVEMBER, 1989
C
*/
	      for (i = 0; i < count_cations; i++) {
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
      s_h2o->la=log10(AW);
      mu_x = I;
      for (i = 0; i < 2*count_s + count_anions; i++) {
	      if (IPRSNT[i] == FALSE) continue;
	      spec[i]->lg=LGAMMA[i]*CONV;
	      /*n
	      output_msg(OUTPUT_MESSAGE, "%s:\t%e\t%e\t%e \n", spec[i]->name, spec[i]->la, spec[i]->lm, spec[i]->lg);
	      */
      }
      /*
      output_msg(OUTPUT_MESSAGE, "OSMOT: %e\n", OSMOT);
      output_msg(OUTPUT_MESSAGE, "COSMOT: %e\n", COSMOT);
      output_msg(OUTPUT_MESSAGE, "F: %e\n", F);
      */
      /*
      *I_X = I;
      *COSMOT_X = COSMOT;
      */
      return(OK);
}
#ifdef SKIP
/* ---------------------------------------------------------------------- */
int pitzer (double *TK_X, double *I_X, double *COSMOT_X)
/* ---------------------------------------------------------------------- */
{
	int i, j, i0, i1, i2;
	double param, alpha, z0, z1, z2;
	/*
	double CONV, XI, XX, OSUM, BIGZ, DI, F, XXX, GAMCLM, 
		CSUM, PHIMAC, OSMOT, BMXP, ETHEAP, CMX, BMX, PHI,
		BMXPHI, PHIPHI, AW, A, B;
	*/
	double CONV, XI, XX, OSUM, BIGZ, DI, F, XXX, GAMCLM, 
		CSUM, PHIMAC, OSMOT, AW, B;
	double COSMOT;
	double I, TK;
	int IC, ICON;
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
	I = *I_X;
	TK = *TK_X;
	/*	DH_AB(TK, &A, &B); */
	/*
	  C
	  C     TRANSFER DATA FROM MO TO M
	  C
	*/
	for (i = 0; i < 3*count_s; i++) {
		IPRSNT[i] = FALSE;
		if (spec[i] != NULL) {
			M[i] = under(spec[i]->lm);
			if (M[i] > MIN_TOTAL) IPRSNT[i] = TRUE;
		}
	}
	for (i = count_s; i < count_s + count_neutrals; i++) {
		if (M[i] > MIN_TOTAL) LNEUT = TRUE;
	}
	/*  TESTING !!!!!!!!!!! */
	for (i = 0; i < 3*count_s; i++) {
		M[i] = 0.0;
		IPRSNT[i] = FALSE;
	}
	/*
	 * 67 Cl-
	 * 68 HCO3-
	 * 69 HSO4-
	 * 70 OH-
	 * 71 SO4-2
	 */
	M[1] = 1.0; /* Ca+2 */
	IPRSNT[1] = TRUE;
	M[11] = 0.0; /* Na+ */
	IPRSNT[11] = FALSE;
	M[5] = 1.0; /* K+ */
	IPRSNT[5] = TRUE;
	M[2*count_s + 5] = 1.0; /* Cl */
	IPRSNT[2*count_s + 5] = TRUE;
	IC = 2*count_s + 5;
	M[71] = 1.0; /* SO4-2 */
	IPRSNT[71] = TRUE;
	ICON = 0;
/*
C
C     COMPUTE PITZER COEFFICIENTS' TEMPERATURE DEPENDENCE
C
*/
	PTEMP(TK);
	for (i = 0; i < 2*count_s + count_anions; i++) {
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
	GAMCLM=F+I*2.0e0*(mcb0->U.b0 + mcb1->U.b1*XXX) + 1.5e0*mcc0->U.c0*I*I;
/*
 *
C
C     EQUATION (3) PART 1
C
 */
	/* F=F+M(J)*M(K)*BMXP() */
	/* BMXP=(BCX(2,J,K)*GP(ALPHA(2)*DSQRT(I))+BCX(3,J,K)*GP(ALPHA(3)*DSQRT(I)))/I */
	for (i = 0; i < count_pitz_param; i++) {
		i0 = pitz_params[i]->ispec[0];
		i1 = pitz_params[i]->ispec[1];
		if (IPRSNT[i0] == FALSE || IPRSNT[i1] == FALSE) continue;
		param = pitz_params[i]->p;
		alpha = pitz_params[i]->alpha;
		if (pitz_params[i]->type == TYPE_B1) {
			F += M[i0]*M[i1]*param*GP(alpha*sqrt(I))/I; 
		} else if (pitz_params[i]->type == TYPE_B2) {
			F += M[i0]*M[i1]*param*GP(alpha*sqrt(I))/I; 
		}
	}
/*
C
C     EQUATION (3) PART 2
C
*/
	/* F=F+M(J)*M(K)*ETHEAP() */
	for (i = 0; i < count_cations - 1; i++) {
		z0 = cations[i]->z;
		if (IPRSNT[i] == FALSE) continue;
		for (j = i+1; j < count_cations; j++) {
			if (IPRSNT[j] == FALSE) continue;
			z1 = cations[j]->z;
			F=F+M[i]*M[j]*ETHETAP(z0, z1, I);
		}
	}
	for (i = 2*count_s; i < 2*count_s + count_anions - 1; i++) {
		z0 = spec[i]->z;
		for (j = i+1; j < 2*count_s + count_anions; j++) {
			z1 = spec[j]->z;
			F=F+M[i]*M[j]*ETHETAP(z0, z1, I);
		}
	}

/*
C
C     EQUATION (2B) PART 4
C
*/
	CSUM=0.0e0;
	for (i = 0; i < count_pitz_param; i++) {
		/* CSUM=CSUM+M(J)*M(K)*CMX()*/
		/* CMX=BCX(4,J,K)/(2.0D0*DSQRT(DABS(Z(J)*Z(K)))) */
		if (pitz_params[i]->type == TYPE_C0) {
			i0 = pitz_params[i]->ispec[0];
			i1 = pitz_params[i]->ispec[1];
			z0 = spec[pitz_params[i]->ispec[0]]->z;
			z1 = spec[pitz_params[i]->ispec[1]]->z;
			CSUM += M[i0]*M[i1]*pitz_params[i]->p/(2.0e0*sqrt(fabs(z0*z1)));
		}
	}
	for (i = 0; i < count_cations; i++) {
		/*
		  C
		  C     EQUATION (2C) PART 1
		  C
		  LGAMMA(K)=Z(K)*Z(K)*F+DABS(Z(K))*CSUM
		*/
		z0 = spec[i]->z;
		LGAMMA[i]=z0*z0*F+fabs(z0)*CSUM;
	}
	for (i = 2*count_s; i < 2*count_s + count_anions; i++) {
		/*
		  C
		  C     EQUATION (2C) PART 1
		  C
		  LGAMMA(K)=Z(K)*Z(K)*F+DABS(Z(K))*CSUM
		*/
		z0 = spec[i]->z;
		LGAMMA[i]=z0*z0*F+fabs(z0)*CSUM;
	}

	/*
	 *  Add etheta term to activity coefficients
	 */
	for (i = 0; i < count_cations - 1; i++) {
		z0 = cations[i]->z;
		if (IPRSNT[i] == FALSE) continue;
		for (j = i+1; j < count_cations; j++) {
			if (IPRSNT[j] == FALSE) continue;
			z1 = cations[j]->z;
			LGAMMA[i] += 2.0*M[j]*(ETHETA(z0, z1, I) ); 
			LGAMMA[j] += 2.0*M[i]*(ETHETA(z0, z1, I) ); 
		}
	}
	for (i = 2*count_s; i < 2*count_s + count_anions - 1; i++) {
		z0 = spec[i]->z;
		for (j = i+1; j < 2*count_s + count_anions; j++) {
			z1 = spec[j]->z;
			LGAMMA[i] += 2.0*M[j]*(ETHETA(z0, z1, I) ); 
			LGAMMA[j] += 2.0*M[i]*(ETHETA(z0, z1, I) ); 
		}
	}
/*
 *    Sums of Pitzer parameters for activity coefficients
 */
	for (i = 0; i < count_pitz_param; i++) {
		i0 = pitz_params[i]->ispec[0];
		i1 = pitz_params[i]->ispec[1];
		if (IPRSNT[i0] == FALSE || IPRSNT[i1] == FALSE) continue;
		z0 = spec[pitz_params[i]->ispec[0]]->z;
		z1 = spec[pitz_params[i]->ispec[1]]->z;
		param = pitz_params[i]->p;
		alpha = pitz_params[i]->alpha;
		/*
		  C
		  C     EQUATION (2B) PART 3
		  C
		  LGAMMA(J)=LGAMMA(J)+M(K)*(2.0D0*BMX()+BIGZ*CMX())
		  BMX=BCX(1,J,K)+BCX(2,J,K)*G(ALPHA(2)*DSQRT(I))+BCX(3,J,K)*G(ALPHA(3)*DSQRT(I))
		  CMX=BCX(4,J,K)/(2.0D0*DSQRT(DABS(Z(J)*Z(K)))) 
		*/
		if (pitz_params[i]->type == TYPE_B0) {
			LGAMMA[i0] += M[i1]*2.0*param;
			LGAMMA[i1] += M[i0]*2.0*param;
		} else if (pitz_params[i]->type == TYPE_B1) {
			LGAMMA[i0] += M[i1]*2.0*param*G(alpha*sqrt(I));
			LGAMMA[i1] += M[i0]*2.0*param*G(alpha*sqrt(I));
		} else if (pitz_params[i]->type == TYPE_B2) {
			LGAMMA[i0] += M[i1]*2.0*param*G(alpha*sqrt(I));
			LGAMMA[i1] += M[i0]*2.0*param*G(alpha*sqrt(I));
		} else if (pitz_params[i]->type == TYPE_C0) {
			LGAMMA[i0] += M[i1]*BIGZ*param/(2.0*sqrt(fabs(z0*z1))); 
			LGAMMA[i1] += M[i0]*BIGZ*param/(2.0*sqrt(fabs(z0*z1))); 
		} else if (pitz_params[i]->type == TYPE_THETA) {
			/*
			  C
			  C     EQUATION (2B) PART 2
			  C
			  LGAMMA(J)=LGAMMA(J)+2.0D0*M(K)*PHI()
			  PHI=THETA(J,K)+ETHETA()
			  etheta added above in case theta is not defined for a pair
			*/
			LGAMMA[i0] += 2.0*M[i1]*(param /*+ ETHETA(z0, z1, I) */ ); 
			LGAMMA[i1] += 2.0*M[i0]*(param /*+ ETHETA(z0, z1, I) */ ); 
		} else if (pitz_params[i]->type == TYPE_PSI) {
			i2 = pitz_params[i]->ispec[2];
			if (IPRSNT[i2] == FALSE) continue;
			z2 = spec[pitz_params[i]->ispec[2]]->z;
			/*
			  C
			  C     EQUATION (2B) PART 2
			  C
			  LGAMMA(J)=LGAMMA(J)+M(KK)*M(K)*PSI(J,K,KK)
			*/
			LGAMMA[i0] += M[i1]*M[i2]*param;
			LGAMMA[i1] += M[i0]*M[i2]*param;
			LGAMMA[i2] += M[i0]*M[i1]*param;

		} else if (pitz_params[i]->type == TYPE_LAMDA) {
			/*LGAMMA(J)=LGAMMA(J)+2.0D0*MN(K)*LAM(J,K)*/
			LGAMMA[i0] += 2.0*M[i1]*param;
			LGAMMA[i1] += 2.0*M[i0]*param;
		} else if (pitz_params[i]->type == TYPE_ZETA) {
			i2 = pitz_params[i]->ispec[2];
			if (IPRSNT[i2] == FALSE) continue;
			/*
			  C
			  C      EQUATION A.2B (FELMY AND WEARE, 1986) FOR ZETA
			  C
			  LGAMMA(J)=LGAMMA(J)+M(K)*MN(KK)*ZETA(J,K,KK)
			*/
			LGAMMA[i0] += M[i1]*M[i2]*param;
			LGAMMA[i1] += M[i0]*M[i2]*param;
			LGAMMA[i2] += M[i0]*M[i1]*param;
		}
	}

/*
C
C     CONVERT TO MACINNES CONVENTION
C
*/
      if (ICON != 0) {
	      PHIMAC=LGAMMA[IC]-GAMCLM;
/*
C
C     CORRECTED ERROR IN PHIMAC, NOVEMBER, 1989
C
*/
	      for (i = 0; i < count_cations; i++) {
		      if (IPRSNT[i] == TRUE) {
			      LGAMMA[i]=LGAMMA[i]+spec[i]->z*PHIMAC;
		      }
	      }
      }
/*
C
C     CALCULATE THE OSMOTIC COEFFICIENT
C
C     EQUATION (2A) PART 1
C
*/
      OSMOT=-(A0)*pow(I,1.5e0)/(1.0e0+B*DI);
      /*
       *  Add etheta terms to osmotic coefficient
       */
      /* OSMOT=OSMOT+M(J)*M(K)*PHIPHI() */
      /* PHIPHI=THETA(J,K)+ETHETA()+I*ETHEAP() */
      /* saving THETA for later */
      for (i = 0; i < count_cations - 1; i++) {
	      z0 = cations[i]->z;
	      if (IPRSNT[i] == FALSE) continue;
	      for (j = i+1; j < count_cations; j++) {
		      if (IPRSNT[j] == FALSE) continue;
		      z1 = cations[j]->z;
		      OSMOT += M[i]*M[j]*(ETHETA(z0, z1, I) + I*ETHETAP(z0, z1, I) ); 
	      }
      }
      for (i = 2*count_s; i < 2*count_s + count_anions - 1; i++) {
	      z0 = spec[i]->z;
	      for (j = i+1; j < 2*count_s + count_anions; j++) {
		      z1 = spec[j]->z;
		      OSMOT += M[i]*M[j]*(ETHETA(z0, z1, I) + I*ETHETAP(z0, z1, I) ); 
	      }
      }
      /*
       *    Sums of Pitzer parameters for osmotic coefficient
       */
      for (i = 0; i < count_pitz_param; i++) {
	      i0 = pitz_params[i]->ispec[0];
	      i1 = pitz_params[i]->ispec[1];
	      if (IPRSNT[i0] == FALSE || IPRSNT[i1] == FALSE) continue;
	      z0 = spec[pitz_params[i]->ispec[0]]->z;
	      z1 = spec[pitz_params[i]->ispec[1]]->z;
	      param = pitz_params[i]->p;
	      alpha = pitz_params[i]->alpha;
	      /*
		C
		C     EQUATION (2B) PART 3
		C
		OSMOT=OSMOT+M(J)*M(K)*(BMXPHI()+BIGZ*CMX())
		BMXPHI=BCX(1,J,K)+BCX(2,J,K)*DEXP(-ALPHA(2)*DSQRT(I))+BCX(3,J,K)*
		       DEXP(-ALPHA(3)*DSQRT(I))
		CMX=BCX(4,J,K)/(2.0D0*DSQRT(DABS(Z(J)*Z(K)))) 
	      */
	      if (pitz_params[i]->type == TYPE_B0) {
		      OSMOT += M[i0]*M[i1]*param;
	      } else if (pitz_params[i]->type == TYPE_B1) {
		      OSMOT += M[i0]*M[i1]*param*exp(-alpha*DI);
	      } else if (pitz_params[i]->type == TYPE_B2) {
		      OSMOT += M[i0]*M[i1]*param*exp(-alpha*DI);
	      } else if (pitz_params[i]->type == TYPE_C0) {
		      OSMOT += M[i0]*M[i1]*BIGZ*param/(2.0*sqrt(fabs(z0*z1))); 
	      } else if (pitz_params[i]->type == TYPE_THETA) {
		      /*
			C
			C     EQUATION (2B) PART 2
			C
			OSMOT=OSMOT+M(J)*M(K)*PHIPHI() 
			PHIPHI=THETA(J,K)+ETHETA()+I*ETHEAP() 
			etheta and ethetap added above in case theta is not defined for a pair
		      */
		      OSMOT += M[i0]*M[i1]*param;
	      } else if (pitz_params[i]->type == TYPE_PSI) {
		      i2 = pitz_params[i]->ispec[2];
		      if (IPRSNT[i2] == FALSE) continue;
		      z2 = spec[pitz_params[i]->ispec[2]]->z;
		      /*
			C
			C     EQUATION (2B) PART 2
			C
			OSMOT=OSMOT+M(J)*M(K)*M(KK)*PSI(J,K,KK)
		      */
		      OSMOT += M[i0]*M[i1]*M[i2]*param;
	      } else if (pitz_params[i]->type == TYPE_LAMDA) {
		      /* OSMOT=OSMOT+MN(K)*M(J)*LAM(J,K) */
		      OSMOT += M[i0]*M[i1]*param;
	      } else if (pitz_params[i]->type == TYPE_ZETA) {
		      i2 = pitz_params[i]->ispec[2];
		      if (IPRSNT[i2] == FALSE) continue;
		      /*
			C
			C      EQUATION A.2B (FELMY AND WEARE, 1986) FOR ZETA
			C
			OSMOT=OSMOT+MN(K)*M(J)*M(KK)*ZETA(J,KK,K)
		      */
		      OSMOT += M[i0]*M[i1]*M[i2]*param;
	      }
      }
      COSMOT = 1.0e0 + 2.0e0*OSMOT/OSUM;
/*
C
C     CALCULATE THE ACTIVITY OF WATER
C
*/
      AW=exp(-OSUM*COSMOT/55.50837e0);
      for (i = 0; i < 2*count_s + count_anions; i++) {
	      if (IPRSNT[i] == FALSE) continue;
	      spec[i]->lg=LGAMMA[i]*CONV;
	      fprintf(stderr, "%s: %e\n", spec[i]->name, LGAMMA[i]*CONV);
      }
      fprintf(stderr, "COSMOT: %e\n", COSMOT);
      
      *I_X = I;
      *COSMOT_X = COSMOT;
      return(OK);
}
#endif
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
#ifdef SKIP
/* ---------------------------------------------------------------------- */
int DH_AB (double TK, double *A, double *B)
/* ---------------------------------------------------------------------- */
{
/*
C
C          COMPUTE TEMPERATURE DEPENDENCE OF A AND B FOR DEBYE-HUCKEL
C
*/
	double S1, S2, S3, C1, TC;
	TC = TK -273.15;
	S1=374.11e0-TC;
	S2=pow(S1,0.33333333e0);
	S3=1.0e0+0.1342489e0*S2-3.946263e-03*S1;
	S3=S3/(3.1975e0-0.3151548e0*S2-1.203374e-03*S1+7.48908e-13*pow(S1,4.0e0));
	S3=sqrt(S3);
	if (TK < 373.15e0) {
		C1=87.74e0-TC*(TC*(1.41e-06*TC-9.398e-04)+0.4008e0);
	} else {
		C1=5321.0e0/TK+233.76e0-TK*(TK*(8.292e-07*TK-1.417e-03)+0.9297e0);
	}
	C1=sqrt(C1*TK);
	*A=1824600.0e0*S3/pow(C1,3.0e0);
	*B=50.29e0*S3/C1;
	return (OK);
}
#endif
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
		pitz_params[i] = free_check_null(pitz_params[i]);
	}
	pitz_params = free_check_null(pitz_params);
	LGAMMA = free_check_null(LGAMMA);
	IPRSNT = free_check_null(IPRSNT);
	spec = free_check_null(spec);
	M = free_check_null(M);

	return OK;
}
