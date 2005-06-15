#define EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "output.h"
#include "phrqproto.h"
#define PITZER
#include "pitzer.h"

static char const svnid[] = "$Id: pitzer.c 248 2005-04-14 17:10:53Z dlpark $";

int BB (double T);
double PS (double T);
double VLEST (double T);
int DFIND (double DOUT, double P, double D, double T);


/* ---------------------------------------------------------------------- */
int pitzer_init (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Initialization for pitzer
 */
	max_pitz_param = 100;
	count_pitz_param = 0;
	space ((void *) &pitz_params, INIT, &max_pitz_param, sizeof(struct pitz_param *));
	return OK;
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
	return(return_value);
}
/* ---------------------------------------------------------------------- */
int DW (double T)
/* ---------------------------------------------------------------------- */
/*
C
C      SUBROUTINE TO CALCULATE THE DENSITY OF WATER AS A FUNCTION OF
C      TEMPERATURE.  T IS IN KELVIN, P IS IN PASCALS, DW0 IS IN G/CM^3
C
C     FROM L. HAAR, J. S. GALLAGHER, AND G. S. KELL, (1984)
C
*/
{
	double FP = 9.869232667e0, P, DGSS, D;
	
	BB (T);
	P=1.0e0/FP;
	if (T > 373.149e0) P=PS(T);
	DGSS=P/T/.4e0;
	if (T < TZ) {
		DGSS=1.0e0/(VLEST(T));
	}
	DFIND (D,P,DGSS,T);
	DW0=D;
	VP=P*FP;
	return OK;
}

/* ---------------------------------------------------------------------- */
int BB (double T)
/* ---------------------------------------------------------------------- */
/*

C
C     THIS SUBROUTINE CALCULATES THE B'S NEEDED FOR FUNCTION DW.
C     THE B'S CALCULATED HERE ARE IN CM3/G.
C
C     FROM L. HAAR, J. S. GALLAGHER, AND G. S. KELL, (1984)
C
*/
{
	double V[11];
	int I;

	V[1]=1.0;
	for (I = 2; I <= 10; I++) {
		V[I]=V[I-1]*TZ/T;
	}
	B1=P[1]+P[2]*log(1.e0/V[2]);
	B2=Q[1];
	B1T=P[2]*V[2]/TZ;
	B2T=0.e0;
	B1TT=0.e0;
	B2TT=0.e0;
	for (I = 3; I <= 10; I++) {
		B1=B1+P[I]*V[I-1];
		B2=B2+Q[I]*V[I-1];
		B1T=B1T-(I-2)*P[I]*V[I-1]/T;
		B2T=B2T-(I-2)*Q[I]*V[I-1]/T;
		B1TT=B1TT+P[I]*(I-2)*(I-2)*V[I-1]/T/T;
		B2TT=B2TT+Q[I]*(I-2)*(I-2)*V[I-1]/T/T;
	}
	B1TT=B1TT-B1T/T;
	B2TT=B2TT-B2T/T;
	return OK;
}
/* ---------------------------------------------------------------------- */
double PS (double T)
/* ---------------------------------------------------------------------- */
/*
C
C     THIS FUNCTION CALCULATES AN APPROXIMATION TO THE VAPOR PRESSURE, P
C     AS A FUNCTION OF THE INPUT TEMPERATURE. THE VAPOR PRESSURE
C     CALCULATED AGREES WITH THE VAPOR PRESSURE PREDICTED BY THE SURFACE
C     TO WITHIN .02% TO WITHIN A DEGREE OR SO OF THE CRITICAL TEMPERATUR
C     AND CAN SERVE AS AN INITIAL GUESS FOR FURTHER REFINEMENT BY
C     IMPOSING THE CONDITION THAT GL=GV.
C
C     FROM L. HAAR, J. S. GALLAGHER, AND G. S. KELL, (1984)
C
*/
{
	double A[9]={0, -7.8889166e0,2.5514255e0,-6.716169e0,
		     33.239495e0,-105.38479e0,174.35319e0,-148.39348e0,
		     48.631602e0};
	double PL, V, W, B, Z, Q;
	int I;
	if(T <= 314.e0) {
		PL=6.3573118e0-8858.843e0/T+607.56335e0*pow(T,-.6e0);
		return (.1e0*exp(PL));
	}
	V=T/647.25e0;
	W=fabs(1.e0-V);
	B=0.e0;
	for (I = 1; I <= 8; I++) {
		Z=I;
		B=B+A[I]*pow(W,((Z+1.e0)/2.e0));
	}
	Q=B/V;
	return (22.093e0*exp(Q));
}
/* ---------------------------------------------------------------------- */
double VLEST (double T)
/* ---------------------------------------------------------------------- */
/*
C
C     FROM L. HAAR, J. S. GALLAGHER, AND G. S. KELL, (1984)
C
*/
{
	double A=-1.59259e1,B=6.57886e-2,C=-1.12666e-4,D=7.33191e-8,
		E=1.60229e3,F=2.88572e0,G=650.0e0;

	return (A+B*T+C*T*T+D*T*T*T+E/T+F/(G-T));
}
/* ---------------------------------------------------------------------- */
int DFIND (double DOUT, double P, double D, double T)
/* ---------------------------------------------------------------------- */
/*
C
C     ROUTINE TO FIND DENSITY CORRESPONDING TO INPUT PRESSURE P(MPA), AN
C     TEMPERATURE T(K), USING INITIAL GUESS DENSITY D(G/CM3). THE OUTPUT
C     DENSITY IS IN G/CM3.
C
C     FROM L. HAAR, J. S. GALLAGHER, AND G. S. KELL, (1984)
C
*/
{
	int L;
	double DD, RT, PP_dfind, DPD, DPDX, DP, X; 
	/*	double DD, RT, PP, DPD, DPDX, DP, X; */

	DD=D;
	RT=GASCON*T;
	if(DD <= 0.e0) DD=1.e-8;
	if(DD > 1.9e0) DD=1.9e0;
	L=0;
	for (L=1; L <= 30; L++) {
		if(DD <= 0.e0) DD=1.e-8;
		if(DD > 1.9e0) DD=1.9e0;
		QQ(T,DD);
		PP_dfind = RT*DD*BASE(DD)+Q0;
		DPD=RT*(Z+Y*DZ)+Q5;
		/*
C
C  THE FOLLOWING 3 LINES CHECK FOR NEGATIVE DP/DRHO, AND IF SO ASSUME
C    GUESS TO BE IN 2-PHASE REGION, AND CORRECT GUESS ACCORDINGLY.
C
		*/
		if(DPD <= 0.e0) {
			if(D >= .2967e0) DD=DD*1.02e0;
			if(D < .2967e0) DD=DD*.98e0;
			if(L <= 10) continue;
		} else {
/* 13 */
			DPDX=DPD*1.1e0;
			if(DPDX < 0.1e0) DPDX=0.1e0;
			DP=fabs(1.e0-PP_dfind/P);
			if(DP < 1.e-8) break;
			if(D > .3e0 && DP < 1.e-7) break;
			if(D > .7e0 && DP < 1.e-6) break;
			X=(P-PP_dfind)/DPDX;
			if(fabs(X) > .1e0) X=X*.1e0/fabs(X);
			DD=DD+X;
			if(DD < 0.e0) DD=1.e-8;
		}
	}
	if (L > 30) error_msg("In subroutine DFIND", STOP);
/*   20 CONTINUE */
	DOUT=DD;
	return OK;
}
/*
C
C     ROUTINE TO FIND DENSITY CORRESPONDING TO INPUT PRESSURE P(MPA), AN
C     TEMPERATURE T(K), USING INITIAL GUESS DENSITY D(G/CM3). THE OUTPUT
C     DENSITY IS IN G/CM3.
C
C     FROM L. HAAR, J. S. GALLAGHER, AND G. S. KELL, (1984)
C
      SUBROUTINE DFIND(DOUT,P,D,T)
C
      DOUBLE PRECISION Q0, Q5
      COMMON /QQQQ/ Q0,Q5
C
      DOUBLE PRECISION GASCON, TZ, AA, Z, DZ, Y
      COMMON /ACONST/ GASCON,TZ,AA,Z,DZ,Y
C
      DOUBLE PRECISION BASE, DD, RT, PP, DPD, DPDX, DP, X, DOUT, P, D, T
      INTEGER L
      EXTERNAL QQ, BASE
      INTRINSIC DABS
C
      DD=D
      RT=GASCON*T
      IF(DD.LE.0.D0) DD=1.D-8
      IF(DD.GT.1.9D0) DD=1.9D0
      L=0
    9 L=L+1
      IF(DD.LE.0.D0) DD=1.D-8
      IF(DD.GT.1.9D0) DD=1.9D0
      CALL QQ(T,DD)
      PP = RT*DD*BASE(DD)+Q0
      DPD=RT*(Z+Y*DZ)+Q5
C
C  THE FOLLOWING 3 LINES CHECK FOR NEGATIVE DP/DRHO, AND IF SO ASSUME
C    GUESS TO BE IN 2-PHASE REGION, AND CORRECT GUESS ACCORDINGLY.
C
      IF(DPD.GT.0.D0) GO TO 13
      IF(D.GE..2967D0) DD=DD*1.02D0
      IF(D.LT..2967D0) DD=DD*.98D0
      IF(L.LE.10) GO TO 9
  13  DPDX=DPD*1.1D0
      IF(DPDX.LT.0.1D0) DPDX=0.1D0
      DP=DABS(1.D0-PP/P)
      IF(DP.LT.1.D-8) GO TO 20
      IF(D.GT..3D0 .AND. DP.LT.1.D-7) GO TO 20
      IF(D.GT..7D0 .AND. DP.LT.1.D-6) GO TO 20
      X=(P-PP)/DPDX
      IF(DABS(X).GT..1D0) X=X*.1D0/DABS(X)
      DD=DD+X
      IF(DD.LE.0.D0) DD=1.D-8
      IF(L.LE.30) GO TO 9
      STOP 1
   20 CONTINUE
      DOUT=DD
      RETURN
      END
*/
