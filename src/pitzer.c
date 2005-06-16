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
struct species **ions, **cations, **anions, **neutrals;
int count_cations, count_anions, count_neutrals;
int MAXCATIONS, FIRSTANION, MAXNEUTRAL;
/* routines */
int calc_pitz_param (struct pitz_param *pz_ptr, double TK, double TR);

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
int pitzer_tidy (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Make lists of species for cations, anions, neutral
 */
	int i;

	if (ions != NULL) ions = free_check_null(ions);
	ions = PHRQ_malloc((size_t) (2*count_s*sizeof(struct species *)));
	cations = ions;
	anions = &(ions[count_s]);
	MAXCATIONS = count_s;
	FIRSTANION = count_s;
	MAXNEUTRAL = count_s;
	if (neutrals != NULL) neutrals = free_check_null(neutrals);
	neutrals = PHRQ_malloc((size_t) (count_s*sizeof(struct species *)));

	count_cations = 0;
	count_anions = 0;
	count_neutrals = 0;
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
	DW(298.15);
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
	    case TYPE_PSI:
		pz_ptr->U.psi = param;
		break;
	    case TYPE_Other:
	        error_msg("Should not be TYPE_Other in function calc_pitz_param", STOP);
	        break;
	}
	return OK;
}


/*
      SUBROUTINE PITZER (MO,LG,AW,MU,ICON)
      INCLUDE 'phrqpitz.inc'
      INCLUDE 'pitz.inc'
C  dlp 11/5/98
C      DOUBLE PRECISION MU,MO(MAXS),LG(MAXS),MN(20),LGN(20),LGAMMA(40),
C     &                 M(40)
C      LOGICAL IPRSNT(40),LNEUT
      DOUBLE PRECISION MU,MO(MAXS),LG(MAXS),MN(MAXNEUTRAL),
     &     LGN(MAXNEUTRAL), LGAMMA(MAXIONS), M(MAXIONS)
      LOGICAL IPRSNT(MAXIONS),LNEUT
C
      DOUBLE PRECISION BC
C dlp 11/5/98
C      COMMON / MX1 / BC(4,20,21:40,5)
      COMMON / MX1 / BC(4,MAXCATIONS,FIRSTANION:MAXIONS,5)
C
      DOUBLE PRECISION I
      INTEGER J, K
      COMMON / MX2 / I,J,K
C
      DOUBLE PRECISION BCX, OTEMP
C dlp 11/5/98
C      COMMON / MX3 / BCX(4,20,21:40), OTEMP
      COMMON / MX3 / BCX(4,MAXCATIONS,FIRSTANION:MAXIONS), OTEMP
C
      DOUBLE PRECISION A0
      COMMON / MX4 / A0
C
      DOUBLE PRECISION Z
C dlp 11/5/98
C      COMMON / MX5 / Z(40)
      COMMON / MX5 / Z(MAXIONS)
C
      DOUBLE PRECISION LAM, PSI, ZETA
      INTEGER TRANS, IN
C dlp 11/5/98
C      COMMON / MX6 / LAM(40,20),PSI(40,40,40),ZETA(20,21:40,20),
C     &               TRANS(40),IN(20)
      COMMON / MX6 / LAM(MAXIONS,MAXNEUTRAL),
     &     PSI(MAXIONS,MAXIONS,MAXIONS),
     &     ZETA(MAXCATIONS,FIRSTANION:MAXIONS,MAXNEUTRAL),
     &     TRANS(MAXIONS),IN(MAXNEUTRAL)
C
      DOUBLE PRECISION THETA
C dlp 11/5/98
C      COMMON / MX7 / THETA(40,40)
      COMMON / MX7 / THETA(MAXIONS,MAXIONS)

      DOUBLE PRECISION B
      COMMON / MX0 / B
C
      INTEGER M1, M2, M3
      COMMON / PI1 / M1,M2,M3
C
      DOUBLE PRECISION COSMOT
      COMMON / COS / COSMOT
C
      DOUBLE PRECISION CONV, XI, XX, OSUM, BIGZ, DI, F, XXX, GAMCLM, 
     &                 CSUM, PHIMAC, OSMOT, BMXP, ETHEAP, CMX, BMX, PHI,
     &                 BMXPHI, PHIPHI, AW
      INTEGER N, IK, IC, I1, I2, I3, KK, ISPEC, ICON
      EXTERNAL PTEMP, ISPEC, BMXP, ETHEAP, CMX, BMX, PHI, BMXPHI, PHIPHI
      INTRINSIC DABS, DSQRT, DLOG, DEXP
      DATA CONV / 0.4342944819D0 /
C
C     THIS IS THE MAIN ROUTINE FOR CALCULATING GAMMAS USING HARVIE AND
C     WEARE (1980) MODEL.  VALUES OF CONCENTRATION ARE PASSED INTO
C     SUBROUTINE PITZER VIA THE VARIABLE 'M0'.  'M0' CONTAINS MOLALITY
C     OF SPECIES, NUMBERED BY PHREEQE'S CONVENTION.  ONCE INSIDE PITZER,
C     VALUES OF CONCENTRATION OF SPECIES USED IN PITZER'S MODEL ARE
C     TRANSFERED TO THE VARIABLES 'M' AND 'MN'. FOR A LIST OF SPECIES,
C     SEE SPECIES INPUT IN 'PITZER.DATA'.
C
C     THE VARIABLE 'M' CONTAINS CONCENTRATIONS OF IONIZED SPECIES:
C                    M(1:20)  => CATIONS
C                    M(21:40) => ANIONS
C     'MN' CONTAINS CONCENTRATION OF NEUTRAL SPECIES
C
C     VARIABLES 'SPECIES' AND 'NEUT' CONTAIN NAMES OF IONIZED AND
C     NEUTRAL SPECIES RESPECTIVELY, USING THE SAME NUMBERING SYSTEM AS
C     'M' AND 'MN'.
C
C     OTHER IMPORTANT VARIABLES:
C           M1 => MAX INDEX # FOR CATIONS (<21)
C           M2 => MAX INDEX # FOR ANIONS (20 < M2 < 41)
C           M3 => MAX INDEX # FOR NEUTRALS (<21)
C           IPRESENT,NEUTRAL => LOGICAL VARIABLE ARRAYS TO INDICATE
C                               WHICH SPECIES HAVE NON-ZERO
C                               CONCENTRATION.
C           I => IONIC STRENGTH
C           BCX(1,-,-) => BETA (0)  |  SECOND SUBSCRIPT IS FOR CATION
C           BCX(2,-,-) => BETA (1)  |  INDEX;
C           BCX(3,-,-) => BETA (2)  |  THIRD SUBSCRIPT IS FOR ANION
C           BCX(4,-,-) => C PHI     |  INDEX.
C           LAM => LAMBDA, IONIZED-NEUTRAL INTERACTION PARAMETER
C           ZETA => ZETA, CATION-ANION-NEUTRAL INTERACTION PARAMETER
C           PSI => PSI, 1 CATION - 2 ANION INTERACTION OR 2 CATION -
C                  1 ANION INTERACTION PARAMETER.
C           CONV => CONVERSION FACTOR FROM LOG BASE E TO LOG BASE 10
C           LGAMMA => LOG OF GAMMAS
C           PHIMAC => CONVERSION FACTOR FOR MACINNES CONVENTION
C           COSMOTIC => OSMOTIC COEFFICIENT
C           AW => ACTIVITY OF WATER
C           GAMCLM => LOG GAMMA OF KCL SYSTEM AT SAME IONIC STRENGTH
C                     AND TEMPERATURE.  USED BY 'PHIMAC'
C           Z => CHANGE OF SPECIES
C
C     EQUATION NUMBERS REFER TO HARVIE AND WEARE (1980)
C
C     INITIALIZE
C
      XI=0.0D0
      XX=0.0D0
      OSUM=0.0D0
      LNEUT=.FALSE.
C
C     TRANSFER DATA FROM MO TO M
C
      DO 20 N=1,M2
      M(N)=0.0D0
C dlp 11/5/98
C      IF (N.GT.M1.AND.N.LT.21) GO TO 20
      IF (N.GT.M1.AND.N.LT.FIRSTANION) GO TO 20
      M(N)=MO(TRANS(N))
   20 CONTINUE
      DO 30 N=1,M3
      MN(N)=MO(IN(N))
      IF (MN(N).GT.1.0D-40) LNEUT=.TRUE.
   30 CONTINUE
C
      DO 40 N=1,M2
      IF (M(N).GT.1.0D-40) GO TO 50
      IPRSNT(N)=.FALSE.
      GO TO 40
   50 IPRSNT(N)=.TRUE.
   40 CONTINUE
C
C     COMPUTE PITZER COEFFICIENTS' TEMPERATURE DEPENDENCE
C
      CALL PTEMP
C
      DO 10 N=1,M2
      IF (.NOT.IPRSNT(N)) GO TO 10
      XX=XX+M(N)*DABS(Z(N))
      XI=XI+M(N)*Z(N)*Z(N)
      OSUM=OSUM+M(N)
   10 CONTINUE
      DO 15 N=1,M3
      OSUM=OSUM+MN(N)
   15 CONTINUE
      I=XI/2.0D0
C
C     EQUATION (8)
C
      BIGZ=XX
      DI=DSQRT(I)
C
C     CALCULATE F & GAMCLM
C
      F=-A0*(DI/(1.0D0+B*DI)+2.0D0*DLOG(1.0D0+B*DI)/B)
      XXX=2.0D0*DI
      XXX=(1.0D0-(1.0D0+XXX-XXX*XXX*0.5D0)*DEXP(-XXX))/(XXX*XXX)
      IK=ISPEC('K+      ')
      IC=ISPEC('CL-     ')
      GAMCLM=F+I*2.0D0*(BCX(1,IK,IC)+BCX(2,IK,IC)*XXX)+1.5D0
     1*BCX(4,IK,IC)*I*I
      DO 75 J=1,M1
      IF (.NOT.IPRSNT(J)) GO TO 75
C dlp 11/5/98
C      DO 70 K=21,M2
      DO 70 K=FIRSTANION,M2
      IF (.NOT.IPRSNT(K)) GO TO 70
C
C     EQUATION (3) PART 1
C
      F=F+M(J)*M(K)*BMXP()
   70 CONTINUE
   75 CONTINUE
      DO 85 J=1,M1-1
      IF (.NOT.IPRSNT(J)) GO TO 85
      DO 80 K=J+1,M1
      IF (.NOT.IPRSNT(K)) GO TO 80
C
C     EQUATION (3) PART 2
C
      F=F+M(J)*M(K)*ETHEAP()
   80 CONTINUE
   85 CONTINUE
C dlp 11/5/98
C      DO 90 J=21,M2-1
      DO 90 J=FIRSTANION,M2-1
      IF (.NOT.IPRSNT(J)) GO TO 90
      DO 110 K=J+1,M2
      IF (.NOT.IPRSNT(K)) GO TO 110
C
C     EQUATION (3) PART 3
C
      F=F+M(J)*M(K)*ETHEAP()
  110 CONTINUE
   90 CONTINUE
C
      CSUM=0.0D0
      DO 125 J=1,M1
      IF (.NOT.IPRSNT(J)) GO TO 125
C dlp 11/5/98
C      DO 120 K=21,M2
      DO 120 K=FIRSTANION,M2
      IF (.NOT.IPRSNT(K)) GO TO 120
C
C     EQUATION (2B) PART 4
C
      CSUM=CSUM+M(J)*M(K)*CMX()
  120 CONTINUE
  125 CONTINUE
C
C     CALCULATE LGAMMA FOR CATIONS
C
      DO 130 J=1,M1
      LGAMMA(J)=Z(J)*Z(J)*F+DABS(Z(J))*CSUM
C dlp 11/5/98
C      DO 140 K=21,M2
      DO 140 K=FIRSTANION,M2
      IF (.NOT.IPRSNT(K)) GO TO 140
C
C     EQUATION (2B) PART 3
C
      LGAMMA(J)=LGAMMA(J)+M(K)*(2.0D0*BMX()+BIGZ*CMX())
  140 CONTINUE
      DO 150 K=1,M1
      IF (.NOT.IPRSNT(K)) GO TO 150
      LGAMMA(J)=LGAMMA(J)+2.0D0*M(K)*PHI()
C dlp 11/5/98
C      DO 160 KK=21,M2
      DO 160 KK=FIRSTANION,M2
      IF (.NOT.IPRSNT(KK)) GO TO 160
C
C     EQUATION (2B) PART 2
C
      LGAMMA(J)=LGAMMA(J)+M(KK)*M(K)*PSI(J,K,KK)
  160 CONTINUE
  150 CONTINUE
C dlp 11/5/98
C      DO 170 K=21,M2-1
      DO 170 K=FIRSTANION,M2-1
      IF (.NOT.IPRSNT(K)) GO TO 170
      DO 180 KK=K+1,M2
      IF (.NOT.IPRSNT(KK)) GO TO 180
C
C     EQUATION (2B) PART 3
C
      LGAMMA(J)=LGAMMA(J)+M(K)*M(KK)*PSI(K,KK,J)
  180 CONTINUE
  170 CONTINUE
      IF (.NOT.LNEUT) GO TO 130
      DO 190 K=1,M3
      LGAMMA(J)=LGAMMA(J)+2.0D0*MN(K)*LAM(J,K)
  190 CONTINUE
C
C      EQUATION A.2B (FELMY AND WEARE, 1986) FOR ZETA
C
C dlp 11/5/98
C      DO 192 K=21,M2
      DO 192 K=FIRSTANION,M2
      IF (.NOT.IPRSNT(K)) GO TO 192
      DO 194 KK=1,M3
      LGAMMA(J)=LGAMMA(J)+M(K)*MN(KK)*ZETA(J,K,KK)
  194 CONTINUE
  192 CONTINUE
  130 CONTINUE
      DO 1200 I1=1,M1
C dlp 11/5/98
C      DO 1300 I2=21,M2
      DO 1300 I2=FIRSTANION,M2
      DO 1400 I3=1,M3
 1400 CONTINUE
 1300 CONTINUE
 1200 CONTINUE
C
C     CALCULATE LGAMMA OF ANIONS
C
C dlp 11/5/98
C      DO 230 K=21,M2
      DO 230 K=FIRSTANION,M2
C
C     EQUATION (2C) PART 1
C
      LGAMMA(K)=Z(K)*Z(K)*F+DABS(Z(K))*CSUM
      DO 240 J=1,M1
      IF (.NOT.IPRSNT(J)) GO TO 240
C
C     EQUATION (2C) PART 2
C
      LGAMMA(K)=LGAMMA(K)+M(J)*(2.0D0*BMX()+BIGZ*CMX())
  240 CONTINUE
C dlp 11/5/98
C      DO 250 J=21,M2
      DO 250 J=FIRSTANION,M2
      IF (.NOT.IPRSNT(J)) GO TO 250
      LGAMMA(K)=LGAMMA(K)+2.0D0*M(J)*PHI()
      DO 260 KK=1,M1
      IF (.NOT.IPRSNT(KK)) GO TO 260
C
C     EQUATION (2C) PART 3
C
      LGAMMA(K)=LGAMMA(K)+M(KK)*M(J)*PSI(K,J,KK)
  260 CONTINUE
  250 CONTINUE
      DO 270 J=1,M1-1
      IF (.NOT.IPRSNT(J)) GO TO 270
      DO 280 KK=J+1,M1
      IF (.NOT.IPRSNT(KK)) GO TO 280
C
C     EQUATION (2C) PART 4
C
      LGAMMA(K)=LGAMMA(K)+M(J)*M(KK)*PSI(J,KK,K)
  280 CONTINUE
  270 CONTINUE
      IF (.NOT.LNEUT) GO TO 230
      DO 290 J=1,M3
      LGAMMA(K)=LGAMMA(K)+2.0D0*MN(J)*LAM(K,J)
  290 CONTINUE
C
C      EQUATION (A.2C) (FELMY AND WEARE, 1986) FOR ZETA
C
      DO 292 J=1,M1
      IF (.NOT.IPRSNT(J)) GO TO 292
      DO 294 KK=1,M3
      LGAMMA(K)=LGAMMA(K)+MN(KK)*M(J)*ZETA(J,K,KK)
  294 CONTINUE
  292 CONTINUE
  230 CONTINUE
C
C     CONVERT TO MACINNES CONVENTION
C
      IF (ICON.EQ.0) GO TO 300
      PHIMAC=LGAMMA(IC)-GAMCLM
C
C     CORRECTED ERROR IN PHIMAC, NOVEMBER, 1989
C
      DO 220 K=1,M2
      IF (.NOT.IPRSNT(K)) GO TO 220
      LGAMMA(K)=LGAMMA(K)+Z(K)*PHIMAC
  220 CONTINUE
C
  300 IF (.NOT.LNEUT) GO TO 860
C
C     CALCULATE THE GAMMA OF NEUTRAL IONS
C
      DO 800 K=1,M3
      LGN(K)=0.0D0
      DO 870 J=1,M2
      IF (.NOT.IPRSNT(J)) GO TO 870
      LGN(K)=LGN(K)+2.0D0*M(J)*LAM(J,K)
  870 CONTINUE
C
C     EQUATION (A.2D) (FELMY AND WEARE, 1986) FOR ZETA
C
      DO 350 J=1,M1
      IF (.NOT.IPRSNT(J)) GO TO 350
C dlp 11/5/98
C      DO 360 KK=21,M2
      DO 360 KK=FIRSTANION,M2
      IF (.NOT.IPRSNT(KK)) GO TO 360
      LGN(K)=LGN(K)+M(J)*M(KK)*ZETA(J,KK,K)
  360 CONTINUE
  350 CONTINUE
  800 CONTINUE
C
C     CALCULATE THE OSMOTIC COEFFICIENT
C
C     EQUATION (2A) PART 1
C
  860 OSMOT=-(A0)*I**1.5D0/(1.0D0+B*DI)
      DO 420 J=1,M1
      IF (.NOT.IPRSNT(J)) GO TO 420
C dlp 11/5/98
C      DO 430 K=21,M2
      DO 430 K=FIRSTANION,M2
      IF (.NOT.IPRSNT(K)) GO TO 430
C
C     EQUATION (2A) PART 2
C
      OSMOT=OSMOT+M(J)*M(K)*(BMXPHI()+BIGZ*CMX())
  430 CONTINUE
  420 CONTINUE
      DO 440 J=1,M1-1
      IF (.NOT.IPRSNT(J)) GO TO 440
      DO 450 K=J+1,M1
      IF (.NOT.IPRSNT(K)) GO TO 450
      OSMOT=OSMOT+M(J)*M(K)*PHIPHI()
C dlp 11/5/98
C      DO 460 KK=21,M2
      DO 460 KK=FIRSTANION,M2
      IF (.NOT.IPRSNT(KK)) GO TO 460
C
C     EQUATION (2A) PART 3
C
      OSMOT=OSMOT+M(J)*M(K)*M(KK)*PSI(J,K,KK)
  460 CONTINUE
  450 CONTINUE
  440 CONTINUE
C dlp 11/5/98
C      DO 470 J=21,M2-1
      DO 470 J=FIRSTANION,M2-1
      IF (.NOT.IPRSNT(J)) GO TO 470
      DO 480 K=J+1,M2
      IF (.NOT.IPRSNT(K)) GO TO 480
      OSMOT=OSMOT+M(J)*M(K)*PHIPHI()
      DO 490 KK=1,M1
      IF (.NOT.IPRSNT(KK)) GO TO 490
C
C     EQUATION (2A) PART 4
C
      OSMOT=OSMOT+M(J)*M(K)*M(KK)*PSI(J,K,KK)
  490 CONTINUE
  480 CONTINUE
  470 CONTINUE
      IF (.NOT.LNEUT) GO TO 850
      DO 810 K=1,M3
      DO 820 J=1,M2
      IF (.NOT.IPRSNT(J)) GO TO 820
C
C     EQUATION (A.3A) PART 5  HARVIE, MOLLER, WEARE (1984)
C
      OSMOT=OSMOT+MN(K)*M(J)*LAM(J,K)
  820 CONTINUE
C
C      EQUATION (A.2A) (FELMY AND WEARE, 1986) FOR ZETA
C
      DO 950 J=1,M1
      IF (.NOT.IPRSNT(J)) GO TO 950
C dlp 11/5/98
C      DO 960 KK=21,M2
      DO 960 KK=FIRSTANION,M2
      IF (.NOT.IPRSNT(KK)) GO TO 960
      OSMOT=OSMOT+MN(K)*M(J)*M(KK)*ZETA(J,KK,K)
  960 CONTINUE
  950 CONTINUE
  810 CONTINUE
  850 COSMOT=1.0D0+2.0D0*OSMOT/OSUM
C
C     CALCULATE THE ACTIVITY OF WATER
C
      AW=DEXP(-OSUM*COSMOT/55.50837D0)
C
C     SET APPROPRIATE VALUES FOR RETURN
C
      MU=I
      DO 900 N=1,M2
      IF (.NOT.IPRSNT(N)) GO TO 900
      LG(TRANS(N))=LGAMMA(N)*CONV
  900 CONTINUE
      IF (.NOT.LNEUT) RETURN
      DO 910 N=1,M3
      LG(IN(N))=LGN(N)*CONV
  910 CONTINUE
      RETURN
      END
C
C
C
      DOUBLE PRECISION FUNCTION BMXPHI()
      INCLUDE 'pitz.inc'
C
C     THESE FUNCTIONS CALCULATE THE BM'S
C
      DOUBLE PRECISION I
      INTEGER J, K
      COMMON / MX2 / I,J,K
C
      DOUBLE PRECISION BCX, OTEMP
C dlp 11/5/98
C      COMMON / MX3 / BCX(4,20,21:40), OTEMP
      COMMON / MX3 / BCX(4,MAXCATIONS,FIRSTANION:MAXIONS), OTEMP
C
C     ALPHA DEFINED IN BLOCK DATA
C
      DOUBLE PRECISION ALPHA
      COMMON / MX9 / ALPHA(3)
C
      DOUBLE PRECISION BMX, BMXP, G, GP, Y
      LOGICAL TWOTWO
      EXTERNAL TWOTWO
      INTRINSIC DEXP, DSQRT
C
C     STATEMENT FUNCTIONS
C
C     EQUATION (5)
C
      G(Y)=2.0D0*(1.0D0-(1.0D0+Y)*DEXP(-Y))/Y**2.0D0
C
C     EQUATION (6)
C
      GP(Y)=-2.0D0*(1.0D0-(1.0D0+Y+Y**2.0D0/2.0D0)*DEXP(-Y))/
     1 Y**2.0D0
C
      IF (TWOTWO()) GO TO 50
C
C     EQUATION (4A)
C
      BMXPHI=BCX(1,J,K)+BCX(2,J,K)*DEXP(-ALPHA(1)*DSQRT(I))
      RETURN
C
C     EQUATION (7A)
C
   50 BMXPHI=BCX(1,J,K)+BCX(2,J,K)*DEXP(-ALPHA(2)*DSQRT(I))+BCX(3,J,K)*
     1 DEXP(-ALPHA(3)*DSQRT(I))
      RETURN
C     *********
      ENTRY BMX
C     *********
      IF (TWOTWO()) GO TO 60
C
C     EQUATION (4B)
C
      BMX=BCX(1,J,K)+BCX(2,J,K)*G(ALPHA(1)*DSQRT(I))
      RETURN
C
C     EQUATION (7B)
C
   60 BMX=BCX(1,J,K)+BCX(2,J,K)*G(ALPHA(2)*DSQRT(I))+BCX(3,J,K)*
     1 G(ALPHA(3)*DSQRT(I))
      RETURN
C     **********
      ENTRY BMXP
C     **********
      IF (TWOTWO()) GO TO 70
C
C     EQUATION (4C)
C
      BMXP=BCX(2,J,K)*GP(ALPHA(1)*DSQRT(I))/I
      RETURN
   70 CONTINUE
C
C     EQUATION (7C)
C
      BMXP=(BCX(2,J,K)*GP(ALPHA(2)*DSQRT(I))+BCX(3,J,K)*GP(ALPHA(3)*
     1 DSQRT(I)))/I
      RETURN
      END
C
C
C
      LOGICAL FUNCTION TWOTWO()
      INCLUDE 'pitz.inc'
C
C     SUBROUTINE TO DETERMINE TWO-TWO ELECTROLYTES
C     TWOTWO = .TRUE.  IF BOTH IONS ARE DOUBLY CHARGED
C
      DOUBLE PRECISION I
      INTEGER J, K
      COMMON / MX2 / I,J,K
C
      DOUBLE PRECISION Z
C dlp 11/5/98
C      COMMON / MX5 / Z(40)
      COMMON / MX5 / Z(MAXIONS)
C
      DOUBLE PRECISION ZDIF
      INTRINSIC DABS
      TWOTWO=.FALSE.
      ZDIF=DABS(DABS(Z(J)*Z(K))-4.0D0)
      IF (ZDIF.LT.1.0D0) TWOTWO=.TRUE.
      RETURN
      END
C
C
C
      DOUBLE PRECISION FUNCTION CMX()
      INCLUDE 'pitz.inc'
C
C     FUNCTION TO CALCULATE C(MX)
C
      DOUBLE PRECISION I
      INTEGER J, K
      COMMON / MX2 / I,J,K
C
      DOUBLE PRECISION BCX, OTEMP
C dlp 11/5/98
C      COMMON / MX3 / BCX(4,20,21:40), OTEMP
      COMMON / MX3 / BCX(4,MAXCATIONS,FIRSTANION:MAXIONS), OTEMP
C
      DOUBLE PRECISION Z
C dlp 11/5/98
C      COMMON / MX5 / Z(40)
      COMMON / MX5 / Z(MAXIONS)
C
      INTRINSIC DSQRT, DABS
C
C     EQUATION (9)
C
      CMX=BCX(4,J,K)/(2.0D0*DSQRT(DABS(Z(J)*Z(K))))
      RETURN
      END
C
C
C
      DOUBLE PRECISION FUNCTION ETHETA()
      INCLUDE 'pitz.inc'
C
C     FUNCTIONS TO CALCULATE ETHETA AND ETHEAP
C
      DOUBLE PRECISION I
      INTEGER J, K
      COMMON / MX2 / I,J,K
C
      DOUBLE PRECISION A0
      COMMON / MX4 / A0
C
      DOUBLE PRECISION Z
C dlp 11/5/98
C      COMMON / MX5 / Z(40)
      COMMON / MX5 / Z(MAXIONS)
C
      DOUBLE PRECISION JAY,JPRIME, XCON, ZZ, XJK, XJJ, XKK, ETHEAP
      INTEGER IFLAG
      EXTERNAL JAY, JPRIME
      INTRINSIC DSQRT, DABS
C
      IFLAG=1
      ETHETA=0.0D0
      GO TO 10
C     *************
      ENTRY ETHEAP
C     *************
      IFLAG=2
      ETHEAP=0.0D0
   10 XCON=6.0D0*A0*DSQRT(I)
      IF (DABS(Z(J)-Z(K)).LE.1.0D-40) RETURN
      ZZ=Z(J)*Z(K)
C
C     NEXT 3 ARE EQUATION (A1)
C
      XJK=XCON*ZZ
      XJJ=XCON*Z(J)*Z(J)
      XKK=XCON*Z(K)*Z(K)
C
C     EQUATION (A2)
C
      ETHETA=ZZ*(JAY(XJK)-JAY(XJJ)/2.0D0-JAY(XKK)/2.0D0)/(4.0D0*I)
      IF (IFLAG.EQ.1) RETURN
C
C     EQUATION (A3)
C
      ETHEAP=ZZ*(JPRIME(XJK)-JPRIME(XJJ)/2.0D0-JPRIME(XKK)/
     1 2.0D0)/(8.0D0*I**2.0D0) - ETHETA/I
      RETURN
      END
C
C     FUNCTION TO CALCULATE JAY AND JPRIME
C
C     J0 AND J1, USED IN CALCULATION OF ETHETA AND ETHEAP
C
      DOUBLE PRECISION FUNCTION JAY(X)
C
      DOUBLE PRECISION AK, BK, DK
      COMMON / MX8 / AK(0:20,2),BK(0:22),DK(0:22)
C
      DOUBLE PRECISION JPRIME, DZ, X, Y
      EXTERNAL BDK
C
      CALL BDK (X)
      JAY=X/4.0D0-1.0D0+0.5D0*(BK(0)-BK(2))
      RETURN
C     **************
      ENTRY JPRIME(Y)
C     **************
      CALL BDK (Y)
      IF (Y.GT.1.0D0) GO TO 10
      DZ=0.8D0*Y**(-0.8D0)
      GO TO 20
   10 DZ=-4.0D0*Y**(-1.1D0)/9.0D0
   20 JPRIME=Y*(.25D0+DZ*(DK(0)-DK(2))/2.0D0)
      RETURN
      END
C
C
C
      SUBROUTINE BDK (X)
C
C     NUMERICAL APPROXIMATION TO THE INTEGRALS IN THE EXPRESSIONS FOR J0
C     AND J1.  CHEBYSHEV APPROXIMATION IS USED.  THE CONSTANTS 'AK' ARE
C     DEFINED IN BLOCK COMMON.
C
      DOUBLE PRECISION AK, BK, DK
      COMMON / MX8 / AK(0:20,2),BK(0:22),DK(0:22)
C
      DOUBLE PRECISION Z,X
      INTEGER II, M
C
      II=1
      IF (X.LE.1.0D0) GO TO 10
      II=2
      Z=40.0D0*X**(-1.0D-1)/9.0D0-22.0D0/9.0D0
      GO TO 20
   10 Z=4.0D0*X**(0.2D0)-2.0D0
   20 DO 30 M=20,0,-1
      BK(M)=Z*BK(M+1)-BK(M+2)+AK(M,II)
      DK(M)=BK(M+1)+Z*DK(M+1)-DK(M+2)
   30 CONTINUE
      RETURN
      END
C
C     FUNCTION TO CALCULATE PHI   (PHI' IS EQUAL TO ETHEAP)
C
      DOUBLE PRECISION FUNCTION PHI()
      INCLUDE 'pitz.inc'
C
      DOUBLE PRECISION I
      INTEGER J, K
      COMMON / MX2 / I,J,K
C
      DOUBLE PRECISION THETA
C dlp 11/5/98
C      COMMON / MX7 / THETA(40,40)
      COMMON / MX7 / THETA(MAXIONS,MAXIONS)
C
      DOUBLE PRECISION ETHETA, ETHEAP, PHIPHI
      EXTERNAL ETHETA, ETHEAP
C
C     EQUATION (10B)
C
      PHI=THETA(J,K)+ETHETA()
      RETURN
C     ************
      ENTRY PHIPHI
C     ************
C
C     EQUATION (10A)
C
      PHIPHI=THETA(J,K)+ETHETA()+I*ETHEAP()
      RETURN
      END
C
C     SUBROUTINE TO CALUCLATE TEMPERATURE DEPENDENCE OF PITZER PARAMETER
C
      SUBROUTINE PTEMP
      INCLUDE 'phrqpitz.inc'
      INCLUDE 'phrqpitz.cb'
      INCLUDE 'pitz.inc'
C
      DOUBLE PRECISION BC
C dlp 11/5/98
C      COMMON / MX1 / BC(4,20,21:40,5)
      COMMON / MX1 / BC(4,MAXCATIONS,FIRSTANION:MAXIONS,5)
C
      DOUBLE PRECISION BCX, OTEMP
C dlp 11/5/98
C      COMMON / MX3 / BCX(4,20,21:40), OTEMP
      COMMON / MX3 / BCX(4,MAXCATIONS,FIRSTANION:MAXIONS), OTEMP
C
      DOUBLE PRECISION A0
      COMMON / MX4 / A0
C
      DOUBLE PRECISION VP, DW0
      COMMON / MX10 / VP,DW0
C
      INTEGER M1, M2, M3
      COMMON / PI1 / M1,M2,M3
C
      DOUBLE PRECISION DC, TR, DC0
      INTEGER L, K, N
      EXTERNAL DW, DC
      INTRINSIC DABS, DLOG
      DATA TR /298.15/
C
      IF (DABS(TK-OTEMP).LT.0.01D0) RETURN
      OTEMP=TK
C     Set DW0
      CALL DW(TK)
      IF (DABS(TK-TR).LT.0.01D0) GO TO 10
      DC0=DC(TK)
      A0=1.400684D6*(DW0/(DC0*TK)**3.0D0)**0.5D0
      DO 30 L=1,4
      DO 40 K=1,M1
C dlp 11/5/98
C      DO 50 N=21,M2
      DO 50 N=FIRSTANION,M2
      BCX(L,K,N)=BC(L,K,N,1)+BC(L,K,N,2)*(1.0D0/TK-1.0D0/TR)+BC(L,K,N,3)
     1*DLOG(TK/TR)+BC(L,K,N,4)*(TK-TR)+BC(L,K,N,5)*(TK*TK-TR*TR)
   50 CONTINUE
   40 CONTINUE
   30 CONTINUE
      RETURN
   10 DO 20 L=1,4
      DO 60 K=1,M1
C dlp 11/5/98
C      DO 70 N=21,M2
      DO 70 N=FIRSTANION,M2
      BCX(L,K,N)=BC(L,K,N,1)
   70 CONTINUE
   60 CONTINUE
   20 CONTINUE
      A0=0.392D0
      RETURN
      END
C
C     FUNCTION TO TRANSFORM SPECIES NAME INTO SPECIES NUMBER
C
      INTEGER FUNCTION ISPEC(SPEC)
      INCLUDE 'phrqpitz.inc'
      INCLUDE 'pitz.inc'
      CHARACTER*8 SPEC
C
      CHARACTER*8 SPECS, NEUTRL
C dlp 11/5/98
C      COMMON / PI1C / SPECS(40),NEUTRL(20)
      COMMON / PI1C / SPECS(MAXIONS),NEUTRL(MAXNEUTRAL)
C
      INTEGER M1, M2, M3
      COMMON / PI1 / M1,M2,M3
C
      INTEGER I
C
      DO 10 I=1,M2
C dlp 11/5/98
C      IF (I.GT.M1.AND.I.LT.21) GO TO 10
      IF (I.GT.M1.AND.I.LT.FIRSTANION) GO TO 10
      IF (SPEC.EQ.SPECS(I)) GO TO 20
   10 CONTINUE
      DO 40 I=1,M3
      IF (SPEC.EQ.NEUTRL(I)) GO TO 20
   40 CONTINUE
      WRITE(IW,30) SPEC
   30 FORMAT (' ******* ERROR IN DATA BASE ******   ',A8,' IS NOT A '
     1,'VALID ION.')
      ENDFILE (UNIT=IW)
      CLOSE (UNIT=IW)
      CLOSE (UNIT=IR)
      STOP
   20 ISPEC=I
      RETURN
      END
*/
