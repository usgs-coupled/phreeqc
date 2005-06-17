static double BK[23], DK[23];
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
#ifdef SKIP
/* ---------------------------------------------------------------------- */
double PHI (int J, int K)
/* ---------------------------------------------------------------------- */
/*
C
C     FUNCTION TO CALCULATE PHI   (PHI' IS EQUAL TO ETHEAP)
C
*/
{
	double PHI;
/*
C dlp 11/5/98
C      COMMON / MX7 / THETA(40,40)
      COMMON / MX7 / THETA(MAXIONS,MAXIONS)
C
*/
/*
C
C     EQUATION (10B)
C
*/
	PHI=THETA[J,K]+ETHETA();
	return (PHI);
}
/* ---------------------------------------------------------------------- */
double PHIPHI (int J, int K)
/* ---------------------------------------------------------------------- */
/*
C
C     EQUATION (10A)
C
*/
{
	double PHIPHI;
	PHIPHI=THETA[J,K]+ETHETA()+I*ETHEAP();
	return PHIPHI;
}
#endif
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
