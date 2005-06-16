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
	double X, JAY;
	BDK (X);
	JAY=X/4.0e0-1.0e0+0.5e0*(BK[0]-BK[2]);
	return JAY;
}
/* ---------------------------------------------------------------------- */
double JPRIME (double Y)
/* ---------------------------------------------------------------------- */
{
	double DZ, Y, JPRIME;
	BDK (Y);
	if (Y > 1.0e0) {
		DZ=-4.0e0*pow(Y,-1.1e0)/9.0e0;
	} else {
		DZ=0.8e0*pow(Y,-0.8e0);
	}			   
	JPRIME=Y*(.25e0+DZ*(DK[0]-DK[2])/2.0e0);
	return JPRIME;
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
     double AKX[42] = { 
         1.925154014814667D0, -.060076477753119D0, -.029779077456514D0,
         -.007299499690937D0, 0.000388260636404D0, 0.000636874599598D0,
         0.000036583601823D0, -.000045036975204D0, -.000004537895710D0,
         0.000002937706971D0, 0.000000396566462D0, -.000000202099617D0,
         -.000000025267769D0, 0.000000013522610D0, 0.000000001229405D0,
         -.000000000821969D0, -.000000000050847D0, 0.000000000046333D0,
         0.000000000001943D0, -.000000000002563D0, -.000000000010991D0,
         0.628023320520852D0, 0.462762985338493D0, 0.150044637187895D0,
         -.028796057604906D0, -.036552745910311D0, -.001668087945272D0,
         0.006519840398744D0, 0.001130378079086D0, -.000887171310131D0,
         -.000242107641309D0, 0.000087294451594D0, 0.000034682122751D0,
         -.000004583768938D0, -.000003548684306D0, -.000000250453880D0,
         0.000000216991779D0, 0.000000080779570D0, 0.000000004558555D0,
         -.000000006944757D0, -.000000002849257D0, 0.000000000237816D0}
/*
      DOUBLE PRECISION AK, BK, DK
      COMMON / MX8 / AK(0:20,2),BK(0:22),DK(0:22)
*/
{
	double Z,X;
	int II, M;
	if (X <= 1.0e0) {
		II=1;
		Z=4.0e0*pow(X,0.2e0)-2.0D0;
		AK = &AKX[0];
	} else {
		II=2;
		Z=40.0e0*pow(X,-1.0e-1)/9.0e0-22.0e0/9.0e0;
		AK = &AKX[21];
	}
	for (i = 20; i >= 0; i--) {
		BK(i)=Z*BK(i+1)-BK(i+2)+AK(i);
		DK(i)=BK(i+1)+Z*DK(i+1)-DK(i+2);
	}
	return OK;
}
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
	PHI=THETA(J,K)+ETHETA();
	return (PHI);
}
/* ---------------------------------------------------------------------- */
double PHI (int J, int K)
/* ---------------------------------------------------------------------- */
/*
C
C     EQUATION (10A)
C
*/
{
	double PHIPHI;
	PHIPHI=THETA(J,K)+ETHETA()+I*ETHEAP();
	return PHIPHI;
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
	double XCON, ZJ, ZK, ZZ, ETHETA;

	if (ZJ == ZK) return(0.0);
	XCON=6.0D0*A0*sqrt(I);
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
	ETHETA=ZZ*(JAY(XJK)-JAY(XJJ)/2.0e0-JAY(XKK)/2.0e0)/(4.0e0*I);
	return (ETHETA);
}
/* ---------------------------------------------------------------------- */
double ETHETAP (double ZJ, double ZK, double I)
/* ---------------------------------------------------------------------- */
{
	double XCON, ZJ, ZK, ZZ, ETHETA, ETHETAP;

	if (ZJ == ZK) return(0.0);
	XCON=6.0D0*A0*sqrt(I);
	ZZ=ZJ*ZK;
/*
C
C     NEXT 3 ARE EQUATION (A1)
C
*/
	XJK=XCON*ZZ;
	XJJ=XCON*ZJ**ZJ;
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

