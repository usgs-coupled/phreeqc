typedef enum { TYPE_B0, TYPE_B1, TYPE_B2, TYPE_C0, TYPE_THETA, TYPE_LAMDA, TYPE_ZETA, TYPE_PSI, TYPE_Other } pitz_param_type;


#ifdef PITZER
double  GASCON=0.461522e0, TZ=647.073e0, AA=1.e0;
double Z, DZ, Y;
double P[11]={0, 0.7478629e0, -.3540782e0, 0.e0, 0e0, .007159876e0, 0.e0,
	      -.003528426e0, 0., 0., 0.};
double Q[11]={0, 1.1278334e0,0.e0,-.5944001e0,-5.010996e0,0.e0,.63684256e0,
	      0., 0., 0., 0.};
double G1=11.e0,G2=44.333333333333e0,GF=3.5e0;
double B1, B2, B1T, B2T, B1TT, B2TT;
double G[41]={0,-.53062968529023e3,.22744901424408e4,.78779333020687e3
	      ,-.69830527374994e2,.17863832875422e5,-.39514731563338e5
	      ,.33803884280753e5,-.13855050202703e5,-.25637436613260e6
	      ,.48212575981415e6,-.34183016969660e6, .12223156417448e6
	      ,.11797433655832e7,-.21734810110373e7, .10829952168620e7
	      ,-.25441998064049e6,-.31377774947767e7,.52911910757704e7
	      ,-.13802577177877e7,-.25109914369001e6, .46561826115608e7
	      ,-.72752773275387e7,.41774246148294e6,.14016358244614e7
	      ,-.31555231392127e7,.47929666384584e7,.40912664781209e6
	      ,-.13626369388386e7, .69625220862664e6,-.10834900096447e7
	      ,-.22722827401688e6,.38365486000660e6,.68833257944332e4
	      ,.21757245522644e5,-.26627944829770e4,-.70730418082074e5
	      ,-.225e0,-1.68e0,.055e0,-93.0e0};
double ATZ[5]={0,64.e1,64.e1,641.6e0,27.e1},ADZ[5]={0,.319e0,.319e0,.319e0,1.55e0},
	AAT[5]={0,2.e4,2.e4,4.e4,25.e0}, AAD[5]={0,34.e0,4.e1,3.e1,1.05e3};
double VP, DW0;
/* COMMON /QQQQ/ */
double Q0,Q5;
#else
/* common ACONST */
extern double GASCON, TZ, AA;
extern double Z, DZ, Y;
extern double P[11], Q[11];
extern double G1,G2,GF;
extern double B1, B2, B1T, B2T, B1TT, B2TT;
extern double G;
extern double ATZ[4],ADZ[4],AAT[4],AAD[4];
/* common MX10 */
extern double VP, DW0;
/* COMMON /QQQQ/ */
extern double Q0,Q5;
#endif


struct pitz_param {
	char * species[3];
	pitz_param_type type;
	union {double b0; 
		double b1; 
		double b2; 
		double c0; 
		double theta; 
		double lamda;
		double zeta;
		double psi;} U;
	double a[5];
};

struct pitz_param **pitz_params;
int count_pitz_param, max_pitz_param;

struct pitz_param *pitz_param_read (char *string, int n);
int pitz_param_search(struct pitz_param *pzp_ptr);
