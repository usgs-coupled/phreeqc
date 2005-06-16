typedef enum { TYPE_B0, TYPE_B1, TYPE_B2, TYPE_C0, TYPE_THETA, TYPE_LAMDA, TYPE_ZETA, TYPE_PSI, TYPE_Other } pitz_param_type;


#ifdef PITZER

/* COMMON /MX10/ */
double VP, DW0;
#else
/* COMMON /MX10/ */
extern double VP, DW0;
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


/* routines define in pitzer.c */
struct pitz_param *pitz_param_read (char *string, int n);
int pitz_param_search(struct pitz_param *pzp_ptr);

/* defined in DW */
int DW (double T);
double DC (double T);
