typedef enum { TYPE_B0, TYPE_B1, TYPE_B2, TYPE_C0, TYPE_THETA, TYPE_LAMDA, TYPE_ZETA, TYPE_PSI, TYPE_ETHETA, TYPE_Other } pitz_param_type;


#ifdef PITZER
/* COMMON /MX10/ */
double VP, DW0;
#else
/* COMMON /MX10/ */
extern double VP, DW0;
#endif

struct pitz_param {
	char * species[3];
	int ispec[3];
	pitz_param_type type;
	double p;
	union {double b0; 
		double b1; 
		double b2; 
		double c0; 
		double theta; 
		double lamda;
		double zeta;
		double psi;} U;
	double a[5];
	double alpha;
	struct theta_param *thetas;
};

struct pitz_param **pitz_params;
int count_pitz_param, max_pitz_param;



/* routines define in pitzer_structures.c */
struct pitz_param *pitz_param_read (char *string, int n);
int pitz_param_search(struct pitz_param *pzp_ptr);
struct theta_param *theta_param_search(double zj, double zk);
struct theta_param *theta_param_alloc (void);
int theta_param_init (struct theta_param *theta_param_ptr);



/* defined in DW */
int DW (double T);
double DC (double T);

struct theta_param {
	double zj;
	double zk;
	double etheta;
	double ethetap;
};
struct theta_param **theta_params;
int count_theta_param, max_theta_param;
