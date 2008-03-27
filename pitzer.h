typedef enum
{ TYPE_B0, TYPE_B1, TYPE_B2, TYPE_C0, TYPE_THETA, TYPE_LAMDA, TYPE_ZETA,
    TYPE_PSI, TYPE_ETHETA, TYPE_ALPHAS, TYPE_MU, TYPE_ETA, TYPE_Other } pitz_param_type;


PITZER_EXTERNAL LDBLE VP, DW0;

struct pitz_param
{
  char *species[3];
  int ispec[3];
  pitz_param_type type;
  LDBLE p;
  union
  {
    LDBLE b0;
    LDBLE b1;
    LDBLE b2;
    LDBLE c0;
    LDBLE theta;
    LDBLE lamda;
    LDBLE zeta;
    LDBLE psi;
    LDBLE alphas;
    LDBLE mu;
    LDBLE eta;
  } U;
  LDBLE a[5];
  LDBLE alpha;
  LDBLE os_coef;
  LDBLE ln_coef[3];
  struct theta_param *thetas;
};

PITZER_EXTERNAL struct pitz_param **pitz_params;
PITZER_EXTERNAL int count_pitz_param, max_pitz_param;



/* routines define in pitzer_structures.c */
PITZER_EXTERNAL struct pitz_param *pitz_param_read (char *string, int n);
PITZER_EXTERNAL int pitz_param_search (struct pitz_param *pzp_ptr);
PITZER_EXTERNAL struct theta_param *theta_param_search (LDBLE zj, LDBLE zk);
PITZER_EXTERNAL struct theta_param *theta_param_alloc (void);
PITZER_EXTERNAL int theta_param_init (struct theta_param *theta_param_ptr);



/* defined in DW */
PITZER_EXTERNAL int DW (LDBLE T);
PITZER_EXTERNAL LDBLE DC (LDBLE T);

struct theta_param
{
  LDBLE zj;
  LDBLE zk;
  LDBLE etheta;
  LDBLE ethetap;
};
PITZER_EXTERNAL struct theta_param **theta_params;
PITZER_EXTERNAL int count_theta_param, max_theta_param;
PITZER_EXTERNAL int use_etheta;
PITZER_EXTERNAL LDBLE OTEMP;
