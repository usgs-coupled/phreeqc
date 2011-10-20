#ifndef _INC_PITZER_H
#define _INC_PITZER_H

 LDBLE VP, DW0;

/* Pitzer parameters */

 struct pitz_param **pitz_params;
 int count_pitz_param, max_pitz_param;
 struct pitz_param **sit_params;
 int count_sit_param, max_sit_param;

/* routines define in pitzer_structures.c */
 struct pitz_param *pitz_param_read(char *string, int n);
 int pitz_param_search(struct pitz_param *pzp_ptr);
 int sit_param_search(struct pitz_param *pzp_ptr);
 struct theta_param *theta_param_search(LDBLE zj, LDBLE zk);
 struct theta_param *theta_param_alloc(void);
 int theta_param_init(struct theta_param *theta_param_ptr);

/* defined in DW */
 int DW(LDBLE T);
 LDBLE DC(LDBLE T);

 struct theta_param **theta_params;
 int count_theta_param, max_theta_param;
 int use_etheta;
 LDBLE OTEMP;
#endif /* _INC_PITZER_H */
