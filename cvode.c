/*#define DEBUG_CVODE*/
extern char *error_string;
/*******************************************************************
 *                                                                 *
 * File          : cvode.c                                         *
 * Programmers   : Scott D. Cohen, Alan C. Hindmarsh, Radu Serban, *
 *                 and Dan Shumaker @ LLNL                         *
 * Version of    : 24 July 2002                                    *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California * 
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/cvode/LICENSE                         *
 *-----------------------------------------------------------------*
 * This is the implementation file for the main CVODE integrator.  *
 * It is independent of the CVODE linear solver in use.            *
 *                                                                 *
 *******************************************************************/


/************************************************************/
/******************* BEGIN Imports **************************/
/************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "cvode.h"
#include "sundialstypes.h"
#include "nvector.h"
#include "sundialsmath.h"
#include "output.h"
#define KINETICS_EXTERNAL 
#include "kinetics.h"
#include "phqalloc.h"
/* WARNING don't include any headers below here */
#define malloc PHRQ_malloc
static char const svnid[] = "$Id$";

/************************************************************/
/******************** END Imports ***************************/
/************************************************************/


/***************************************************************/
/*********************** BEGIN Macros **************************/
/***************************************************************/

/* Macro: loop */

#define loop for(;;)

/***************************************************************/
/************************ END Macros ***************************/
/***************************************************************/



/************************************************************/
/************** BEGIN CVODE Private Constants ***************/
/************************************************************/

#define FOURTH RCONST(0.25)  /* real 0.25    */
#define THREE RCONST(3.0)    /* real 3.0     */
#define FOUR RCONST(4.0)     /* real 4.0     */
#define HUN RCONST(100.0)    /* real 100.0   */
#define TINY RCONST(1.0e-10) /* small number */
#define HALF   RCONST(0.5)  /* real 0.5   */
#define ZERO   RCONST(0.0)  /* real 0.0   */
#define ONE    RCONST(1.0)  /* real 1.0   */
#define TWO    RCONST(2.0)  /* real 2.0   */
#define TWELVE RCONST(12.0) /* real 12.0  */

/***************************************************************/
/************** BEGIN Default Constants ************************/
/***************************************************************/

#define HMIN_DEFAULT     ZERO    /* hmin default value     */
#define HMAX_INV_DEFAULT ZERO    /* hmax_inv default value */
#define MXHNIL_DEFAULT   10      /* mxhnil default value   */
#define MXSTEP_DEFAULT   50      /* mxstep default value   */


/***************************************************************/
/*************** END Default Constants *************************/
/***************************************************************/


/***************************************************************/
/************ BEGIN Routine-Specific Constants *****************/
/***************************************************************/

/* CVodeDky */

#define FUZZ_FACTOR RCONST(100.0)

/* CVHin */

#define HLB_FACTOR RCONST(100.0)
#define HUB_FACTOR RCONST(0.1)
#define H_BIAS     HALF
#define MAX_ITERS  4

/* CVSet */

#define CORTES RCONST(0.1)

/* CVStep return values */

#define SUCCESS_STEP      0
#define REP_ERR_FAIL     -1
#define REP_CONV_FAIL    -2
#define SETUP_FAILED     -3
#define SOLVE_FAILED     -4

/* CVStep control constants */

#define PREDICT_AGAIN    -5
#define DO_ERROR_TEST     1

/* CVStep */

#define THRESH RCONST(1.5)
#define ETAMX1 RCONST(10000.0) 
#define ETAMX2 RCONST(10.0)
#define ETAMX3 RCONST(10.0)
#define ETAMXF RCONST(0.2)
#define ETAMIN RCONST(0.1)
#define ETACF  RCONST(0.25)
#define ADDON  RCONST(0.000001)
#define BIAS1  RCONST(6.0)
#define BIAS2  RCONST(6.0)
#define BIAS3  RCONST(10.0)
#define ONEPSM RCONST(1.000001)

#define SMALL_NST    10   /* nst > SMALL_NST => use ETAMX3          */
#define MXNCF        10   /* max no. of convergence failures during */
                          /* one step try                           */
#define MXNEF         7   /* max no. of error test failures during  */
                          /* one step try                           */
#define MXNEF1        3   /* max no. of error test failures before  */
                          /* forcing a reduction of order           */
#define SMALL_NEF     2   /* if an error failure occurs and         */
                          /* SMALL_NEF <= nef <= MXNEF1, then       */
                          /* reset eta =  MIN(eta, ETAMXF)          */
#define LONG_WAIT    10   /* number of steps to wait before         */
                          /* considering an order change when       */
                          /* q==1 and MXNEF1 error test failures    */
                          /* have occurred                          */

/* CVnls return values */

#define SOLVED            0
#define CONV_FAIL        -1 
#define SETUP_FAIL_UNREC -2
#define SOLVE_FAIL_UNREC -3

/* CVnls input flags */

#define FIRST_CALL      0
#define PREV_CONV_FAIL -1
#define PREV_ERR_FAIL  -2

/* CVnls other constants */

#define FUNC_MAXCOR 3  /* maximum no. of corrector iterations   */
                       /* for iter == FUNCTIONAL                */
#define NEWT_MAXCOR 3  /* maximum no. of corrector iterations   */
                       /* for iter == NEWTON                    */

#define CRDOWN RCONST(0.3) /* constant used in the estimation of the   */
                           /* convergence rate (crate) of the          */
                           /* iterates for the nonlinear equation      */
#define DGMAX  RCONST(0.3) /* iter == NEWTON, |gamma/gammap-1| > DGMAX */
                           /* => call lsetup                           */

#define RDIV      TWO  /* declare divergence if ratio del/delp > RDIV  */
#define MSBP       20  /* max no. of steps between lsetup calls        */

#define TRY_AGAIN  99  /* control constant for CVnlsNewton - should be */
                       /* distinct from CVnls return values            */


/***************************************************************/
/*************** END Routine-Specific Constants  ***************/
/***************************************************************/


/***************************************************************/
/***************** BEGIN Error Messages ************************/
/***************************************************************/

/* CVodeMalloc/CVReInit Error Messages */

#define CVM             "CVodeMalloc/CVReInit-- "

#define MSG_Y0_NULL     CVM "y0=NULL illegal.\n\n"

#define MSG_BAD_N       CVM "N=%ld < 1 illegal.\n\n"

#define MSG_BAD_LMM_1   CVM "lmm=%d illegal.\n"
#define MSG_BAD_LMM_2   "The legal values are ADAMS=%d and BDF=%d.\n\n"
#define MSG_BAD_LMM     MSG_BAD_LMM_1 MSG_BAD_LMM_2

#define MSG_BAD_ITER_1  CVM "iter=%d illegal.\n"
#define MSG_BAD_ITER_2  "The legal values are FUNCTIONAL=%d "
#define MSG_BAD_ITER_3  "and NEWTON=%d.\n\n"
#define MSG_BAD_ITER    MSG_BAD_ITER_1 MSG_BAD_ITER_2 MSG_BAD_ITER_3

#define MSG_BAD_ITOL_1  CVM "itol=%d illegal.\n"
#define MSG_BAD_ITOL_2  "The legal values are SS=%d and SV=%d.\n\n"
#define MSG_BAD_ITOL    MSG_BAD_ITOL_1 MSG_BAD_ITOL_2

#define MSG_F_NULL       CVM "f=NULL illegal.\n\n"

#define MSG_RELTOL_NULL  CVM "reltol=NULL illegal.\n\n"
 
#define MSG_BAD_RELTOL   CVM "*reltol=%g < 0 illegal.\n\n"

#define MSG_ABSTOL_NULL  CVM "abstol=NULL illegal.\n\n"

#define MSG_BAD_ABSTOL   CVM "Some abstol component < 0.0 illegal.\n\n"

#define MSG_BAD_OPTIN_1  CVM "optIn=%d illegal.\n"
#define MSG_BAD_OPTIN_2  "The legal values are FALSE=%d and TRUE=%d.\n\n"
#define MSG_BAD_OPTIN    MSG_BAD_OPTIN_1 MSG_BAD_OPTIN_2

#define MSG_BAD_OPT      CVM "optIn=TRUE, but iopt=ropt=NULL.\n\n"

#define MSG_MEM_FAIL    CVM "A memory request failed.\n\n"

#define MSG_BAD_EWT     CVM "Some initial ewt component = 0.0 illegal.\n\n"

#define MSG_REI_NO_MEM  "CVReInit-- cvode_mem = NULL illegal.\n\n"

#define MSG_REI_MAXORD1 "CVReInit-- Illegal attempt to increase "
#define MSG_REI_MAXORD2 "maximum method order from %d to %d.\n\n"
#define MSG_REI_MAXORD  MSG_REI_MAXORD1 MSG_REI_MAXORD2 


/* CVode error messages */

#define CVODE            "CVode-- "

#define NO_MEM           "cvode_mem=NULL illegal.\n\n"

#define MSG_CVODE_NO_MEM CVODE NO_MEM
 
#define MSG_LINIT_NULL   CVODE "The linear solver's init routine is NULL.\n\n"

#define MSG_LSETUP_NULL  CVODE "The linear solver's setup routine is NULL.\n\n"

#define MSG_LSOLVE_NULL  CVODE "The linear solver's solve routine is NULL.\n\n"

#define MSG_LFREE_NULL   CVODE "The linear solver's free routine is NULL.\n\n"

#define MSG_LINIT_FAIL   CVODE "The linear solver's init routine failed.\n\n"

#define MSG_YOUT_NULL    CVODE "yout=NULL illegal.\n\n"

#define MSG_T_NULL       CVODE "t=NULL illegal.\n\n"

#define MSG_BAD_ITASK_1   CVODE "itask=%d illegal.\nThe legal values are"
#define MSG_BAD_ITASK_2   " NORMAL=%d and ONE_STEP=%d.\n\n"
#define MSG_BAD_ITASK     MSG_BAD_ITASK_1 MSG_BAD_ITASK_2

#define MSG_BAD_HMIN_HMAX_1 CVODE "Inconsistent step size limits:\n"
#define MSG_BAD_HMIN_HMAX_2 "ropt[HMIN]=%g > ropt[HMAX]=%g.\n\n"
#define MSG_BAD_HMIN_HMAX   MSG_BAD_HMIN_HMAX_1 MSG_BAD_HMIN_HMAX_2

#define MSG_BAD_H0        CVODE "h0=%g and tout-t0=%g inconsistent.\n\n"

#define MSG_BAD_TOUT_1    CVODE "Trouble interpolating at tout = %g.\n"
#define MSG_BAD_TOUT_2    "tout too far back in direction of integration.\n\n"
#define MSG_BAD_TOUT      MSG_BAD_TOUT_1 MSG_BAD_TOUT_2

#define MSG_MAX_STEPS_1   CVODE "At t=%g, mxstep=%d steps taken on "
#define MSG_MAX_STEPS_2   "this call before\nreaching tout=%g.\n\n"
#define MSG_MAX_STEPS     MSG_MAX_STEPS_1 MSG_MAX_STEPS_2

#define MSG_EWT_NOW_BAD_1  CVODE "At t=%g, "
#define MSG_EWT_NOW_BAD_2  "some ewt component has become <= 0.0.\n\n"
#define MSG_EWT_NOW_BAD    MSG_EWT_NOW_BAD_1 MSG_EWT_NOW_BAD_2

#define MSG_TOO_MUCH_ACC  CVODE "At t=%g, too much accuracy requested.\n\n"

#define MSG_HNIL_1  CVODE "Warning.. internal t=%g and step size h=%g\n"
#define MSG_HNIL_2  "are such that t + h == t on the next step.\n"
#define MSG_HNIL_3  "The solver will continue anyway.\n\n"
#define MSG_HNIL    MSG_HNIL_1 MSG_HNIL_2 MSG_HNIL_3

#define MSG_HNIL_DONE_1   CVODE "The above warning has been issued %d times "
#define MSG_HNIL_DONE_2   "and will not be\nissued again for this problem.\n\n"
#define MSG_HNIL_DONE     MSG_HNIL_DONE_1 MSG_HNIL_DONE_2

#define MSG_ERR_FAILS_1   CVODE "At t=%g and step size h=%g, the error test\n"
#define MSG_ERR_FAILS_2   "failed repeatedly or with |h| = hmin.\n\n"
#define MSG_ERR_FAILS     MSG_ERR_FAILS_1 MSG_ERR_FAILS_2

#define MSG_CONV_FAILS_1  CVODE "At t=%g and step size h=%g, the corrector\n"
#define MSG_CONV_FAILS_2  "convergence failed repeatedly or "
#define MSG_CONV_FAILS_3  "with |h| = hmin.\n\n"
#define MSG_CONV_FAILS    MSG_CONV_FAILS_1 MSG_CONV_FAILS_2 MSG_CONV_FAILS_3

#define MSG_SETUP_FAILED_1 CVODE "At t=%g, the setup routine failed in an "
#define MSG_SETUP_FAILED_2 "unrecoverable manner.\n\n"
#define MSG_SETUP_FAILED   MSG_SETUP_FAILED_1 MSG_SETUP_FAILED_2

#define MSG_SOLVE_FAILED_1 CVODE "At t=%g, the solve routine failed in an "
#define MSG_SOLVE_FAILED_2 "unrecoverable manner.\n\n"
#define MSG_SOLVE_FAILED   MSG_SOLVE_FAILED_1 MSG_SOLVE_FAILED_2

#define MSG_TOO_CLOSE_1    CVODE "tout=%g too close to t0=%g to start"
#define MSG_TOO_CLOSE_2    " integration.\n\n"
#define MSG_TOO_CLOSE      MSG_TOO_CLOSE_1 MSG_TOO_CLOSE_2


/* CVodeDky Error Messages */

#define DKY         "CVodeDky-- "

#define MSG_DKY_NO_MEM  DKY NO_MEM

#define MSG_BAD_K   DKY "k=%d illegal.\n\n"

#define MSG_BAD_T_1 DKY "t=%g illegal.\n"
#define MSG_BAD_T_2 "t not in interval tcur-hu=%g to tcur=%g.\n\n"
#define MSG_BAD_T   MSG_BAD_T_1 MSG_BAD_T_2

#define MSG_BAD_DKY DKY "dky=NULL illegal.\n\n"

/***************************************************************/
/****************** END Error Messages *************************/
/***************************************************************/


/************************************************************/
/*************** END CVODE Private Constants ****************/
/************************************************************/


/**************************************************************/
/********* BEGIN Private Helper Functions Prototypes **********/
/**************************************************************/

static booleantype CVAllocVectors(CVodeMem cv_mem, integertype neq, int maxord,
                                  M_Env machEnv);
static void CVFreeVectors(CVodeMem cv_mem, int maxord);

static booleantype CVEwtSet(CVodeMem cv_mem, N_Vector ycur);
static booleantype CVEwtSetSS(CVodeMem cv_mem, N_Vector ycur);
static booleantype CVEwtSetSV(CVodeMem cv_mem, N_Vector ycur);

static booleantype CVHin(CVodeMem cv_mem, realtype tout);
static realtype CVUpperBoundH0(CVodeMem cv_mem, realtype tdist);
static realtype CVYddNorm(CVodeMem cv_mem, realtype hg);

static int  CVStep(CVodeMem cv_mem);

static int CVsldet(CVodeMem cv_mem);

static void CVAdjustParams(CVodeMem cv_mem);
static void CVAdjustOrder(CVodeMem cv_mem, int deltaq);
static void CVAdjustAdams(CVodeMem cv_mem, int deltaq);
static void CVAdjustBDF(CVodeMem cv_mem, int deltaq);
static void CVIncreaseBDF(CVodeMem cv_mem);
static void CVDecreaseBDF(CVodeMem cv_mem);

static void CVRescale(CVodeMem cv_mem);

static void CVPredict(CVodeMem cv_mem);

static void CVSet(CVodeMem cv_mem);
static void CVSetAdams(CVodeMem cv_mem);
static realtype CVAdamsStart(CVodeMem cv_mem, realtype m[]);
static void CVAdamsFinish(CVodeMem cv_mem, realtype m[], realtype M[],
                          realtype hsum);
static realtype CVAltSum(int iend, realtype a[], int k);
static void CVSetBDF(CVodeMem cv_mem);
static void CVSetTqBDF(CVodeMem cv_mem, realtype hsum, realtype alpha0,
                     realtype alpha0_hat, realtype xi_inv, realtype xistar_inv);

static int CVnls(CVodeMem cv_mem, int nflag);
static int CVnlsFunctional(CVodeMem cv_mem);
static int CVnlsNewton(CVodeMem cv_mem, int nflag);
static int CVNewtonIteration(CVodeMem cv_mem);

static int  CVHandleNFlag(CVodeMem cv_mem, int *nflagPtr, realtype saved_t,
                          int *ncfPtr);

static void CVRestore(CVodeMem cv_mem, realtype saved_t);

static booleantype CVDoErrorTest(CVodeMem cv_mem, int *nflagPtr, int *kflagPtr,
                               realtype saved_t, int *nefPtr, realtype *dsmPtr);

static void CVCompleteStep(CVodeMem cv_mem);

static void CVPrepareNextStep(CVodeMem cv_mem, realtype dsm);
static void CVSetEta(CVodeMem cv_mem);
static realtype CVComputeEtaqm1(CVodeMem cv_mem);
static realtype CVComputeEtaqp1(CVodeMem cv_mem);
static void CVChooseEta(CVodeMem cv_mem);
static void CVBDFStab(CVodeMem cv_mem);

static int  CVHandleFailure(CVodeMem cv_mem,int kflag);


/**************************************************************/
/********** END Private Helper Functions Prototypes ***********/
/**************************************************************/


/**************************************************************/
/**************** BEGIN Readability Constants *****************/
/**************************************************************/


#define uround (cv_mem->cv_uround)  
#define zn     (cv_mem->cv_zn) 
#define ewt    (cv_mem->cv_ewt)  
#define y      (cv_mem->cv_y)
#define acor   (cv_mem->cv_acor)
#define tempv  (cv_mem->cv_tempv)
#define ftemp  (cv_mem->cv_ftemp) 
#define q      (cv_mem->cv_q)
#define qprime (cv_mem->cv_qprime)
#define qwait  (cv_mem->cv_qwait)
#define L      (cv_mem->cv_L)
#define h      (cv_mem->cv_h)
#define hprime (cv_mem->cv_hprime)
#define eta    (cv_mem-> cv_eta) 
#define etaqm1 (cv_mem-> cv_etaqm1) 
#define etaq   (cv_mem-> cv_etaq) 
#define etaqp1 (cv_mem-> cv_etaqp1) 
#define nscon  (cv_mem->cv_nscon)
#define ssdat  (cv_mem->cv_ssdat)
#define hscale (cv_mem->cv_hscale) 
#define tn     (cv_mem->cv_tn)
#define tau    (cv_mem->cv_tau)
#define tq     (cv_mem->cv_tq)
#define l      (cv_mem->cv_l)
#define rl1    (cv_mem->cv_rl1)
#define gamma  (cv_mem->cv_gamma) 
#define gammap (cv_mem->cv_gammap) 
#define gamrat (cv_mem->cv_gamrat)
#define crate  (cv_mem->cv_crate)
#define acnrm  (cv_mem->cv_acnrm)
#define mnewt  (cv_mem->cv_mnewt)
#define qmax   (cv_mem->cv_qmax) 
#define mxstep (cv_mem->cv_mxstep)
#define maxcor (cv_mem->cv_maxcor)
#define mxhnil (cv_mem->cv_mxhnil)
#define hmin   (cv_mem->cv_hmin)
#define hmax_inv (cv_mem->cv_hmax_inv)
#define etamax (cv_mem->cv_etamax)
#define nst    (cv_mem->cv_nst)
#define nfe    (cv_mem->cv_nfe)
#define ncfn   (cv_mem->cv_ncfn)
#define netf   (cv_mem->cv_netf)
#define nni    (cv_mem-> cv_nni)
#define nsetups (cv_mem->cv_nsetups)
#define nhnil  (cv_mem->cv_nhnil)
#define lrw    (cv_mem->cv_lrw)
#define liw    (cv_mem->cv_liw)
#define linit  (cv_mem->cv_linit)
#define lsetup (cv_mem->cv_lsetup)
#define lsolve (cv_mem->cv_lsolve) 
#define lfree  (cv_mem->cv_lfree) 
#define lmem   (cv_mem->cv_lmem) 
#define qu     (cv_mem->cv_qu)          
#define nstlp  (cv_mem->cv_nstlp)  
#define hu     (cv_mem->cv_hu)         
#define saved_tq5 (cv_mem->cv_saved_tq5)  
#define jcur   (cv_mem->cv_jcur)         
#define tolsf  (cv_mem->cv_tolsf)      
#define setupNonNull (cv_mem->cv_setupNonNull) 
#define machenv (cv_mem->cv_machenv)
#define sldeton (cv_mem->cv_sldeton)

/**************************************************************/
/***************** END Readability Constants ******************/
/**************************************************************/


/***************************************************************/
/************* BEGIN CVODE Implementation **********************/
/***************************************************************/


/***************************************************************/
/********* BEGIN Exported Functions Implementation *************/
/***************************************************************/


/******************** CVodeMalloc *******************************

 CVodeMalloc allocates and initializes memory for a problem. All
 problem specification inputs are checked for errors. If any
 error occurs during initialization, it is reported to the file
 whose file pointer is errfp and NULL is returned. Otherwise, the
 pointer to successfully initialized problem memory is returned.
 
*****************************************************************/

void *CVodeMalloc(integertype N, RhsFn f, realtype t0, N_Vector y0, 
                  int lmm, int iter, int itol, 
                  realtype *reltol, void *abstol,
                  void *f_data, FILE *errfp, booleantype optIn, 
                  long int iopt[], realtype ropt[], M_Env machEnv)
{
  booleantype allocOK, ioptExists, roptExists, neg_abstol, ewtsetOK;
  int maxord;
  CVodeMem cv_mem;
  FILE *fp;
  int i,k;
  
  if (svnid == NULL) fprintf(stderr," ");
  /* Check for legal input parameters */
  
  fp = (errfp == NULL) ? stdout : errfp;

  if (y0==NULL) {
	  output_msg(OUTPUT_CVODE, MSG_Y0_NULL);
	  return(NULL);
  }
  
  if (N <= 0) {
	  output_msg(OUTPUT_CVODE, MSG_BAD_N, N);
	  return(NULL);
  }

  if ((lmm != ADAMS) && (lmm != BDF)) {
	  output_msg(OUTPUT_CVODE, MSG_BAD_LMM, lmm, ADAMS, BDF);
	  return(NULL);
  }

  if ((iter != FUNCTIONAL) && (iter != NEWTON)) {
	  output_msg(OUTPUT_CVODE, MSG_BAD_ITER, iter, FUNCTIONAL, NEWTON);
	  return(NULL);
  }

  if ((itol != SS) && (itol != SV)) {
	  output_msg(OUTPUT_CVODE, MSG_BAD_ITOL, itol, SS, SV);
	  return(NULL);
  }

  if (f == NULL) {
	  output_msg(OUTPUT_CVODE, MSG_F_NULL);
	  return(NULL);
  }

  if (reltol == NULL) {
	  output_msg(OUTPUT_CVODE, MSG_RELTOL_NULL);
	  return(NULL);
  }

  if (*reltol < ZERO) {
	  output_msg(OUTPUT_CVODE, MSG_BAD_RELTOL, (double) *reltol);
	  return(NULL);
  }
   
  if (abstol == NULL) {
	  output_msg(OUTPUT_CVODE, MSG_ABSTOL_NULL);
	  return(NULL);
  }

  if (itol == SS) {
    neg_abstol = (*((realtype *)abstol) < ZERO);
  } else {
    neg_abstol = (N_VMin((N_Vector)abstol) < ZERO);
  }
  if (neg_abstol) {
	  output_msg(OUTPUT_CVODE, MSG_BAD_ABSTOL);
	  return(NULL);
  }

  if ((optIn != FALSE) && (optIn != TRUE)) {
	  output_msg(OUTPUT_CVODE, MSG_BAD_OPTIN, optIn, FALSE, TRUE);
	  return(NULL);
  }

  if ((optIn) && (iopt == NULL) && (ropt == NULL)) {
	  output_msg(OUTPUT_CVODE, MSG_BAD_OPT);
	  return(NULL);
  } 

  ioptExists = (iopt != NULL);
  roptExists = (ropt != NULL);

  /* Compute maxord */

  maxord = (lmm == ADAMS) ? ADAMS_Q_MAX : BDF_Q_MAX;

  if (optIn && ioptExists) {
    if (iopt[MAXORD] > 0)  maxord = MIN(maxord, iopt[MAXORD]);
  }

  cv_mem = (CVodeMem) malloc(sizeof(struct CVodeMemRec));
  if (cv_mem == NULL) {
	  output_msg(OUTPUT_CVODE, MSG_MEM_FAIL);
	  return(NULL);
  }
 
  /* Allocate the vectors */

  allocOK = CVAllocVectors(cv_mem, N, maxord, machEnv);
  if (!allocOK) {
	  output_msg(OUTPUT_CVODE, MSG_MEM_FAIL);
	  free(cv_mem);
	  return(NULL);
  }

  /* Copy tolerances into memory, and set the ewt vector */

  cv_mem->cv_itol = itol;
  cv_mem->cv_reltol = reltol;      
  cv_mem->cv_abstol = abstol;
  ewtsetOK = CVEwtSet(cv_mem, y0);
  if (!ewtsetOK) {
	  output_msg(OUTPUT_CVODE, MSG_BAD_EWT);
	  CVFreeVectors(cv_mem, maxord);
	  free(cv_mem);
	  return(NULL);
  }
  
  /* All error checking is complete at this point */

  /* Copy the remaining input parameters into CVODE memory */

  cv_mem->cv_N = N;
  cv_mem->cv_f = f;
  cv_mem->cv_f_data = f_data;
  cv_mem->cv_lmm = lmm;    
  cv_mem->cv_iter = iter;
  cv_mem->cv_optIn = optIn;
  cv_mem->cv_iopt = iopt;
  cv_mem->cv_ropt = ropt;
  cv_mem->cv_errfp = fp;
  tn = t0;
  machenv = machEnv;

  /* Set step parameters */

  q = 1;
  L = 2;
  qwait = L;
  qmax = maxord;
  etamax = ETAMX1;

  /* Set uround */

  uround = UnitRoundoff();

  /* Set the linear solver addresses to NULL.
     (We check != NULL later, in CVode, if using NEWTON.) */

  linit = NULL;
  lsetup = NULL;
  lsolve = NULL;
  lfree = NULL;
  lmem = NULL;

  /* Initialize zn[0] in the history array */
  
  N_VScale(ONE, y0, zn[0]);
 
  /* Handle the remaining optional inputs (CVode checks ropt[HMAX]) */

  hmax_inv = HMAX_INV_DEFAULT;
  hmin = HMIN_DEFAULT;
  if (optIn && roptExists) {
    if (ropt[HMIN] > ZERO) hmin = ropt[HMIN];
  }

  mxhnil = MXHNIL_DEFAULT;
  mxstep = MXSTEP_DEFAULT;
  if (optIn && ioptExists) {
    if (iopt[MXHNIL] != 0) mxhnil = iopt[MXHNIL];
    if (iopt[MXSTEP] > 0) mxstep = iopt[MXSTEP];
  }
 
  if ((!optIn) && roptExists) ropt[H0] = ZERO;
 
  /* Set maxcor */

  maxcor = (iter==NEWTON) ? NEWT_MAXCOR : FUNC_MAXCOR;
  
  /* Initialize all the counters */
 
  nst = nfe = ncfn = netf = nni = nsetups = nhnil = nstlp = 0;
  
  /* Initialize all other variables corresponding to optional outputs */
  
  qu = 0;
  hu = ZERO;
  tolsf = ONE;

  /* Initialize optional output locations in iopt, ropt */
  /* and  Stablilty Limit Detection data.               */

  nscon = 0;
  sldeton = FALSE;
  if (ioptExists) {
    iopt[NST] = iopt[NFE] = iopt[NSETUPS] = iopt[NNI] = 0;
    iopt[NCFN] = iopt[NETF] = 0;
    iopt[QU] = qu;
    iopt[QCUR] = 0;
    iopt[LENRW] = lrw;
    iopt[LENIW] = liw;
    if(optIn && iopt[SLDET] && (lmm == BDF)) {
      sldeton = TRUE;
      iopt[NOR] = 0;
      for (i = 1; i <= 5; i++) {
        for (k = 1; k <= 3; k++) ssdat[i-1][k-1] = ZERO;}
    }
  }
  
  if (roptExists) {
    ropt[HU] = hu;
    ropt[HCUR] = ZERO;
    ropt[TCUR] = t0;
    ropt[TOLSF] = tolsf;
  }
      

  /* Problem has been successfully initialized */

  return((void *)cv_mem);
}


/******************** CVReInit **********************************

 CVReInit re-initializes CVODE's memory for a problem, assuming
 it has already been allocated in a prior CVodeMalloc call.
 All problem specification inputs are checked for errors.
 The problem size N is assumed to be unchanged since the call to
 CVodeMalloc, and the maximum order maxord must not be larger.
 If any error occurs during initialization, it is reported to the
 file whose file pointer is errfp.
 The return value is SUCCESS = 0 if no errors occurred, or
 a negative value otherwise.
 
*****************************************************************/

int CVReInit(void *cvode_mem, RhsFn f, realtype t0, N_Vector y0,
             int lmm, int iter, int itol, 
             realtype *reltol, void *abstol,
             void *f_data,  FILE *errfp, booleantype optIn, 
             long int iopt[], realtype ropt[], M_Env machEnv) 
{
  booleantype ioptExists, roptExists, neg_abstol, ewtsetOK;
  int maxord, i,k;
  CVodeMem cv_mem;
  FILE *fp;
 
  /* Check for legal input parameters */
  
  fp = (errfp == NULL) ? stdout : errfp;

  if (cvode_mem == NULL) {
	  output_msg(OUTPUT_CVODE, MSG_REI_NO_MEM);
	  return(CVREI_NO_MEM);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (y0 == NULL) {
	  output_msg(OUTPUT_CVODE, MSG_Y0_NULL);
	  return(CVREI_ILL_INPUT);
  }
  
  if ((lmm != ADAMS) && (lmm != BDF)) {
	  output_msg(OUTPUT_CVODE, MSG_BAD_LMM, lmm, ADAMS, BDF);
	  return(CVREI_ILL_INPUT);
  }

  if ((iter != FUNCTIONAL) && (iter != NEWTON)) {
	  output_msg(OUTPUT_CVODE, MSG_BAD_ITER, iter, FUNCTIONAL, NEWTON);
	  return(CVREI_ILL_INPUT);
  }

  if ((itol != SS) && (itol != SV)) {
	  output_msg(OUTPUT_CVODE, MSG_BAD_ITOL, itol, SS, SV);
	  return(CVREI_ILL_INPUT);
  }

  if (f == NULL) {
	  output_msg(OUTPUT_CVODE, MSG_F_NULL);
	  return(CVREI_ILL_INPUT);
  }

  if (reltol == NULL) {
	  output_msg(OUTPUT_CVODE, MSG_RELTOL_NULL);
	  return(CVREI_ILL_INPUT);
  }

  if (*reltol < ZERO) {
	  output_msg(OUTPUT_CVODE, MSG_BAD_RELTOL,(double)  *reltol);
	  return(CVREI_ILL_INPUT);
  }
   
  if (abstol == NULL) {
	  output_msg(OUTPUT_CVODE, MSG_ABSTOL_NULL);
	  return(CVREI_ILL_INPUT);
  }

  if (itol == SS) {
    neg_abstol = (*((realtype *)abstol) < ZERO);
  } else {
    neg_abstol = (N_VMin((N_Vector)abstol) < ZERO);
  }
  if (neg_abstol) {
	  output_msg(OUTPUT_CVODE, MSG_BAD_ABSTOL);
	  return(CVREI_ILL_INPUT);
  }

  if ((optIn != FALSE) && (optIn != TRUE)) {
	  output_msg(OUTPUT_CVODE, MSG_BAD_OPTIN, optIn, FALSE, TRUE);
	  return(CVREI_ILL_INPUT);
  }

  if ((optIn) && (iopt == NULL) && (ropt == NULL)) {
	  output_msg(OUTPUT_CVODE, MSG_BAD_OPT);
	  return(CVREI_ILL_INPUT);
  } 

  ioptExists = (iopt != NULL);
  roptExists = (ropt != NULL);

  /* Compute new maxord and check against old value */

  maxord = (lmm == ADAMS) ? ADAMS_Q_MAX : BDF_Q_MAX;
  if (optIn && ioptExists)
    { if (iopt[MAXORD] > 0)  maxord = MIN(maxord, iopt[MAXORD]); }
  if (maxord > qmax) {
	  output_msg(OUTPUT_CVODE, MSG_REI_MAXORD, qmax, maxord);
	  return(CVREI_ILL_INPUT);
  }

  /* Copy tolerances into memory, and set the ewt vector */

  cv_mem->cv_itol = itol;
  cv_mem->cv_reltol = reltol;      
  cv_mem->cv_abstol = abstol;
  ewtsetOK = CVEwtSet(cv_mem, y0);
  if (!ewtsetOK) {
	  output_msg(OUTPUT_CVODE, MSG_BAD_EWT);
	  return(CVREI_ILL_INPUT);
  }
  
  /* All error checking is complete at this point */
  
  /* Copy the remaining input parameters into CVODE memory */

  cv_mem->cv_f = f;
  cv_mem->cv_f_data = f_data;
  cv_mem->cv_lmm = lmm;    
  cv_mem->cv_iter = iter;
  cv_mem->cv_optIn = optIn;
  cv_mem->cv_iopt = iopt;
  cv_mem->cv_ropt = ropt;
  cv_mem->cv_errfp = fp;
  tn = t0;
  machenv = machEnv;

  /* Set step parameters */

  q = 1;
  L = 2;
  qwait = L;
  qmax = maxord;
  etamax = ETAMX1;

  /* Set uround */

  uround = UnitRoundoff();

  /* Initialize zn[0] in the history array */
  
  N_VScale(ONE, y0, zn[0]);
 
  /* Handle the remaining optional inputs (CVode checks ropt[HMAX]) */

  hmax_inv = HMAX_INV_DEFAULT;
  hmin = HMIN_DEFAULT;
  if (optIn && roptExists) {
    if (ropt[HMIN] > ZERO) hmin = ropt[HMIN];
  }

  mxhnil = MXHNIL_DEFAULT;
  mxstep = MXSTEP_DEFAULT;
  if (optIn && ioptExists) {
    if (iopt[MXHNIL] != 0) mxhnil = iopt[MXHNIL];
    if (iopt[MXSTEP] > 0) mxstep = iopt[MXSTEP];
  }
 
  if ((!optIn) && roptExists) ropt[H0] = ZERO;
 
  /* Set maxcor */

  maxcor = (iter==NEWTON) ? NEWT_MAXCOR : FUNC_MAXCOR;
  
  /* Initialize all the counters */
 
  nst = nfe = ncfn = netf = nni = nsetups = nhnil = nstlp = 0;
  
  /* Initialize all other vars corresponding to optional outputs */
  
  qu = 0;
  hu = ZERO;
  tolsf = ONE;

  /* Initialize optional output locations in iopt, ropt */
  /* and  Stablilty Limit Detection data.               */

  nscon = 0;
  sldeton = FALSE;
  if (ioptExists) {
    iopt[NST] = iopt[NFE] = iopt[NSETUPS] = iopt[NNI] = 0;
    iopt[NCFN] = iopt[NETF] = 0;
    iopt[QU] = qu;
    iopt[QCUR] = 0;
    iopt[LENRW] = lrw;
    iopt[LENIW] = liw;
    if(optIn && iopt[SLDET] && (lmm == BDF)) {
      sldeton = TRUE;
      iopt[NOR] = 0;
      for (i = 1; i <= 5; i++) {
        for (k = 1; k <= 3; k++) ssdat[i-1][k-1] = ZERO;}
    }
  }
  
  if (roptExists) {
    ropt[HU] = hu;
    ropt[HCUR] = ZERO;
    ropt[TCUR] = t0;
    ropt[TOLSF] = tolsf;
  }
      
  /* Problem has been successfully re-initialized */

  return(SUCCESS);
}


/**************************************************************/
/************** BEGIN More Readability Constants **************/
/**************************************************************/

#define N      (cv_mem->cv_N)
#define f      (cv_mem->cv_f)      
#define f_data (cv_mem->cv_f_data)    
#define lmm    (cv_mem->cv_lmm) 
#define iter   (cv_mem->cv_iter)        
#define itol   (cv_mem->cv_itol)         
#define reltol (cv_mem->cv_reltol)       
#define abstol (cv_mem->cv_abstol)     
#define optIn  (cv_mem->cv_optIn)
#define iopt   (cv_mem->cv_iopt)
#define ropt   (cv_mem->cv_ropt)
#define errfp  (cv_mem->cv_errfp)

/**************************************************************/
/*************** END More Readability Constants ***************/
/**************************************************************/


/********************* CVode ****************************************

 This routine is the main driver of the CVODE package. 

 It integrates over a time interval defined by the user, by calling
 CVStep to do internal time steps.

 The first time that CVode is called for a successfully initialized
 problem, it computes a tentative initial step size h.

 CVode supports two modes, specified by itask: NORMAL and ONE_STEP.
 In the NORMAL mode, the solver steps until it reaches or passes tout
 and then interpolates to obtain y(tout).
 In the ONE_STEP mode, it takes one internal step and returns.

********************************************************************/

int CVode(void *cvode_mem, realtype tout, N_Vector yout, 
          realtype *t, int itask)
{
  int nstloc, kflag, istate, next_q, ier;
  realtype rh, next_h;
  booleantype hOK, ewtsetOK;
  CVodeMem cv_mem;
  realtype t0;

  /* Check for legal inputs in all cases */

  cv_mem = (CVodeMem) cvode_mem;
  if (cvode_mem == NULL) {
	  output_msg(OUTPUT_CVODE, MSG_CVODE_NO_MEM);
	  return(CVODE_NO_MEM);
  }
  
  if ((y = yout) == NULL) {
	  output_msg(OUTPUT_CVODE, MSG_YOUT_NULL);       
	  return(ILL_INPUT);
  }
  
  if (t == NULL) {
	  output_msg(OUTPUT_CVODE, MSG_T_NULL);
	  return(ILL_INPUT);
  }
  t0 = tn;
  *t = tn;

  if ((itask != NORMAL) && (itask != ONE_STEP)) {
	  output_msg(OUTPUT_CVODE, MSG_BAD_ITASK, itask, NORMAL, ONE_STEP);
	  return(ILL_INPUT);
  }

  /* Set hmax_inv from ropt[HMAX] and test for hmin > hmax */

  if (optIn && ropt != NULL) {
    if (ropt[HMAX] > ZERO) hmax_inv = ONE/ropt[HMAX];
    if (hmin*hmax_inv > ONE) {
	    output_msg(OUTPUT_CVODE, MSG_BAD_HMIN_HMAX, (double) hmin, (double) ropt[HMAX]);
	    return(ILL_INPUT);
    }
  }

  /* On first call, check solver functions and call linit function */
  
  if (nst == 0) {
    if (iter == NEWTON) {
      if (linit == NULL) {
	      output_msg(OUTPUT_CVODE, MSG_LINIT_NULL);
	      return(ILL_INPUT);
      }
      if (lsetup == NULL) {
	      output_msg(OUTPUT_CVODE, MSG_LSETUP_NULL);
	      return(ILL_INPUT);
      }
      if (lsolve == NULL) {
	      output_msg(OUTPUT_CVODE, MSG_LSOLVE_NULL);
	      return(ILL_INPUT);
      }
      if (lfree == NULL) {
	      output_msg(OUTPUT_CVODE, MSG_LFREE_NULL);
	      return(ILL_INPUT);
      }
      ier = linit(cv_mem);
      if (ier != LINIT_OK) {
	      output_msg(OUTPUT_CVODE, MSG_LINIT_FAIL);
	      return(ILL_INPUT);
      }
    }
    
    /* On the first call, call f at (t0,y0), set zn[1] = y'(t0), 
       set initial h (from H0 or CVHin), and scale zn[1] by h   */
    cvode_rate_sim_time = cvode_rate_sim_time_start + tn;
    cvode_step_fraction = 0;

    f(N, tn, zn[0], zn[1], f_data); 
    nfe = 1;
    h = ZERO;
    if (ropt != NULL) h = ropt[H0];
    if ( (h != ZERO) && ((tout-tn)*h < ZERO) ) {
	    output_msg(OUTPUT_CVODE, MSG_BAD_H0, (double) h, (double) (tout-tn));
	    return(ILL_INPUT);
    }
    if (h == ZERO) {
      hOK = CVHin(cv_mem, tout);
      if (!hOK) {
	      output_msg(OUTPUT_CVODE, MSG_TOO_CLOSE, (double) tout, (double) tn);
	      return(ILL_INPUT);
      }
    }
    rh = ABS(h)*hmax_inv;
    if (rh > ONE) h /= rh;
    if (ABS(h) < hmin) h *= hmin/ABS(h);
    hscale = h; 
    N_VScale(h, zn[1], zn[1]);

  } /* end of first call block */

  /* If not the first call, check if tout already reached */

  if ( (itask == NORMAL) && (nst > 0) && ((tn-tout)*h >= ZERO) ) {
    *t = tout;
    ier =  CVodeDky(cv_mem, tout, 0, yout);
    if (ier != OKAY) {  /* ier must be == BAD_T */
	    output_msg(OUTPUT_CVODE, MSG_BAD_TOUT, (double) tout);
	    return(ILL_INPUT);
    }
    return(SUCCESS);
  }

  /* Looping point for internal steps */

  nstloc = 0;
  loop {
   
    next_h = h;
    next_q = q;
    
    /* Reset and check ewt */

    if (nst > 0) {
      ewtsetOK = CVEwtSet(cv_mem, zn[0]);
      if (!ewtsetOK) {
	      output_msg(OUTPUT_CVODE, MSG_EWT_NOW_BAD, (double) tn);
	      istate = ILL_INPUT;
	      *t = tn;
	      N_VScale(ONE, zn[0], yout);
	      break;
      }
    }
    
    /* Check for too many steps */
    
    if (nstloc >= mxstep) {
	    /* output_msg(OUTPUT_CVODE, MSG_MAX_STEPS, tn, mxstep, tout); */
	    istate = TOO_MUCH_WORK;
	    *t = tn;
	    N_VScale(ONE, zn[0], yout);
	    break;
    }

    /* Check for too much accuracy requested */

    if ((tolsf = uround * N_VWrmsNorm(zn[0], ewt)) > ONE) {
	    output_msg(OUTPUT_CVODE, MSG_TOO_MUCH_ACC, (double) tn);
	    istate = TOO_MUCH_ACC;
	    *t = tn;
	    N_VScale(ONE, zn[0], yout);
	    tolsf *= TWO;
	    break;
    }

    /* Check for h below roundoff level in tn */

    if (tn + h == tn) {
      nhnil++;
      if (nhnil <= mxhnil) output_msg(OUTPUT_CVODE, MSG_HNIL, (double) tn, (double) h);
      if (nhnil == mxhnil) output_msg(OUTPUT_CVODE, MSG_HNIL_DONE, mxhnil);
    }

    /* Call CVStep to take a step */

    kflag = CVStep(cv_mem);
#ifdef DEBUG_CVODE
    cvode_test = TRUE;
    f(N, tn, y, ftemp, f_data);
    cvode_test = FALSE;
    if (cvode_error == TRUE) {
	    output_msg(OUTPUT_CVODE,"After CVStep, y Fail\n");
    } else {
	    output_msg(OUTPUT_CVODE,"After CVStep, y OK\n");
    }
    cvode_test = TRUE;
    f(N, tn, zn[0], ftemp, f_data);
    cvode_test = FALSE;
    if (cvode_error == TRUE) {
	    output_msg(OUTPUT_CVODE,"After CVStep, zn Fail\n");
    } else {
	    output_msg(OUTPUT_CVODE,"After CVStep, zn OK\n");
    }
#endif
    /* Process failed step cases, and exit loop */
   
    if (kflag != SUCCESS_STEP) {
      istate = CVHandleFailure(cv_mem, kflag);
      *t = tn;
      N_VScale(ONE, zn[0], yout);
      break;
    }
    
    nstloc++;

    /* Check if in one-step mode, and if so copy y and exit loop */
    
    if (itask == ONE_STEP) {
      istate = SUCCESS;
      *t = tn;
      N_VScale(ONE, zn[0], yout);
      next_q = qprime;
      next_h = hprime;
      break;
    }
    cvode_rate_sim_time = cvode_rate_sim_time_start + tn;
    cvode_step_fraction = (tn - t0)/(tout - t0);
    /*
    output_msg(OUTPUT_CVODE, "ODE: tn %e, t0 %e, tout %e, step_frac %e\n", (double) tn, (double) t0, (double) tout, (double) cvode_step_fraction);
    */
    /* Check if tout reached, and if so interpolate and exit loop */

    if ((tn-tout)*h >= ZERO) {
	    /*
	      output_msg(OUTPUT_CVODE, "*tn %e, t0 %e, tout %e, h %e\n", tn, t0, tout,h);
	    */
	    cvode_rate_sim_time = cvode_rate_sim_time_start + tout;
	    cvode_step_fraction = 1.0;
      istate = SUCCESS;
      *t = tout;
      (void) CVodeDky(cv_mem, tout, 0, yout);
      next_q = qprime;
      next_h = hprime;
      break;
    }
  }

  /* End of step loop; load optional outputs and return */

  if (iopt != NULL) {
    iopt[NST] = nst;
    iopt[NFE] = nfe;
    iopt[NSETUPS] = nsetups;
    iopt[NNI] = nni;
    iopt[NCFN] = ncfn;
    iopt[NETF] = netf;
    iopt[QU] = q;
    iopt[QCUR] = next_q;
  }
  
  if (ropt != NULL) {
    ropt[HU] = h;
    ropt[HCUR] = next_h;
    ropt[TCUR] = tn;
    ropt[TOLSF] = tolsf;
  }
#ifdef DEBUG_CVODE
  /*
   * check interpolation
   */
    cvode_test = TRUE;
    f(N, tn, y, ftemp, f_data);
    cvode_test = FALSE;
    if (cvode_error == TRUE) {
	    output_msg(OUTPUT_CVODE,"End of cvode, Interpolated y Fail\n");
	    return(-1);
    } else {
	    output_msg(OUTPUT_CVODE,"End of cvode, Interpolated y OK\n");
    }
#endif
  return(istate);
}

/*************** CVodeDky ********************************************

 This routine computes the k-th derivative of the interpolating
 polynomial at the time t and stores the result in the vector dky.
 The formula is:
          q 
   dky = SUM c(j,k) * (t - tn)^(j-k) * h^(-j) * zn[j] , 
         j=k 
 where c(j,k) = j*(j-1)*...*(j-k+1), q is the current order, and
 zn[j] is the j-th column of the Nordsieck history array.

 This function is called by CVode with k = 0 and t = tout, but
 may also be called directly by the user.

**********************************************************************/

int CVodeDky(void *cvode_mem, realtype t, int k, N_Vector dky)
{
  realtype s, c, r;
  realtype tfuzz, tp, tn1;
  int i, j;
  CVodeMem cv_mem;
  
  cv_mem = (CVodeMem) cvode_mem;

  /* Check all inputs for legality */
 
  if (cvode_mem == NULL) {
	  output_msg(OUTPUT_CVODE, MSG_DKY_NO_MEM);
	  return(DKY_NO_MEM);
  }
  
  if (dky == NULL) {
	  output_msg(OUTPUT_CVODE, MSG_BAD_DKY);
	  return(BAD_DKY);
  }

  if ((k < 0) || (k > q)) {
	  output_msg(OUTPUT_CVODE, MSG_BAD_K, k);
	  return(BAD_K);
  }
  
  tfuzz = FUZZ_FACTOR * uround * (ABS(tn) + ABS(hu));
  if (hu < ZERO) tfuzz = -tfuzz;
  tp = tn - hu - tfuzz;
  tn1 = tn + tfuzz;
  if ((t-tp)*(t-tn1) > ZERO) {
	  output_msg(OUTPUT_CVODE, MSG_BAD_T, (double) t, (double) (tn-hu), (double) tn);
	  return(BAD_T);
  }

  /* Sum the differentiated interpolating polynomial */

  s = (t - tn) / h;
  for (j=q; j >= k; j--) {
    c = ONE;
    for (i=j; i >= j-k+1; i--) c *= i;
    if (j == q) {
      N_VScale(c, zn[q], dky);
    } else {
      N_VLinearSum(c, zn[j], s, dky, dky);
    }
  }
  if (k == 0) return(OKAY);
  r = RPowerI(h,-k);
  N_VScale(r, dky, dky);
  return(OKAY);
}
 
/********************* CVodeFree **********************************

 This routine frees the problem memory allocated by CVodeMalloc.
 Such memory includes all the vectors allocated by CVAllocVectors,
 and the memory lmem for the linear solver (deallocated by a call
 to lfree).

*******************************************************************/

void CVodeFree(void *cvode_mem)
{
  CVodeMem cv_mem;

  cv_mem = (CVodeMem) cvode_mem;
  
  if (cvode_mem == NULL) return;

  CVFreeVectors(cv_mem, qmax);
  if (iter == NEWTON) lfree(cv_mem);
  free(cv_mem);
}


/***************************************************************/
/********** END Exported Functions Implementation **************/
/***************************************************************/


/*******************************************************************/
/******** BEGIN Private Helper Functions Implementation ************/
/*******************************************************************/
 
/****************** CVAllocVectors ***********************************

 This routine allocates the CVODE vectors ewt, acor, tempv, ftemp, and
 zn[0], ..., zn[maxord]. The length of the vectors is the input
 parameter neq and the maximum order (needed to allocate zn) is the
 input parameter maxord. If all memory allocations are successful,
 CVAllocVectors returns TRUE. Otherwise all allocated memory is freed
 and CVAllocVectors returns FALSE.
 This routine also sets the optional outputs lrw and liw, which are
 (respectively) the lengths of the real and integer work spaces
 allocated here.

**********************************************************************/

static booleantype CVAllocVectors(CVodeMem cv_mem, integertype neq, 
                                  int maxord, M_Env machEnv)
{
  int i, j;

  /* Allocate ewt, acor, tempv, ftemp */
  
  ewt = N_VNew(neq, machEnv);
  if (ewt == NULL) return(FALSE);
  acor = N_VNew(neq, machEnv);
  if (acor == NULL) {
    N_VFree(ewt);
    return(FALSE);
  }
  tempv = N_VNew(neq, machEnv);
  if (tempv == NULL) {
    N_VFree(ewt);
    N_VFree(acor);
    return(FALSE);
  }
  ftemp = N_VNew(neq, machEnv);
  if (ftemp == NULL) {
    N_VFree(tempv);
    N_VFree(ewt);
    N_VFree(acor);
    return(FALSE);
  }

  /* Allocate zn[0] ... zn[maxord] */

  for (j=0; j <= maxord; j++) {
    zn[j] = N_VNew(neq, machEnv);
    if (zn[j] == NULL) {
      N_VFree(ewt);
      N_VFree(acor);
      N_VFree(tempv);
      N_VFree(ftemp);
      for (i=0; i < j; i++) N_VFree(zn[i]);
      return(FALSE);
    }
  }

  /* Set solver workspace lengths  */

  lrw = (maxord + 5)*neq;
  liw = 0;

  return(TRUE);
}

/***************** CVFreeVectors *********************************
  
 This routine frees the CVODE vectors allocated in CVAllocVectors.

******************************************************************/

static void CVFreeVectors(CVodeMem cv_mem, int maxord)
{
  int j;
  
  N_VFree(ewt);
  N_VFree(acor);
  N_VFree(tempv);
  N_VFree(ftemp);
  for(j=0; j <= maxord; j++) N_VFree(zn[j]);
}

/*********************** CVEwtSet **************************************
  
 This routine is responsible for setting the error weight vector ewt,
 according to tol_type, as follows:

 (1) ewt[i] = 1 / (*reltol * ABS(ycur[i]) + *abstol), i=0,...,neq-1
     if tol_type = SS
 (2) ewt[i] = 1 / (*reltol * ABS(ycur[i]) + abstol[i]), i=0,...,neq-1
     if tol_type = SV

  CVEwtSet returns TRUE if ewt is successfully set as above to a
  positive vector and FALSE otherwise. In the latter case, ewt is
  considered undefined after the FALSE return from CVEwtSet.

  All the real work is done in the routines CVEwtSetSS, CVEwtSetSV.
 
***********************************************************************/

static booleantype CVEwtSet(CVodeMem cv_mem, N_Vector ycur)
{
  switch(itol) {
  case SS: return(CVEwtSetSS(cv_mem, ycur));
  case SV: return(CVEwtSetSV(cv_mem, ycur));
  }
  return(-99);
}

/*********************** CVEwtSetSS *********************************

 This routine sets ewt as decribed above in the case tol_type = SS.
 It tests for non-positive components before inverting. CVEwtSetSS
 returns TRUE if ewt is successfully set to a positive vector
 and FALSE otherwise. In the latter case, ewt is considered
 undefined after the FALSE return from CVEwtSetSS.

********************************************************************/

static booleantype CVEwtSetSS(CVodeMem cv_mem, N_Vector ycur)
{
  realtype rtoli, atoli;
  
  rtoli = *reltol;
  atoli = *((realtype *)abstol);
  N_VAbs(ycur, tempv);
  N_VScale(rtoli, tempv, tempv);
  N_VAddConst(tempv, atoli, tempv);
  if (N_VMin(tempv) <= ZERO) return(FALSE);
  N_VInv(tempv, ewt);
  return(TRUE);
}

/*********************** CVEwtSetSV *********************************

 This routine sets ewt as decribed above in the case tol_type = SV.
 It tests for non-positive components before inverting. CVEwtSetSV
 returns TRUE if ewt is successfully set to a positive vector
 and FALSE otherwise. In the latter case, ewt is considered
 undefined after the FALSE return from CVEwtSetSV.

********************************************************************/

static booleantype CVEwtSetSV(CVodeMem cv_mem, N_Vector ycur)
{
  realtype rtoli;
  rtoli = *reltol;
  N_VAbs(ycur, tempv);
  N_VLinearSum(rtoli, tempv, ONE, (N_Vector) abstol, tempv);
  if (N_VMin(tempv) <= ZERO) return(FALSE);
  N_VInv(tempv, ewt);
  return(TRUE);
}

/******************* CVHin ***************************************

 This routine computes a tentative initial step size h0. 
 If tout is too close to tn (= t0), then CVHin returns FALSE and
 h remains uninitialized. Otherwise, CVHin sets h to the chosen 
 value h0 and returns TRUE.

 The algorithm used seeks to find h0 as a solution of
       (WRMS norm of (h0^2 ydd / 2)) = 1, 
 where ydd = estimated second derivative of y.

*****************************************************************/

static booleantype CVHin(CVodeMem cv_mem, realtype tout)
{
  int sign, count;
  realtype tdiff, tdist, tround, hlb, hub;
  realtype hg, hgs, hnew, hrat, h0, yddnrm;

  /* Test for tout too close to tn */
  
  if ((tdiff = tout-tn) == ZERO) return(FALSE);
  
  sign = (tdiff > ZERO) ? 1 : -1;
  tdist = ABS(tdiff);
  tround = uround * MAX(ABS(tn), ABS(tout));
  if (tdist < TWO*tround) return(FALSE);
  
  /* Set lower and upper bounds on h0, and take geometric mean 
     Exit with this value if the bounds cross each other       */

  hlb = HLB_FACTOR * tround;
  hub = CVUpperBoundH0(cv_mem, tdist);
  hg  = RSqrt(hlb*hub);
  if (hub < hlb) {
    if (sign == -1) hg = -hg;
    h = hg;
    return(TRUE);
  }
  
  /* Loop up to MAX_ITERS times to find h0.
     Stop if new and previous values differ by a factor < 2.
     Stop if hnew/hg > 2 after one iteration, as this probably means
     that the ydd value is bad because of cancellation error.        */

  count = 0;
  loop {
    hgs = hg*sign;
    yddnrm = CVYddNorm(cv_mem, hgs);
    if (cvode_error == TRUE) {
	    hg /= 2.;
#ifdef DEBUG_CVODE
	    output_msg(OUTPUT_CVODE, "halving step in CVHin\n");
#endif
	    continue;
    }

    hnew =  (yddnrm*hub*hub > TWO) ? RSqrt(TWO/yddnrm) : RSqrt(hg*hub);
    count++;
    if (count >= MAX_ITERS) break;
    hrat = hnew/hg;
    if ((hrat > HALF) && (hrat < TWO)) break;
    if ((count >= 2) && (hrat > TWO)) {
      hnew = hg;
      break;
    }
    hg = hnew;
  }
  
  /* Apply bounds, bias factor, and attach sign */

  h0 = H_BIAS*hnew;
  if (h0 < hlb) h0 = hlb;
  if (h0 > hub) h0 = hub;
  if (sign == -1) h0 = -h0;
  h = h0;
  return(TRUE);
}

/******************** CVUpperBoundH0 ******************************

 This routine sets an upper bound on abs(h0) based on
 tdist = abs(tout - t0) and the values of y[i]/y'[i].

******************************************************************/

static realtype CVUpperBoundH0(CVodeMem cv_mem, realtype tdist)
{
  realtype atoli, hub_inv, hub;
  booleantype vectorAtol;
  N_Vector temp1, temp2;

  atoli = 0;
  vectorAtol = (itol == SV);
  if (!vectorAtol) atoli = *((realtype *) abstol);
  temp1 = tempv;
  temp2 = acor;
  N_VAbs(zn[0], temp1);
  N_VAbs(zn[1], temp2);
  if (vectorAtol) {
    N_VLinearSum(HUB_FACTOR, temp1, ONE, (N_Vector)abstol, temp1);
  } else {
    N_VScale(HUB_FACTOR, temp1, temp1);
    N_VAddConst(temp1, atoli, temp1);
  }
  N_VDiv(temp2, temp1, temp1);
  hub_inv = N_VMaxNorm(temp1);
  hub = HUB_FACTOR*tdist;
  if (hub*hub_inv > ONE) hub = ONE/hub_inv;
  return(hub);
}

/****************** CVYddNorm *************************************

 This routine computes an estimate of the second derivative of y
 using a difference quotient, and returns its WRMS norm.

******************************************************************/

static realtype CVYddNorm(CVodeMem cv_mem, realtype hg)
{
  realtype yddnrm;
  
  N_VLinearSum(hg, zn[1], ONE, zn[0], y);
  f(N, tn+hg, y, tempv, f_data);
#ifdef DEBUG_CVODE
  if (cvode_error == TRUE) {
	  output_msg(OUTPUT_CVODE,"CVYddNorm error\n");
  }
#endif
  nfe++;
  N_VLinearSum(ONE, tempv, -ONE, zn[1], tempv);
  N_VScale(ONE/hg, tempv, tempv);

  yddnrm = N_VWrmsNorm(tempv, ewt);
  return(yddnrm);
}

/********************* CVStep **************************************
 
 This routine performs one internal cvode step, from tn to tn + h.
 It calls other routines to do all the work.

 The main operations done here are as follows:
  * preliminary adjustments if a new step size was chosen;
  * prediction of the Nordsieck history array zn at tn + h;
  * setting of multistep method coefficients and test quantities;
  * solution of the nonlinear system;
  * testing the local error;
  * updating zn and other state data if successful;
  * resetting stepsize and order for the next step.
  * if SLDET is on, check for stability, reduce order if necessary.
 On a failure in the nonlinear system solution or error test, the
 step may be reattempted, depending on the nature of the failure.

********************************************************************/

static int CVStep(CVodeMem cv_mem)
{
  realtype saved_t, dsm;
  int ncf, nef, nflag;
  booleantype passed;

  int kflag;
  
  saved_t = tn;
  ncf = nef = 0;
  nflag = FIRST_CALL;

  
  if ((nst > 0) && (hprime != h)) CVAdjustParams(cv_mem);
  
  /* Looping point for attempts to take a step */
  loop {  
    cvode_test = TRUE;
    f(N, tn, y, ftemp, f_data);
    cvode_test = FALSE;
    if (cvode_error == TRUE) {
#ifdef DEBUG_CVODE
	    output_msg(OUTPUT_CVODE,"Before predict, y Fail, time %e\n", tn);
#endif
    } else {
	  cvode_prev_good_time = cvode_last_good_time;
	  N_VScale(1.0, cvode_last_good_y, cvode_prev_good_y);
	  cvode_last_good_time = tn;
	  N_VScale(1.0, y, cvode_last_good_y);
#ifdef DEBUG_CVODE
	  output_msg(OUTPUT_CVODE,"Before predict, y OK, time %e\n", tn);
#endif
    }
#ifdef DEBUG_CVODE
    cvode_test = TRUE;
    f(N, tn, zn[0], ftemp, f_data);
    cvode_test = FALSE;
    if (cvode_error == TRUE) {
	    output_msg(OUTPUT_CVODE,"Before predict, zn Fail\n");
    } else {
	    output_msg(OUTPUT_CVODE,"Before predict, zn OK\n");
    }
    saved_t = tn;
#endif
    CVPredict(cv_mem);  
#ifdef DEBUG_CVODE
    cvode_test = TRUE;
    f(N, tn, y, ftemp, f_data);
    cvode_test = FALSE;
    if (cvode_error == TRUE) {
	    output_msg(OUTPUT_CVODE,"After predict, y Fail\n");
    } else {
	    output_msg(OUTPUT_CVODE,"After predict, y OK\n");
    }
    cvode_test = TRUE;
    f(N, tn, zn[0], ftemp, f_data);
    cvode_test = FALSE;
    if (cvode_error == TRUE) {
	    output_msg(OUTPUT_CVODE,"After predict, zn Fail\n");

    } else {
	    output_msg(OUTPUT_CVODE,"After predict, zn OK\n");
    }
#endif
    CVSet(cv_mem);

    nflag = CVnls(cv_mem, nflag);
#ifdef DEBUG_CVODE
    cvode_test = TRUE;
    f(N, tn, y, ftemp, f_data);
    cvode_test = FALSE;
    if (cvode_error == TRUE) {
	    output_msg(OUTPUT_CVODE,"After CVnls, y Fail\n");
    } else {
	    output_msg(OUTPUT_CVODE,"After CVnls, y OK\n");
    }
    cvode_test = TRUE;
    f(N, tn, zn[0], ftemp, f_data);
    cvode_test = FALSE;
    if (cvode_error == TRUE) {
	    output_msg(OUTPUT_CVODE,"After CVnls, zn Fail\n");
    } else {
	    output_msg(OUTPUT_CVODE,"After CVnls, zn OK\n");
    }
#endif
    kflag = CVHandleNFlag(cv_mem, &nflag, saved_t, &ncf);
    if (kflag == PREDICT_AGAIN) continue;
    if (kflag != DO_ERROR_TEST) return(kflag);
    /* Return if nonlinear solve failed and recovery not possible. */
#ifdef DEBUG_CVODE
    cvode_test = TRUE;
    f(N, tn, y, ftemp, f_data);
    cvode_test = FALSE;
    if (cvode_error == TRUE) {
	    output_msg(OUTPUT_CVODE,"Before error test, y Fail\n");
    } else {
	    output_msg(OUTPUT_CVODE,"Before error test, y OK\n");
    }
    cvode_test = TRUE;
    f(N, tn, zn[0], ftemp, f_data);
    cvode_test = FALSE;
    if (cvode_error == TRUE) {
	    output_msg(OUTPUT_CVODE,"Before error test, zn Fail\n");
    } else {
	    output_msg(OUTPUT_CVODE,"Before error test, zn OK\n");
    }
#endif
    passed = CVDoErrorTest(cv_mem, &nflag, &kflag, saved_t, &nef, &dsm);
#ifdef DEBUG_CVODE
    cvode_test = TRUE;
    f(N, tn, y, ftemp, f_data);
    cvode_test = FALSE;
    if (cvode_error == TRUE) {
	    output_msg(OUTPUT_CVODE,"After error test, y Fail, passed %d\n", passed);
    } else {
	    output_msg(OUTPUT_CVODE,"After error test, y OK, passed %d\n", passed);
    }
    cvode_test = TRUE;
    f(N, tn, zn[0], ftemp, f_data);
    cvode_test = FALSE;
    if (cvode_error == TRUE) {
	    output_msg(OUTPUT_CVODE,"After error test, zn Fail\n");
    } else {
	    output_msg(OUTPUT_CVODE,"After error test, zn OK\n");
    }
#endif
    /* Return if error test failed and recovery not possible. */
    if ((!passed) && (kflag == REP_ERR_FAIL)) return(kflag);
    if (passed) break;
    /* Retry step if error test failed, nflag == PREV_ERR_FAIL */
  }
#ifdef DEBUG_CVODE
  output_msg(OUTPUT_CVODE, "Finished step in CVStep\n");
#endif
  /* Nonlinear system solve and error test were both successful.
     Update data, and consider change of step and/or order.       */


  CVCompleteStep(cv_mem); 
#ifdef DEBUG_CVODE
    cvode_test = TRUE;
    f(N, tn, y, ftemp, f_data);
    cvode_test = FALSE;
    if (cvode_error == TRUE) {
	    output_msg(OUTPUT_CVODE,"After complete step, y Fail\n");
    } else {
	    output_msg(OUTPUT_CVODE,"After complete step, y OK\n");
    }
    cvode_test = TRUE;
    f(N, tn, zn[0], ftemp, f_data);
    cvode_test = FALSE;
    if (cvode_error == TRUE) {
	    output_msg(OUTPUT_CVODE,"After complete step, zn Fail\n");
    } else {
	    output_msg(OUTPUT_CVODE,"After complete step, zn OK\n");
    }
#endif
  CVPrepareNextStep(cv_mem, dsm); 

  /* If Stablilty Limit Detection is turned on, call stability limit
     detection routine for possible order reduction. */

  if (sldeton) CVBDFStab(cv_mem);
#ifdef DEBUG_CVODE
    cvode_test = TRUE;
    f(N, tn, y, ftemp, f_data);
    cvode_test = FALSE;
    if (cvode_error == TRUE) {
	    output_msg(OUTPUT_CVODE,"After cvbfdstab, y Fail\n");
    } else {
	    output_msg(OUTPUT_CVODE,"After cvbfdstab, y OK\n");
    }
    cvode_test = TRUE;
    f(N, tn, zn[0], ftemp, f_data);
    cvode_test = FALSE;
    if (cvode_error == TRUE) {
	    output_msg(OUTPUT_CVODE,"After cvbfdstab, zn Fail\n");
    } else {
	    output_msg(OUTPUT_CVODE,"After cvbfdstab, zn OK\n");
    }
#endif
  etamax = (nst <= SMALL_NST) ? ETAMX2 : ETAMX3;

  /*  Finally, we rescale the acor array to be the 
      estimated local error vector. */

  N_VScale(ONE/tq[2], acor, acor);
  return(SUCCESS_STEP);
      
}


/********************* CVAdjustParams ********************************

 This routine is called when a change in step size was decided upon,
 and it handles the required adjustments to the history array zn.
 If there is to be a change in order, we call CVAdjustOrder and reset
 q, L = q+1, and qwait.  Then in any case, we call CVRescale, which
 resets h and rescales the Nordsieck array.

**********************************************************************/

static void CVAdjustParams(CVodeMem cv_mem)
{
  if (qprime != q) {
    CVAdjustOrder(cv_mem, qprime-q);
    q = qprime;
    L = q+1;
    qwait = L;
  }
  CVRescale(cv_mem);
}

/********************* CVAdjustOrder *****************************

  This routine is a high level routine which handles an order
  change by an amount deltaq (= +1 or -1). If a decrease in order
  is requested and q==2, then the routine returns immediately.
  Otherwise CVAdjustAdams or CVAdjustBDF is called to handle the
  order change (depending on the value of lmm).

******************************************************************/

static void CVAdjustOrder(CVodeMem cv_mem, int deltaq)
{
  if ((q==2) && (deltaq != 1)) return;
  
  switch(lmm){
    case ADAMS: CVAdjustAdams(cv_mem, deltaq);
                break;
    case BDF:   CVAdjustBDF(cv_mem, deltaq);
                break;
  }
}

/*************** CVAdjustAdams ***********************************

 This routine adjusts the history array on a change of order q by
 deltaq, in the case that lmm == ADAMS.

*****************************************************************/

static void CVAdjustAdams(CVodeMem cv_mem, int deltaq)
{
  int i, j;
  realtype xi, hsum;

  /* On an order increase, set new column of zn to zero and return */
  
  if (deltaq==1) {
    N_VConst(ZERO, zn[L]);
    return;
  }

  /* On an order decrease, each zn[j] is adjusted by a multiple
     of zn[q].  The coefficients in the adjustment are the 
     coefficients of the polynomial x*x*(x+xi_1)*...*(x+xi_j),
     integrated, where xi_j = [t_n - t_(n-j)]/h.               */

  for (i=0; i <= qmax; i++) l[i] = ZERO;
  l[1] = ONE;
  hsum = ZERO;
  for (j=1; j <= q-2; j++) {
    hsum += tau[j];
    xi = hsum / hscale;
    for (i=j+1; i >= 1; i--) l[i] = l[i]*xi + l[i-1];
  }
  
  for (j=1; j <= q-2; j++) l[j+1] = q * (l[j] / (j+1));
  
  for (j=2; j < q; j++)
    N_VLinearSum(-l[j], zn[q], ONE, zn[j], zn[j]);
}

/********************** CVAdjustBDF *******************************

 This is a high level routine which handles adjustments to the
 history array on a change of order by deltaq in the case that 
 lmm == BDF.  CVAdjustBDF calls CVIncreaseBDF if deltaq = +1 and 
 CVDecreaseBDF if deltaq = -1 to do the actual work.

******************************************************************/

static void CVAdjustBDF(CVodeMem cv_mem, int deltaq)
{
  switch(deltaq) {
    case 1 : CVIncreaseBDF(cv_mem);
             return;
    case -1: CVDecreaseBDF(cv_mem);
             return;
  }
}

/******************** CVIncreaseBDF **********************************

 This routine adjusts the history array on an increase in the 
 order q in the case that lmm == BDF.  
 A new column zn[q+1] is set equal to a multiple of the saved 
 vector (= acor) in zn[qmax].  Then each zn[j] is adjusted by
 a multiple of zn[q+1].  The coefficients in the adjustment are the 
 coefficients of the polynomial x*x*(x+xi_1)*...*(x+xi_(q-1)),
 where xi_j = [t_n - t_(n-j)]/h.

*********************************************************************/

static void CVIncreaseBDF(CVodeMem cv_mem)
{
  realtype alpha0, alpha1, prod, xi, xiold, hsum, A1;
  int i, j;
  
  for (i=0; i <= qmax; i++) l[i] = ZERO;
  l[2] = alpha1 = prod = xiold = ONE;
  alpha0 = -ONE;
  hsum = hscale;
  if (q > 1) {
    for (j=1; j < q; j++) {
      hsum += tau[j+1];
      xi = hsum / hscale;
      prod *= xi;
      alpha0 -= ONE / (j+1);
      alpha1 += ONE / xi;
      for (i=j+2; i >= 2; i--) l[i] = l[i]*xiold + l[i-1];
      xiold = xi;
    }
  }
  A1 = (-alpha0 - alpha1) / prod;
  N_VScale(A1, zn[qmax], zn[L]);
  for (j=2; j <= q; j++) {
    N_VLinearSum(l[j], zn[L], ONE, zn[j], zn[j]);
  }  
}

/********************* CVDecreaseBDF ******************************

 This routine adjusts the history array on a decrease in the 
 order q in the case that lmm == BDF.  
 Each zn[j] is adjusted by a multiple of zn[q].  The coefficients
 in the adjustment are the coefficients of the polynomial
 x*x*(x+xi_1)*...*(x+xi_(q-2)), where xi_j = [t_n - t_(n-j)]/h.

******************************************************************/

static void CVDecreaseBDF(CVodeMem cv_mem)
{
  realtype hsum, xi;
  int i, j;
  
  for (i=0; i <= qmax; i++) l[i] = ZERO;
  l[2] = ONE;
  hsum = ZERO;
  for(j=1; j <= q-2; j++) {
    hsum += tau[j];
    xi = hsum /hscale;
    for (i=j+2; i >= 2; i--) l[i] = l[i]*xi + l[i-1];
  }
  
  for(j=2; j < q; j++)
    N_VLinearSum(-l[j], zn[q], ONE, zn[j], zn[j]);
}

/**************** CVRescale ***********************************

  This routine rescales the Nordsieck array by multiplying the
  jth column zn[j] by eta^j, j = 1, ..., q.  Then the value of
  h is rescaled by eta, and hscale is reset to h.

***************************************************************/

static void CVRescale(CVodeMem cv_mem)
{
  int j;
  realtype factor;
  
  factor = eta;
  for (j=1; j <= q; j++) {
    N_VScale(factor, zn[j], zn[j]);
    factor *= eta;
  }
  h = hscale * eta;
  hscale = h;
  nscon = 0;
}

/********************* CVPredict *************************************

 This routine advances tn by the tentative step size h, and computes
 the predicted array z_n(0), which is overwritten on zn.  The
 prediction of zn is done by repeated additions.

*********************************************************************/

static void CVPredict(CVodeMem cv_mem)
{
  int j, k;
  
  tn += h;
  for (k = 1; k <= q; k++)
    for (j = q; j >= k; j--) 
      N_VLinearSum(ONE, zn[j-1], ONE, zn[j], zn[j-1]); 
}

/************************** CVSet *********************************

 This routine is a high level routine which calls CVSetAdams or
 CVSetBDF to set the polynomial l, the test quantity array tq, 
 and the related variables  rl1, gamma, and gamrat.

******************************************************************/

static void CVSet(CVodeMem cv_mem)
{
  switch(lmm) {
    case ADAMS: CVSetAdams(cv_mem);
                break;
    case BDF  : CVSetBDF(cv_mem);
                break;
  }
  rl1 = ONE / l[1];
  gamma = h * rl1;
  if (nst == 0) gammap = gamma;
  gamrat = (nst > 0) ? gamma / gammap : ONE;  /* protect x / x != 1.0 */
}

/******************** CVSetAdams *********************************

 This routine handles the computation of l and tq for the
 case lmm == ADAMS.

 The components of the array l are the coefficients of a
 polynomial Lambda(x) = l_0 + l_1 x + ... + l_q x^q, given by
                          q-1
 (d/dx) Lambda(x) = c * PRODUCT (1 + x / xi_i) , where
                          i=1
 Lambda(-1) = 0, Lambda(0) = 1, and c is a normalization factor.
 Here xi_i = [t_n - t_(n-i)] / h.

 The array tq is set to test quantities used in the convergence
 test, the error test, and the selection of h at a new order.

*****************************************************************/

static void CVSetAdams(CVodeMem cv_mem)
{
  realtype m[L_MAX], M[3], hsum;
  
  if (q == 1) {
    l[0] = l[1] = tq[1] = tq[5] = ONE;
    tq[2] = TWO;
    tq[3] = TWELVE;
    tq[4] = CORTES * tq[2];       /* = 0.1 * tq[2] */
    return;
  }
  
  hsum = CVAdamsStart(cv_mem, m);
  
  M[0] = CVAltSum(q-1, m, 1);
  M[1] = CVAltSum(q-1, m, 2);
  
  CVAdamsFinish(cv_mem, m, M, hsum);
}

/****************** CVAdamsStart ********************************

 This routine generates in m[] the coefficients of the product
 polynomial needed for the Adams l and tq coefficients for q > 1.
  
******************************************************************/

static realtype CVAdamsStart(CVodeMem cv_mem, realtype m[])
{
  realtype hsum, xi_inv, sum;
  int i, j;
  
  hsum = h;
  m[0] = ONE;
  for (i=1; i <= q; i++) m[i] = ZERO;
  for (j=1; j < q; j++) {
    if ((j==q-1) && (qwait == 1)) {
      sum = CVAltSum(q-2, m, 2);
      tq[1] = m[q-2] / (q * sum);
    }
    xi_inv = h / hsum;
    for (i=j; i >= 1; i--) m[i] += m[i-1] * xi_inv;
    hsum += tau[j];
    /* The m[i] are coefficients of product(1 to j) (1 + x/xi_i) */
  }
  return(hsum);
}

/****************** CVAdamsFinish  *******************************

 This routine completes the calculation of the Adams l and tq.

******************************************************************/

static void CVAdamsFinish(CVodeMem cv_mem, realtype m[], realtype M[], 
                          realtype hsum)
{
  int i;
  realtype M0_inv, xi, xi_inv;
  
  M0_inv = ONE / M[0];
  
  l[0] = ONE;
  for (i=1; i <= q; i++) l[i] = M0_inv * (m[i-1] / i);
  xi = hsum / h;
  xi_inv = ONE / xi;
  
  tq[2] = xi * M[0] / M[1];
  tq[5] = xi / l[q];

  if (qwait == 1) {
    for (i=q; i >= 1; i--) m[i] += m[i-1] * xi_inv;
    M[2] = CVAltSum(q, m, 2);
    tq[3] = L * M[0] / M[2];
  }

  tq[4] = CORTES * tq[2];
}

/****************** CVAltSum **************************************
  
 CVAltSum returns the value of the alternating sum
   sum (i= 0 ... iend) [ (-1)^i * (a[i] / (i + k)) ].
 If iend < 0 then CVAltSum returns 0.
 This operation is needed to compute the integral, from -1 to 0,
 of a polynomial x^(k-1) M(x) given the coefficients of M(x).

******************************************************************/

static realtype CVAltSum(int iend, realtype a[], int k)
{
  int i, sign;
  realtype sum;
  
  if (iend < 0) return(ZERO);
  
  sum = ZERO;
  sign = 1;
  for (i=0; i <= iend; i++) {
    sum += sign * (a[i] / (i+k));
    sign = -sign;
  }
  return(sum);
}

/***************** CVSetBDF **************************************

 This routine computes the coefficients l and tq in the case
 lmm == BDF.  CVSetBDF calls CVSetTqBDF to set the test
 quantity array tq. 

 The components of the array l are the coefficients of a
 polynomial Lambda(x) = l_0 + l_1 x + ... + l_q x^q, given by
                                 q-1
 Lambda(x) = (1 + x / xi*_q) * PRODUCT (1 + x / xi_i) , where
                                 i=1
 xi_i = [t_n - t_(n-i)] / h.

 The array tq is set to test quantities used in the convergence
 test, the error test, and the selection of h at a new order.


*****************************************************************/

static void CVSetBDF(CVodeMem cv_mem)
{
  realtype alpha0, alpha0_hat, xi_inv, xistar_inv, hsum;
  int i,j;
  
  l[0] = l[1] = xi_inv = xistar_inv = ONE;
  for (i=2; i <= q; i++) l[i] = ZERO;
  alpha0 = alpha0_hat = -ONE;
  hsum = h;
  if (q > 1) {
    for (j=2; j < q; j++) {
      hsum += tau[j-1];
      xi_inv = h / hsum;
      alpha0 -= ONE / j;
      for(i=j; i >= 1; i--) l[i] += l[i-1]*xi_inv;
      /* The l[i] are coefficients of product(1 to j) (1 + x/xi_i) */
    }
    
    /* j = q */
    alpha0 -= ONE / q;
    xistar_inv = -l[1] - alpha0;
    hsum += tau[q-1];
    xi_inv = h / hsum;
    alpha0_hat = -l[1] - xi_inv;
    for (i=q; i >= 1; i--) l[i] += l[i-1]*xistar_inv;
  }

  CVSetTqBDF(cv_mem, hsum, alpha0, alpha0_hat, xi_inv, xistar_inv);
}

/****************** CVSetTqBDF ************************************

 This routine sets the test quantity array tq when lmm == BDF.

******************************************************************/

static void CVSetTqBDF(CVodeMem cv_mem, realtype hsum, realtype alpha0,
                      realtype alpha0_hat, realtype xi_inv, realtype xistar_inv)
{
  realtype A1, A2, A3, A4, A5, A6;
  realtype C, CPrime, CPrimePrime;
  
  A1 = ONE - alpha0_hat + alpha0;
  A2 = ONE )set t    xi = hs"Be hsum, alpha0, alpha0_hatEi=1; i, s in    hsum p)NoecreaseBDle(A( == ] = CVAlt/altype A1, l) {
    for (i=q; i >= 1; CE;
  for (i=2;it == 1) 1; A3NE;
  hsu"Ber_inv = -l[1]A4NE;
  hsum = h+i=q; i >= 1; A1 = Ooecr3e(A( alphaA4N+cr3 M[0] / M[1NoecreasA1 = Oo/ C)q-1];
    xi_inv = / hsum;
    alpha0_hat = -l[1]A5NE;
  hsu"-A( alp
   2; j < q1]A6NE;
  hsum = h-i=q; i >= 1; A1 = O1 = Ooecr2e(A( alphaA6N+cr5 M[0] / M[2];
 reasA1 = O1 = Oo = CVAlt1; iq+2 / xr5 M[0]}* tq[2];
}

/****************** CVAltSum ***********nls***************************

 This routinevances tn by the step */
  lst werem solution or error assocg pol
f the a quot
***ati(1 ET is e the cotion oroefficients and is Dlue of lmm). pro,quiro set thnlstation * liw;
    }Newtut - rk.

***********************************************************/

static void CVAdjustBDF(Cm cv_m;
#int deltaq)
{
  switchUG_CVO{
  case SS: rer: CVSetAdams(FUNCTIONAL V(cv_mem, ynlstation * l        )CVSetBDF(cvem);
 a0_haV(cv_mem, ynlsNewtutifdef DEBUG_CVO9);
}

/*********************** CVEwtSetSS  thnlstation * li*****************

 This routinevances tn by the step */
  lst werem solution or error  quoti
d by CVod li probably (no mait wficievt wed)***************************************************/

static void CVAdjustBDF(Cm cv_m;
#tation * l 
  int j, k;
  
  tn += hmrime, CPrimePdel     p    thr increaIise, CVHicount  = q+1,enatiiderf
***n(0), whicyq==1) {
    pRST_CALLcr   *t  if (q mRST_CALLemp, f_data);
  fdef DEBUG_CVODEarSum(ONE,  TRUE) {
#ifdef DEBUG_CVODE
	    output_msg(OUTPUCVODE,"CVYddNorm error\n");;
#tation * l,lse {
***beginnuotNonlinear syste/********CONVag);
 _VScale(HUB_F
	    output_msg(OUTPUCVODE,"CVYddNorm error\n");;
#tation * l,lOK
***beginnuotNonlinear systetemp1, ;
    returnCCESS_SAX_ITERS timn ll, the error ;typeu q 
ticourreion ***ll Ce 
 d + 5)**sign;
   s fourreiocyqer.

****
 It ct is = PRrnatingNE / q;n[1], ONE, zn[1], tempv);
  N_VScale(ONE/hg, ttempv, tempvl1 N_VAddConst(tempv,, zn[j-1], ONE, zn[j], 0j-1]); 
_VAddCoyf nonlineaGies ydd / 2)) = 
 zn[j] ourreion *
  luq in t the error ay tq NE / q;n[1], ONE, znzn[1], tempv);
  Nurn(SUCCESS_STE
    pv, ewt);
  return(SUacor);
    recor, acor1], tempvCCESS_STE
 re coeffilose to  the error CVSIf mR> ric mecond derivative  the error teha0_hat
ticou   anset ton thce
 tcr   n of hgence
 test,y tqnt = 0;
  lo& (qwaitmmmap :cr   *t out)CRDOWNtionr   n    p/    pS_STE
  ******   p* MINacor1]nr   le(Aq[2];ror == TRU ****<=b_inv;u[j-1];
actempv,(mmma;
  ?*   p:, ewt);
  return(SUacor);
   TEP);
    OLVED);  s fouhe error aachieor) /* j = q re coeff2 aftaj;
  e 
  probablys cal TR pro.as amn order,die eruoti/* j = m>= MAX_ITERS(m==
  e 
 outputmmmO)) {
      p> RDIVp*    pS));
   TEP);
   CONVag);
 _VSccoeff2ave / 2)) = 
urreion *,,enatiiderfn of h**sigage he/* j =    pRST   ftemp, f_data);
  fdef DEBUG_CVODEar== TRUE) {
#ifdef DEBUG_CVODE
	    output_msg(OUTPUT_CVODE,"Before predict, y );;
#tation * l,lse {
***ar continue;
    }

 P);
   CONVag);
 _VSccoale(HUB_F
	    output_msg(OUTPUPUCVODE,"CVYddNorm error\n");;
#tation * l,lOK
***ar continue;
   j = q */
 Sum(ONE, }********** CVEwtSetSV *********nlsNewtut*****************

 This routine adjusts the historation of l aNewtut* probably.utines to (HUtuiled,t ar),  [t_n o set thNewtutIprobably ordnal cvom -1 tprobably me if drxit tial y not postep *
***Newtut* probablyss eacaset tar),  [t_ludeeely.
  O[i].

**tingo e the locVnot*******************************************************/

static void CVAdjustParamm cv_m;
#Newtutint deltaq)
{
  switchUG_CVO{
  cmp2;

  avub_inv v1);
  hv1);
3l[0] = l thef passierssed;

  int kflo seSUtuialpha0_vemp2 = a0], tcoeffn[jam the 
  s_vemp2 = to , Cd detectie/* j vVAbs(zn[y;re coeffn[jam ty  s_vemp22= to , Cd detectie ie/* j vVAbs3= acor;
 oeffn[jam t   N_V s_vemp23= to , Cd detecti==1) {
 pper boAIL * thef passiord. ****HUtuilf 
  ps,enatiidblys it s/****/fined ef papv,((FAIL */
 if ((nst > outpuFAIL */
  }
#ifdef DEB)  ?;
    Of DEBUR***:  DEB_OTHER_SAX_ITEDit ha
*** succmpone. ***it
   Utuile histor(s editie  fos  /* j Stab(UtuiNonNulSum(*******SetBDFseSUtuima / AIL */
  }
#iCONVag);
 outpuFAIL */
  }
#ifdef DEB)out;
   TE= gamma;
  utpuFs TWO)Fs  pR+ MSBP  utpureas ? gam-_inv;> DGout _VScale(HUB_***SetBDr   *t  if (q tBDFseSUtuima de_error }oint for attempts to take a at alinear system;
  * testing the , and cEnatiiderf
***y z_n(0), whicyimit
  (HUtuiled,ar),  [t_n nitiaq tBDFse thNewtutIprobably ke a at Newtut* probablysstself.= 0;
  lo&  5)**sign;], ftemp, f_data);
    cvode_test = F */
 Sum(ON Ear== TRUE) {
#ifdef DEBUG_CVODE
	    output_msg(OUTPUT_CVODE,"Before predict, y );;
#Newtut,    rthe
 c*si
#endif
    }
#ifdef DEBUG }

 P);
    CONVag);
 _VSccoale(HUB_F
	    output_msg(OUTPUPUCVODE,"CVYddNorm error\n");;
#Newtut, pass   rthe
 c*si
#endif
    }
#ifdef DEBUGj = q */
    al TRUEFseSUtuiv;u[j-1];
ie(cv (HUtuiifdef DEB thef passta);
    cvode&jLinea     realtype alpha0vub_inv v1);
  hv1);
3r);
   TEnHUtuis>= MAX_ItBDFseSUtuima de_error      ? gamma Dr   *t  if  or      ? rat = (nst > 0)  TEnHt pRSTninv > rror test failed *HUtuilf not p/* j = turn(ZERr;
  
  sum = SETUPag);
_UNRECr);
   TErn(ZERr;>  
  sum = CONVag);
 _VSccoa re coeff2 poscto tn/
  
  if loadis done by ro toty  zn[qmaNE / q;n[1;
    returnCCESS_S
    recor, acor1]tn+hg, y, tre coeffDoa at Newtut* probablys lo& (qwe(cv thNewtutIprobablyag = CVnls(cv_m
    f(N, tn, y, ftemp, f_data);
    cvode_test = FALSE;
    if (cvode_error == TRUE) {
#ifdef DEBUG_CVODE
	    output_msg(OUTPUT_CVODE,"Before predict, y \n");
  NewtutIprobably,", passwe(c  }
   we(tinue;
    }

 P);
   CONVag);
 _VSccoale(HUB_F
	    output_msg(OUTPUPUCVODE,"CVYddNorm error\n"\n");
  NewtutIprobably,"passwe(c  }
   we(tinue;
   j = q */
 mit Dehange in at the error alinear sof h at JCCEbian- rl1, ga     rea;
  *appstisone. ***er,
 zn[j],h**sigage he the a DFse ****HUtui
 M_Env maAdams orthef pa=g);
_BAD_J.stAdams or CP);
  .  realtype alpha0 lo& (qwaitwe(c!EBUGY    if (P);
   we(tin_ItBDFseSUtuima  y, ftemp,ed ef papv,g);
_BAD_JNE, }********** CVEwtSetSV ***** thNewtutIprobably **********

 This routine adjusts the histornal cvode l aNewtut* probably.utfm -1 tprobablyUpdateeds,
quirof the alternating OLVED.utfmne.,quir, de(sumal allocanlsNewtut*
he histor***it
  *HUtuilage he if (dostep *
 -1 tprobably m    of ther;
  * unatingUGY    if. (, ewte theo the nlsNewtut*mu (cties
,ed ef pap***g);
_BAD_J CVEwtSeit
  0 ..Utuilage h) cloets h to thtes the Nordsie the aditi the coepproprg po[i].

**
g OLVEag);
_UNRECiw;
 ONVag);
 backp***e nlsNewtut******************************************************/

static void CVPredict(CVom cv_mNewtutIprobablya
  int j, k;
  
  tn += hminv)trime, CPrimePdel     p    thr  cmp2;

  abalpha0_(q mnmachEnmRST_CALL   pRST_CAt for attempts to take aNewtut* probablys lo& **sign;], ftor Enatiiderles  rsidof tystem;
  * testing the NE / q;n[1], ONE, znvl1 N y);
  f(N,  */
  
  ewempv,, zn[j-1], ONE,

*****  cvode-cor1], tempv= N_VWrmsNo  s foe main olecoveryation *0 lo& (qb= acor;
  N_
 P);cv (Hcoveifdef DEBb);
    cvo)  or   nni>= MAX_Io& (qwaitP);c
  
  sum = SOLVEag);
_UNRECS_STE
 re coeffIfolecoverhadiaible. */ae errinear sof hJCCEbianl;
  * s
 M_Env ne. 
 zn[j],h(sumal aoa rya at alinear sage hetype alpha0 lo& (qwaitP);c>  
 { ;
   TErn(Z(!jLin {
   (UtuiNonNulSu)/*********Y    if ;;
   TEP);
   CONVag);
 _VScco}nonlineaGies ydd / 2)) = 
urreion *;**** ourreion *
  lhe 
  nicyq==1) 
    pv, ewt);
  retbSUacor);
    re], ONE, znzn[1] */
  zn[1]brnCCESS_S
    rej-1], ONE, zn[j], 0j-1]); 
 */
  yS_STE
 re coeffilose to  the error CVSIf mR> ric mecond derivative  the error teha0_hat
ticou   anset ton thce
 tcr   n of hgence
 test,y tqnt = 0;
  lo& (qwaitmmmap :u[j-1];
cr   *t out)CRDOWNtionr   n    /   pS_STE
 }STE
  ******   p* MINacor1]nr   le(Aq[2];ror ==E
    cvode_test = TRUE;eff***ale VEwtSe OLVED  lo& (q
    f(N, tn, y, ftemp, f_data);
    cvode_test = FALSE;
    if (cvode_error == TRUE) {
	    output_msg(OUTPUT_CVODE,"After cvbfdstab, y);;
#Newtut, se {
*** OLVEDcontin }

 P);
   CONVag);
 _VSccoale(HUB_FTPUPUCVODE,"CVYddNorm error\n");;
#Newtut, pa
*** OLVEDcontin* Return if error TRU ****<=b_inv;u[j-1];
actempv,(m==
  ?*   p:, ewt);
  return(SUacor);
   TEjLinma de_error     P);
    OLVED); m solve and error te andecovedto a positive /* j = q */
    almnmachEn++m_STE
 re coeff2 aftaj;
  e 
  probablys cal TR pro.as amn order,die eruot.teha0_haIr.  i mane. 
the errdsof hJCCEbianl;
  * s ne. 
 zn[j],hteha0_ha(sumal aoa rya at alinear sage hetype alpha0 realtype alpha0 lo& (qwait(mmma;
  e 
 outputmmmO)) {
      p> RDIV*   pS));u[j-1];
in(Z(!jLin {
   (UtuiNonNulSu)/*********Y    if ;;
   TEP);
   CONVag);
 _VScco}nonliVSccoeff2ave / 2)) = 
urreion *,,enatiiderfn of h**sigage he/* j =    pRST   ftemp, f_data);
    cvode_test = F */
  TRUE) {
#ifdef DEBUG_CVODE
	    output_msg(OUTPUT_CVODE,"Before predict, y );;
#Newtut, se {
***ar continue;
    }

 P);
   CONVag);
 _VSccoale(HUB_F
	    output_msg(OUTPUPUCVODE,"CVYddNorm error\n");;
#Newtut, pa
***ar continue;
   j = q */
 Sum(ONE, }********** CVEwtSetSV *****mem, &nflag, s**************

 This routine completes the calcoop s****ar syrns imm.
  O[i].

mem, nfla*em, nPtr  of theale is    }
#eliminary adj CVAl**nls*pdateednce
 tecovr;
  * u  * testing the  retur
mem, &nflag, s*of the alterou   anseurn(kflag);
 tten on teset thetep
 ordnal cvom -1 e selectios seIstem;
  * testing the e andne. ecovedto a positive returnncfn nitiaflag =*flaPtrof prordericierdsof h multiplying theated a rsn thcs seIstem;
alinear system;
  * testing the lf not pduor***anthe ble. */ae errinear sby..Utuicor arr
  O[lternating ETUPag);
ED. seIstitlf not pduor***anhe ble. */ae errinear s
 tecove returnr arr
  O
alternating OLVEag);
ED. seets h to thaible. */ae errinear soc
 zn[ step secovr;
  * u

  * testing the l, ynls of thealeFAIL */
 CONVag);
 .iVScc, ewte theo thr arr
  O[lternatingurn( ONVag);
  failag s newVScce of the MXNCF cal|h|lphamin.iVScc,fmne.,qwantita*em, nPtr 
  }
#iCONVag);
  if (deltaqlternatinVScctinue;
    if, teser;
 *******he (dostep *
 -1  is on****************************************************/

static void CVPredict(CVom cv_m, &nflag, sant deltaq)
{
  switch*em, nPtr, dsm;
  int ncf, n                i=1
 xi_i itch*elaPtr
  tn += hn = tn;
  ncem, nfla*em, nPtr0) return(ZFAIL */
  OLVED) P);
   urn(kflag);
  out too cm;
  * testingoln.if not ;rordericienncfn nita rsn theate/* j ncfn(ONE,  cv_mn thifdef DEBt ncf,  point for st failed *HUtuilcallecovery not pe ble. */ae y /* j StabFAIL */
  ETUPag);
_UNRECr  sum = SETUPag);
ED);eturn(ZFAIL */
  OLVEag);
_UNRECr  sum = SOLVEag);
ED point for Atns
    to t,eFAIL */
 CONVag);
;rordericienncf  lo&  5)(*elaPtr
(ONE, MALL_NST) if (q effIfowerhadiMXNCF rinear s cal|h|lphamin, (deltaqurn( ONVag);
 /* j Stabureash)*<=bamin* ifPSM outpu*flaPtro/
 MXNCF));
    sum = urn( ONVag);
 out too Ressaryomputes
 ; (deltaqle (dostep *
 -1  is  /* j MAL*t out)ETACF,bamins ineash) poin*em, nPtr 
  }
#iCONVag);
;m);
}

/**************   sum = tinue;
    if ********** CVAdjustBDF ********v_mn th*******************

 This routine sets the test quan rsn th alternating sumaqle if
    Cof hgndoh alte is done by nt zn OK\nxecnear syst**v_mn th,ay by multiplying theateharray z_sam ti].

**tse VEwtSe the Adlp***e v_mem)
*****************************************************/

static int CVStep(CVodeMedeMem cv_mn thiealtype hsum, realtype alpt ncf,    tn += h;
  for (k = 1=pt ncf, <= q; k++)
    for (j = q; j >= k; j--) 
      N_VLinearSm(ONE, zn[j-1], ONE, zn[j], zn[j-1-]); 
}

/********************* CVAdjustBDF **mem, &nflag, &i*****************

 This routinevances tn by thenal cvode l a. */

  N_VSy tqntset toweighctor. */

  N_VS/ 2))dsmg s loadnce
 ****dsmPtr, e if darray tq dsmg? (j1g s madns seIstem;
y tq  /* Rs,mem, &nflag, &iof the aUG_Cs of Istem;
y tq rines,qwangndo
 -1  is  oeffcoop nRUEFsem cv_mn th),htetita*em, nPtr *** }
#ifdef DEBWRMS norm.

 de_ers of IstMXNEF, nflag == PREV_ar s have oc
 zn[ scal TRreash)*phamin,
qwantita*km, nPtr 
 urn(kflag);
. (Adams or C*km, nPtr haralte inatings = Pof theale is  , &nflam, n.)of IstmwtSe tanhMXNEF1, nflag == PREV_ar s have oc
 zn[ , by an amouf (sldetog s  cvchcs s************************************************/

static int CVStep(CVodeMe;

  int kflem, &nflag, &knt deltaq)
{
  switch*em, nPtr, itch*km, nPtr,                i=1
 xi_i = [t_n dsm;
  int ncf, neitch*eefPtr, dsm;
  in*dsmPtr A2, A3, A4, A5f, nflaALL sm= a0]temp/********(q effIfo tqnt. */

  N_VS/ 2))dsmg /* Rsection orm.

 UG_C /*laALL*dsmPtrRST , n {
    N_smg<=b_inv;**************int for T== PREV_ER;rordericienount  =sf zn tof, &dsnita rsn theateng the/* j (*eefPtr
(ONE, netf(ONE, *em, nPtr 
  }
#ikflag);
NE,  cv_mn thifdef DEBt ncf,  pot for AtnMXNEF,rinear s cal|h|lphamin, (deltaq the lag(cv_murn(kflag);
 /* j Stabureash)*<=bamin* ifPSM outpu*feaPtro/
 MXNEF));u[j-1]*km, nPtr 
 urn(kflag);

  /* On an (de_er)order decrea2 poMALL_NST)i) */vodecienomputes
 torder q i***ar e the loc is  /* j MALL_NST) if (decrea2 poh obabl MAL*
 It _sm,and hscaWRMS norm.

  to ,  ryaer.     /* j Stab*feaPtro<
 MXNEF i >= 1; MAL*t  alp
  RPowerR(BIAS2*dsm,cor)L)N+crDDOf ;;
   MAL*t out)ETAMif, out) reseamins ineash) ODEar== TRU*feaPtro>MX2 : ETAEF) MAL*t oIN) reseally,FODEar==
}

/**************  * On an (de_er)order int for An OK\MXNEF1,rinear s,  cvch, each zn[f (sldetogMS normest fail/* j Stabj=2; j < q; jMAL*t out)ETAMif, amins ineash) poin_mem, qprime-q);
    q -= ONE /  = L; }
  CV-- }
  CVRescale(cv_m==
}

/**************  * On an (de_er)order (q effIfoal, Cdyi***ch zn[1,and   rtV(cvloadiate
 It snr  ch /* j MAL*t out)ETAMif, amins ineash) poinh hscale * en = 0;
}

/***VRescale(ONG_WAIT/**************LLemp, f_data);
  fdef DEBUG_CVODEar TRUE) {
	  output_msg(OUTPUT_CVODE,"CVYddNorm error\n "em, &nflag, &ntin }
/*exit(8);/*  }
  out"CVYd"em, &nflag, &n,)i)rea2TOP /*)order inSum(ONE, temcor, a[1], tempv y);
***   sum = de_er)or* CVAdjustAdams ******m); 
#ifdef *****************

 This routine adjusts the histornal cvodemma,ousther seere are as fw h is realinear 
ray zn.
  * testing the lharaed);
   l a. */

  N_VSy tqntseWprordericie
 -1  is  ount  = nion orcch alternatins hue routu,  her seeem;
yaueng thWRMS napplyalterourreion ***ay zn.
ateng th.set tonv =s of prct is = Pqrnatins r.

,q the nv =1 of thmo= Pofccie.f the arnt  = VRescais ordericierdWRMS n TRfor (i=q;  ( rout hsuL_N)
qwantave he 
  nicm p)Noke a duction. */

  iorder q ***************************************************/

static void CVSetTqBDF(CVodeMem m); 
#ifdef D
  int j, k;
  
  tn += hh;
  m[0] =nio(ONE, n****(ONE, hue}

/***Vu= L; }
; i--) l[i] += l[]*xi + lonv =s o_inv =(j=2; j taq != 1) {
    gammaCVAlnv =2 o_inv = l[q];
v =1 o}

/*; j++) {
 0   N_VScale(facl[j], zn[L], ONE, zn[j], */
  zn[1] factor *= eta;
 for (-- }
 taq !=or (i=q; i 
    ustOrL_N)actor, zn[j], zn[]); 
 */
    for (jn;
  }
  ncf, q5  tq[2] NE, }********** CVEwtSe(cv_mem, dsm); 

 *****************

 This routine adjusts the historation of l azn = 0 .er.    es
 t nic/

  i== ADAMS.nm);  is  -- justPar routustPant zlo0 . the justPa,quirtity arr
 obabl MAL*= justPa*/

 I Colso her ses r succqBDFeamma, and g
  rl1, gahe sa/or order.     es
 t


  CVCom**************************************************/

static void CVSetT ict(CVodeMem cv_mm, dsm); 

  ealtype hsum, realtype alpabli  tn effIfo ALL_NST)in   fuccqB   es
 t


  CVCa/or orsl/* j Stab ALL_NST=b_inv;u[j-1]VRescaleout)VRescL * M[0] / CVAdju L; }
  ChCVAdju Lh;= 1; MAL*t  al
  /* On an order decreaeat 
 ustmentbabl er.****he }
  h
***y z_
 zn[j]   CVCa/*laALLe actua alp
 RPowerR(BIAS2*dsm,cor)L)N+crDDOf ;;
 tn effIfono
  CVCa/or or,e coeff MAL*MS nae 
  nem cv_EAL*MS n(deltaq==1)     for (i!=ap :u[j-1]MAL*t e acM[0] / CVAdju L; }
  Cm cv_EAL**********  * On an order int for ITRfor (i= ric step and each zn[/or or.-1]MALqmBDF ifMALqp1of pr      tmentbabls er.****he }
  h
***ch zns(1 +r rout+1,and peionvVAdjust  Cm ChooseEAL*a new de l a.arorst;em cv_EAL*M array MAL*MS nae 
 ==1) for (i= 2NE, MALqmBD=em m); uteEALqmB**********  MALqp1o=em m); uteEALqp1 DEBUG_CVODE Cm ChooseEALDEBUG_CVODE Cm cv_EAL************* CVSetBDF **************EAL*************************

 This routinevances tn by the  array on anating suMAL*Mccch r;
  o on anaa,ous
 heuris(CVoon rossof h at opCVod li ord. hL_N

 I Colso les th
o ALL_NSstimated lror vector. */

  N_VScale(ON****************************************************/

static int CVSetTqBDF(CVodeMem cv_EAL*
  int j, k;
  
  ttn effIfo ALimaloweem;
yhleshh}
  THRESH,andjeiocsa/or order.     es
 t/* j Stab AL hsTHRESHi >= 1; MAL*t  al }
  ChCVAdju Lh;= 1ale(HUB_FSccoeff turne ALimyo ALL_NSt to L_N returnzn thCVAdju/* j = MAL*t oIN) rese ALL_N ;;
   MAL*/leout)]); 
neash)* L_NVAlt*et = FALSEhCVAdju Lh h;
  nsco{
    CVAdju<Orde***********LL}ut too Res poMALL_NSke a at nm);  is  es
 t/or or,e  to = 0;
ae 
 ==1* CVAdjustAdams ******m); uteEALqmB*****************

 This routine adjusts the historicients lMS n(deltay on anating suMALqmB*ke a 
duction. *  is requested animyo1***************************************************/

static realtype CVAltSum(int iend, rm); uteEALqmB*
  realtype alpha0, alpha0_hatdd){
    cMALqmBD=elpha0_haStabj=2; j < q; jdd)pv, ewt);
  ret zn[j])acorp/**** / hsum;MALqmBD=ecor) RPowerR(BIAS1*dd)  zn[/q)N+crDDOf ;;
 

/********MALqmB)or* CVAdjustAdams ******m); uteEALqp1*****************

 This routine adjusts the historicients lMS n(deltay on anating suMALqpB*ke a 
duction. *order q in ted animyo1***************************************************/

static realtype CVAltSum(int iend, rm); uteEALqpB*
  realtype alpha0, alpha0_hatduicocquot{
    cMALqpBD=elpha0_haStabj=stOrL_N) CVSetAdquotpv,(m p)No/
  ncf, q5 / xRPowerI(h/nv =2 , LS_S
    rej-1], ONE,-dquot;
  for (j=2f(N,  */
  
  ewempv,, duima  ewt);
  ret, tempvacorp/ M[2] hsum;MALqp1*t  alp
  RPowerR(BIAS3*duicocor) L2; jN+crDDOf ;;
 

/********MALqp1*********** CVAdjustBDF **memChooseEAL*****************

 This routine adjusG     MALqmBse ALq,;MALqp1*(lternatins  suMAL******CVAdju 
ut -)in q,;*****+ 1,and peionvVAd)thtes the Nordschooses  * u

r (imumuMAL*natin,rtity MAL* o on***natin,r  to ity *CVAdjuray on aourrespo of lmnating suqCVSIf hange in attih,ay bypr fucror te/

  io**ay ( j ke***h z_sam tdjustOreturn(2)*  is requ at ojustOvariabf the a (3)*order q i at ojustCVSIf han
r (imumuMAL*natin
io**maloweem;
yhleshh}
  THRESH,a at ojustio**kept un/or ord nitiae  * s ties us1***************************************************/

static void CVSetTqBDF(CVodeMem mhooseEALD
  realtype alpha0, alpha0_hatMALL{
    cMALmaleout)MALqmBseout)MALq,;MALqp1 j < q; jStab ALL hsTHRESHi >= 1; MAL*t  al }
  C CVAdju L; }
  COn an order decStab ALL =t e ac :u[j-1]MAL*t e acM[0] / CVAdju L; }
 ale(HUBStab ALL =t e acm1 :u[j-1]MAL*t e acm1M[0] / CVAdju L; -)i;= 1ale(HUB_FSccoMAL*t e acp1M[0] / CVAdju L; + 1DEar== TRU**********)zn[j], zn[]); 
 */
    for (jn;
  }***** CVAltSum ***********, &nflFinear s************

 This routine adjusts the histornrio j)  N_VSmessaorslke a lltheo s  surinear sby
 ******.utin(deltay o**e  re on anatingon***e  re o**ay (deltaqle darraus****************************************************/

static void CVSetAdams(CVm cv_m, &nflFinear knt deltaq)
{
  switchkG_CVO{
  c/orr
  out"striog_lo0 [100E;
  frea2 po zn[qma su abalinetoweighctor. */

  N_Vsl/* j n[jProdturn(SUaco  
  ewempv,n[jAbs(, tempv= N_VWrmsNoITEDiue of lmm).kf, &dsnrio )  N_VSmessaor*MS n(deltaq  N_VSag(cv/* j ase SS (kG_CVO CVSetAdams(urn(kflag);
: q;TPUPUsnrio f(  out"striog_lo0 seoSG(kflag);
S, (doun. ) f_da(doun. ) htin }

 waher;
"CVYd  out"striog_lo0 tin }

 P);
   kflag);
URE)CVSetBDF(cvurn( ONVag);
:q;TPUPUsnrio f(  out"striog_lo0 seoSG( ONVag);
S, (doun. ) f_da(doun. ) htin }

 waher;
"CVYd  out"striog_lo0 tin }

 P);
    ONVag);
URE)CVSetBDF(cvSETUPag);
ED: q;TPUPUsnrio f(  out"striog_lo0 seoSG(SETUPag);
ED, (doun. ) f_tin }

 waher;
"CVYd  out"striog_lo0 tin }

 P);
   SETUPag);
URE)CVSetBDF(cvSOLVEag);
ED: q;TPUPUsnrio f(  out"striog_lo0 seoSG(SOLVEag);
ED, (doun. ) f_tin }

 waher;
"CVYd  out"striog_lo0 tin }

 P);
   SOLVEag);
URE)CVSe

/*********************** CVEwtSetSS *******Stab******************

  This routine rests the historation of l a****Stabetecti turneDeteion *
Algorithm
 STALD

 I Ce the
  orif******a****of h at SLDET opCVod zn.  be a chan ojustio**3;***m th,ay byts to the/ 2))d
  * s   ncfbe a cas it s/***ay (dssaryojustihandne. al, Cdyibeurnmadnn nitiaenough)d
  *htse Vurnz ncf,memslde Ce the
  oCVSIf uirtsumals[q+1stabetection ro violbably mhan ojustio**(dssardn of h at atep
 es
 to********Mccch r;
ly**************************************************/

static void CVSetAodeMem ***StabD
  int j, k;
  
  tn += hh;k, ldf, &ds = [t_ial  
  if (iend q,;sqmBsesqm2;;
   TEtn effIfoojustio**3;***g(dosstOreturntave  hscaleustiv and cd
  ,      push }
  d
  *dowcall iOreturn*** o zn[j] natins ray op.; i <= qmaStabj=2= 3j < q; j++) {)
    for (j3= q; j >= k  { i--) l
  5] += l[]*xi + lssd
 [i][k o}
ssd
 [i-1][k ; q */
  = [t_iali <= ien  i--) l
  i] = M0_i-1;
  alp = [t_iali*= = (-alpsctua = [t_ial*q*  2; *0]tem/q[2] NE,  ;sqmBtua = [t_ial*q* ewt);
  ret zn[j])acorNE,  ;sqm2tua = [t_ial* ewt);
  ret zn[n[j-1acorNE,  ;ssd
 [1][1 o}
sqm2*sqm2;;
   ssd
 [1][2 o}
sqm1*sqm1;;
   ssd
 [1][3 o}
sq*sq********
o{
    CVAdju l[c :u[re coeffIfoojustio**3;***g(dosstOrF ifMnough)ssd
 *htse Vurnz ncf,
 M_Env n***** l[c+5Oreturnit
   tabetection ro deteion *
e histo.i <= qmao{
   abj=2= 3j 
    g**** l[c+5) );u[j-1];
ldf, &o=em slde         break;
  TRU*df, &o> 3j < q; je coeffA1stabetection ro violbablyet tar),  [t_               anorm.

  , &o su4, 5,;***6.            Ressary**********lpha0 realtype alpha0 lo& (q0] / CVAdju L;-1break;
  }MAL*t e acm1M reak;
  }MAL*t oIN) res ALL_N ;;
   
  }MAL*t e a/out)]); neash)* L_NVAlt*et = FALSEALSEhCVAdju Lh*
  nsco{ak;
  opC[NOR o} opC[NOR o+ 1DEar==oeffCVODE,"CVYddNorm error\nsco{ak;
 " Oh zn[f (sl gahe %le is  ***Stab*
 * gamm %l,\hetyp  hs%ora****=f
    }o& (q0] / CVAdj,nionh,h*et = 0 lo& (q0]}(-alpha0 - ale(HUB_FSccoeffets h to thlet*/

  iorder q *htpue , e if d [t_n ds****stabetection ro arnt  =, n****.lpha0 lo& (q***********LL}u******* CVPredict **********slde C******************
  
 CVAltSum returnsts the histordeteios*stabetection robablye quotion thce hscaleALL  tiv and s d
  .***slde C(deltay on amagnitudrivativeALL om tht t/orr= [eris(CVoroot;
r**lT bypr sial Lambd*stabetecto& *n ro t tar),  [t_   
r*o> "sos anuotia *nttlehlessreturn1.0",TEtn MS naductiand ckG_CV. the test quanthou
  only*er,
e
  oriftn ojustio**g(dosste tanh***e of the 3, e ifd
  *htse Vurncolnew ol
f i--)5#endif    entse
  st fai gammtins: d [t_lag(cv_m1 -> Frntd*stabl t/orr= [eris(CVoroot;
or.
 H mait x  CVSet [t_lag(cv_m2 -> Frntd*stabl t/orr= [eris(CVoroot;
 ofr(CVoalinear 
r [t_lag(cv_m3 -> Frntd*stabl t/orr= [eris(CVoroot;
 ofr(CVoalinear ,                i=1 the Newtut*ourreion *
r [t_lag(cv_m4 -> Frntd*stabetectiviolbably mor.
 H mait x  CVSet [t_lag(cv_m5 -> Frntd*stabetectiviolbably m ofr(CVoalinear 
r [t_lag(cv_m6 -> Frntd*stabetectiviolbably m ofr(CVoalinear ,                i=1 the Newtut*ourreion *

r [t_lag(cvretu-> No*stabetection robablyea     realtype alphato  tu
  ne. 
t, fromon robably.

r [t_lag(cv_m-1 -> Min/L_NStbabl er.ssd
 *too*s
 Hl.
r [t_lag(cv_m-2 -> Frrmor.
 H mait x  CVS, vL_NS> vrrt2*vrrt2
r [t_lag(cv_m-3 -> Frrmor.
 H mait x  CVS, Tm;
yhleentbabls                i=1
f prordstepstcie.f  [t_lag(cv_m-4 -> Sm llthduct(1 to /vodeciesle(im thtar syst ofr(CVs.liVScco_lag(cv_m-5 -> Ranating
 It  ofr(CVsane. 
thepstcie.f  [t_lag(cv_m-6u-> No*ourreio garootg /* Rsectiomm).qkammtinsf  [t_lag(cv_m-7u-> Troun. secovr;
 i--)tsusq.f  [t_lag(cv_m-8u-> Troun. secovr;
 i--)B,;***Rania B.f  [t_lag(cv_m-9 -> Rania)tsusq[k odisag(dRse the R*
 It _
  .****************************************************/

static int CVStep(CVodeMem cv_mslde  
  int j, k;
  
  tn += eger(iendiOrk, j,quiOrkmin, lag(cv_m0  
  if (iendtba2] [4    av[4   qkr[4   tsusq[4   tL_N[4   ttL_N[4 rime, CPrimePdrr[4   rrc[4  sqmN[4   qjk[4 [4   vtba2]   qc[6 [4   qco[6 [4   
  if (iendtr  rrcut, vrrtol, vrrt2sesqtol, rrtol  
  if (iend mink  tL_Nk
  
 rv, r 
 r q,;vmin, vL_N rdrrL_N radrr  
  if (iend m llpv= Nsesqm_N rsaqk
 qp rssesqm_Nk rsaqjsesqmi
  
  if (iendrsa, rsb, rsc, rsnt kVS, rd1a, rd1b, rd1c, rd1d  
  if (iendrd2a, rd2b, rd2c, rd3a, rd3b, ctio1,*ourr1M rea if (iendtbap   atNseq = 1seq = 2sebb, rrbout too cm;
minaryuotiatSeiutoffssof h olerntativgencetep se test quan<= qmarrcut    r98  
 vrrtol *****e-4  
 vrrt2tua5**e-4  
 sqtol *****e-3l[0]rrtol *****e-2 < q; jrrD=elpha0_hat too  Indexforourrespo o**ay zn.
deg(dRsystem;
+= erpolbabpts t = l_0 +.0 lo& oo   [t_lv_m1 -> 1 +rpe alpha0 lo& oo   [t_lv_m2 -> qltype alpha0 lo& oo   [t_lv_m3 ->  2;ype alpha0 lo& t too  Indexfi in atbackward-in-endifindex,quv_m1 -> o zn[j] ttPa,q lo& oo   [t_iv_m2 -> vode,oust    ,;MAcpha0 lo& t too get
r (ima, mi
i********mma, tati****** cvom ofr(CVooduct(1 to j)0 lo& t t++) {)[i] k<=3= q; j * (a[i]minko}
ssd
 [1][k ; (a[i]m_NkD=elpha0_haliVSccoi++) l[i] =<=5= sign * (a[i[i]minko}
oIN) mink ssd
 [i][k = FALSEAL]m_NkD=eout)tL_Nk
ssd
 [i][k = FALSEq */
    al TRU]minko< TINY*tL_Nkn * (a[i[ilag(cv_m-1;liVScco_********kG_CVO FALSEq */
 tL_N[k o}
sL_Nk;;
   ssL_N[k o}
sL_Nk*sL_Nk;;
    (a[i] /gamma lpha0_hali 
 r qma lpha0_halii++) l[i] =<=4= sign * (a[i[ir
 [i][k o}
ssd
 [i][k /ssd
 [i+1][k ; (a[i[i] /gamma ] /gamm+ir
 [i][k ; (a[i[i] /g qma ] /g qm+ir
 [i][k *r
 [i][k ; (a[i}  (a[i av[k o}
FOURTH*] /gam; (a[ivtba2k];
 reasFOURTH*] /gs; -) av[k * av[k S_STE
 re coqc[5][k o}
ssd
 [1][k *ssd
 [3][k o-
ssd
 [2][k *ssd
 [2][k ; (a[iqc[4 [k o}
ssd
 [2][k *ssd
 [3][k o-
ssd
 [1][k *ssd
 [4][k ; (a[iqc[3 [k o}
lpha0_haliqc[2 [k o}
ssd
 [2][k *ssd
 [5][k o-
ssd
 [3][k *ssd
 [4][k ; (a[iqc[1 [k o}
ssd
 [4][k *ssd
 [4][k o-
ssd
 [3][k *ssd
 [5][k ; (a[iVSccoi++) l[i] =<=5= sign * (a[i[iqco[i][k o}
qc[i][k ; (a[i}********************************or Er e thkh**sig lo& t too Iso 
ticor.
 H rrmoearly-or.
 H mait x  CVS. th(dRs ofr(CVowill (a[i[have 
t,mon rrmoearly-
t,mon rootse
 tese theo .iVScco_st faila_lag(cv_m1 ithe locprocedar s****s.utfm -leentootgVScco_diffuccmwtSe tanhvrrt2se(deltaq  N_VSlag(cv_m-3.pha0 lo& t tvmino}
oIN)vtba21],oIN)vtba22],vtba23] j < qvL_NST)out)vtba21],out)vtba22],vtba23] j < q
o{
 (vmino< vrrtol*vrrtol1) && (qwaitvL_NS> vrrt2*vrrt2n * (a[i[ilag(cv_m-2;liVScco_********kG_CVO FALSEqle(HUB_FScco jrrD=e( av[1]m+ir
v[2 o+ir
v[3] /THREE;;
   TEtn cco_drrL_No}
lpha0_hali i--) )
    k<=3=kign * (a[i[i radrr;
 reas av[k o-) r= FALSEALSEdrrL_No}
out)drrL_N radrr= FALSEAL}reak;
  TRUdrrL_No> vrrt2n * (a[i[i rlag(cv_m-3;  TEtn cco_};
   TEtn cco_lag(cv_m1 F */
 & oo  can 
t, from/orr= [is(CVoroot;
dro**he nm);  eion *
a0 lo& (q0] (a[i}*****e(HUB_FFSccoeffuq i at  ofr(CVsahe get
r**l lo& (q&& (qwaitreasqco[1][1 )o< TINY*tsL_N[1 )o* (a[i[i]m llt= qco[1][1 ; (a[i[ilag(cv_m-4;  TEtn cco_*******kG_CVO FALSEq */
 [0] / emt= qco[1][2]/qco[1][1 ; (a[i--) l[2] =<=5= sign * (a[i[iqco[i][2]t= qco[i][2]t-/ em*qco[i][1]_VSccoa re coqco[1][2]o}
lpha0_hali emt= qco[1][3]/qco[1][1 ; (a[i--) l[2] =<=5= sign * (a[i[iqco[i][3]t= qco[i][3]t-/ em*qco[i][1]_VSccoa e coqco[1][3]D=elpha0_haliVSccowaitreasqco[2][2])o< TINY*tsL_N[2 )o* (a[i[i]m llt= qco[2][2]; (a[i[ilag(cv_m-4;  TEtn cco_*******kG_CVO FALSEq */
 [0] / emt= qco[2][3]/qco[2][2]; (a[i--) l[3] =<=5= sign * (a[i[iqco[i][3]t= qco[i][3]t-/ em*qco[i][2] FALSEq */
    al TRUreasqco[4][3])o< TINY*tsL_N[3 )o* (a[i[i]m llt= qco[4][3]; (a[i[ilag(cv_m-4;  TEtn cco_*******kG_CVO FALSEq */
 [0] /rr;
 -qco[5][3]/qco[4][3]; (a[io& (qwaitPro< TINYoutpr*o> HUNn * (a[i[ilag(cv_m-5; TEtn cco_*******kG_CVO FALSEq */
 [0] /--) )[i] k<=3= q; j * (a[i  qkr[k o}
qc[5][k o+irr*(qc[4 [k o+irr*rr*(qc[2 [k o+irr*qc[1 [k )O FALSEql ;
    (a[i]qL_No}
lpha0_hali--) )[i] k<=3= q; j * (a[i  saqk;
 reasqkr[k )/ssL_N[k break;
  TRUsaqk;>i]qL_N)i]qL_No}
saqk; (a[i}  (a[isqmi
o}
sqmaxDEar== TRU]qL_No< sqtoln * (a[i[ilag(cv_m2;;
   TEtn cco_oo  can 
t, from/orr= [is(CVoroot;
dro**he "      tr MAc"
a0 lo& (q0] (a[i}*e(HUB_FFSccoo_oo do Newtut*ourreion *sahe improve
r**l0 lo& (q0] (a[i[i--) lt[i] =t<=3= itign * (a[i[i r--) )[i] k<=3= q; j * (a[i  0] / Co}
qc[4 [k o+irr*rr*(THREE*qc[2 [k o+irr*FOUR*qc[1 [k ); (a[i  0] /drr[k o}
lpha0_halieak;
  TRUreasqpv;> TINY*tsL_N[k )/drr[k o}
-qkr[k /qp0_halieak;
 rrc[k o}
r*o+/drr[k 0_halieak;}a     real (a[i[i r--) )[i] k<=3= q; j * (a[i  0] /so}
r*c[k 0_halieak;[i]qL_Nko}
lpha0_halieak;
 --) 
    i<=3= j; j * (a[i  0] /  qjk[j][k o}
qc[5][j o+is*(qc[4 [j o+i                i=1
 xi_i = [t_n -ak;[i]*s*(qc[2 [j o+is*qc[1 [j )O FALSEn -ak;[i]aqj;
 reasqjk[j][k )/ssL_N[j] FALSEn -ak;[i TRUsaqj;>i]qL_Nk)i]qL_Nko}
saqj FALSEn -ak;} _halieak;[i]qLN[k o}
sqL_Nk;;
   -ak;} _
lieak;[i]qLi
o}
sqmN[1 ;rkmini <= ien  [i r--) )[2] k<=3= q; j * (a[i  0] / TRU]qLN[k o<i]qLi
j * (a[i  0] /  kmini <k FALSEn -ak;[i]qLi
o}
sqmN[k 0_halieak;[i}_halieak;}a     realrr;
 r*c[kLi
 0_halieak;
[i  0] / TRU]qLino< sqtoln * (a[i[i [t_lag(cv_m30_halieak;[ioo  can 
t, from/orr= [is(CVorootha0 lo& (q0] /[ioo  breakfCVOe thNewtut*ourreion *h**siga ifdro**he "      tr MAc"
*/ _halieak;[ibreak;;
   -ak;} e(HUB_FScco jk;
 --) 
    i<=3= j; j * (a[i  0] /  qkr[j o}
qjk[j][kLi
 0_halieak;[i}_halieak;}a  TEtn cco_}& (q0] /[ioo  ar e thNewtut*ourreion *h**sig
*/ _haliea
  0] / TRU]qLino> sqtoln * (a[i[i [lag(cv_m-60_halieak;*******kG_CVO FALSE0]}(-alph] /[ioo  ar e th TRU]qL_No< sqtoln e(HUBa0 lo& } ] /[ioo  ar e th T(vmino< vrrtol*vrrtol1)e(HU,  ofr(CVsahe get
r**l lo& t too g     tr bf td)tsusq[k o****merify
r**l0 lo& effAlltuctiand ckG_CVfdro**he e loc eion *
a lo& t t++) )[i] k<=3= q; j * (a[irsao}
ssd
 [1][k ; (a[irsbo}
ssd
 [2][k *rr  
 [irsco}
ssd
 [3][k *rr*rr  
 [irsdo}
ssd
 [4][k *rr*rr*rr  
 [irseo}
ssd
 [5][k *rr*rr*rr*rr  
 [ird1ao}
rsao-irsb  
 [ird1bo}
rsbo-irsc  
 [ird1co}
rsco-irsd  
 [ird1do}
rsdo-irse  
 [ird2ao}
rd1ao-ird1b  
 [ird2bo}
rd1bo-ird1c  
 [ird2co}
rd1co-ird1d  
 [ird3ao}
rd2ao-ird2b  
 [ird3bo}
rd2bo-ird2c0_haliVSccowaitreasrd1b)o< TINY*tL_N[k )/* (a[i[ilag(cv_m-7;tn cco_*******kG_CVO FALSEq */
 [0] /ctio1v_m-rd3a/rd1b  
 [i TRUEtio1v< TINYoutpEtio1v>
FOUR)/* (a[i[ilag(cv_m-7;tn cco_*******kG_CVO FALSEq */
 ourr1D=e( d2b/Etio1)/(rr*rrn;
  }
  rsq[k o}
ssd
 [3][k o+iourr1Mrder int f TRU] rsq[2]v< TINY)/* (a[ilag(cv_m-8**  * On an (kG_CVO FALr int ftbapo}
s rsq[3]/] rsq[2];t ftbamo}
s rsq[1]/] rsq[2];t fq = 1o}
FOURTH*(q*; -)_inv;t fq = 2tuaTWO/(; -)_inv;t fbbo}
rbap*tbamo-+ alphaq = 1*tbap0_ha emt=  alphaq = 2*bb < q; jStabreas em)v< TINY)/* (a[ilag(cv_m-8**  * On an (kG_CVO FALr int ftrbD=ecor) em < q; jStabreastrbD-) r=v>
rrtol1) && (qlag(cv_m-9**  * On an (kG_CVO FALr int fs foheckp***seUBStartio**above
iutoffarrcut t/* j Stabr*o> rrcut1) && (qwaitkAIL */
  j kag(cv_m4;&& (qwaitkAIL */
 2j kag(cv_m5;&& (qwaitkAIL */
 3j kag(cv_m6order int for Alltuctiand ckG_CVfof thealeatns
    to t
a lo& t tOn an (kG_CVO FALu********** CVEwtSetSV ********************************/

static int CVStep****** CVE END Ptiv ae Helper tation *s I; 
#iciebably **********

 /******* CVEwtSetSV ********************************/

static int CVStep******** CVEwtSetSV ********************************/

static int /******* CVEwtSet END rror\ I; 
#iciebably **********

 *******

 /******* CVEwtSetSV ********************************/

static int /
