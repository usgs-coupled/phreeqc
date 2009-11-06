/*******************************************************************
 *                                                                 *
 * File          : nvector.c                                       *
 * Programmers   : Radu Serban, LLNL                               *
 * Version of    : 26 June 2002                                    *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California *
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/shared/LICENSE                        *
 *-----------------------------------------------------------------*
 * This is the implementation file for a generic NVECTOR           *
 * package. It contains the implementation of the N_Vector         *
 * kernels listed in nvector.h.                                    *
 *                                                                 *
 *******************************************************************/
#if !defined(PHREEQC_CLASS)
#define CLASS_QUALIFIER
#else
#include "Phreeqc.h"
#define CLASS_QUALIFIER Pheeqc::
#endif
#include "nvector.h"			/* generic M_Env and N_Vector */
#include "output.h"
#if !defined(PHREEQC_CLASS)
static char const svnid[] =
	"$Id$";
#endif
N_Vector 
N_VNew(integertype n, M_Env machEnv)
{
	N_Vector v_new;
	v_new = machEnv->ops->nvnew(n, machEnv);
	return (v_new);
}

N_Vector_S CLASS_QUALIFIER
N_VNew_S(integertype ns, integertype n, M_Env machEnv)
{
	N_Vector_S vs_new;
	vs_new = machEnv->ops->nvnewS(ns, n, machEnv);
	return (vs_new);
}

void CLASS_QUALIFIER
N_VFree(N_Vector v)
{
	v->menv->ops->nvfree(v);
}

void CLASS_QUALIFIER
N_VFree_S(integertype ns, N_Vector_S vs)
{
	(*vs)->menv->ops->nvfreeS(ns, vs);
}

N_Vector CLASS_QUALIFIER
N_VMake(integertype n, realtype * v_data, M_Env machEnv)
{
	N_Vector v_new;
	v_new = machEnv->ops->nvmake(n, v_data, machEnv);
	return (v_new);
}

void CLASS_QUALIFIER
N_VDispose(N_Vector v)
{
	v->menv->ops->nvdispose(v);
}

realtype * CLASS_QUALIFIER
N_VGetData(N_Vector v)
{
	realtype *data;
	data = v->menv->ops->nvgetdata(v);
	return (data);
}

void CLASS_QUALIFIER
N_VSetData(realtype * v_data, N_Vector v)
{
	v->menv->ops->nvsetdata(v_data, v);
}

void CLASS_QUALIFIER
N_VLinearSum(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
	z->menv->ops->nvlinearsum(a, x, b, y, z);
}

void CLASS_QUALIFIER
N_VConst(realtype c, N_Vector z)
{
	z->menv->ops->nvconst(c, z);
}

void CLASS_QUALIFIER
N_VProd(N_Vector x, N_Vector y, N_Vector z)
{
	z->menv->ops->nvprod(x, y, z);
}

void CLASS_QUALIFIER
N_VDiv(N_Vector x, N_Vector y, N_Vector z)
{
	z->menv->ops->nvdiv(x, y, z);
}

void CLASS_QUALIFIER
N_VScale(realtype c, N_Vector x, N_Vector z)
{
	z->menv->ops->nvscale(c, x, z);
}

void CLASS_QUALIFIER
N_VAbs(N_Vector x, N_Vector z)
{
	z->menv->ops->nvabs(x, z);
}

void CLASS_QUALIFIER
N_VInv(N_Vector x, N_Vector z)
{
	z->menv->ops->nvinv(x, z);
}

void CLASS_QUALIFIER
N_VAddConst(N_Vector x, realtype b, N_Vector z)
{
	z->menv->ops->nvaddconst(x, b, z);
}

realtype CLASS_QUALIFIER
N_VDotProd(N_Vector x, N_Vector y)
{
	realtype prod;
	prod = y->menv->ops->nvdotprod(x, y);
	return (prod);
}

realtype CLASS_QUALIFIER
N_VMaxNorm(N_Vector x)
{
	realtype norm;
	norm = x->menv->ops->nvmaxnorm(x);
	return (norm);
}

realtype CLASS_QUALIFIER
N_VWrmsNorm(N_Vector x, N_Vector w)
{
	realtype norm;
	norm = x->menv->ops->nvwrmsnorm(x, w);
	return (norm);
}

realtype CLASS_QUALIFIER
N_VMin(N_Vector x)
{
	realtype minval;
	minval = x->menv->ops->nvmin(x);
	return (minval);
}

realtype CLASS_QUALIFIER
N_VWL2Norm(N_Vector x, N_Vector w)
{
	realtype norm;
	norm = x->menv->ops->nvwl2norm(x, w);
	return (norm);
}

realtype CLASS_QUALIFIER
N_VL1Norm(N_Vector x)
{
	realtype norm;
	norm = x->menv->ops->nvl1norm(x);
	return (norm);
}

void CLASS_QUALIFIER
N_VOneMask(N_Vector x)
{
	x->menv->ops->nvonemask(x);
}

void CLASS_QUALIFIER
N_VCompare(realtype c, N_Vector x, N_Vector z)
{
	z->menv->ops->nvcompare(c, x, z);
}

booleantype CLASS_QUALIFIER
N_VInvTest(N_Vector x, N_Vector z)
{
	booleantype flag;
	flag = z->menv->ops->nvinvtest(x, z);
	return (flag);
}

booleantype CLASS_QUALIFIER
N_VConstrProdPos(N_Vector c, N_Vector x)
{
	booleantype flag;
	flag = x->menv->ops->nvconstrprodpos(c, x);
	return (flag);
}

booleantype CLASS_QUALIFIER
N_VConstrMask(N_Vector c, N_Vector x, N_Vector m)
{
	booleantype flag;
	flag = x->menv->ops->nvconstrmask(c, x, m);
	return (flag);
}

realtype CLASS_QUALIFIER
N_VMinQuotient(N_Vector num, N_Vector denom)
{
	realtype quotient;
	quotient = num->menv->ops->nvminquotient(num, denom);
	return (quotient);
}

void CLASS_QUALIFIER
N_VPrint(N_Vector x)
{
	x->menv->ops->nvprint(x);
}
