/* GPD.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* *********************************************************************** */
/* HERE WE GATHER SEVERAL FUNCTIONS OF (x,Q2,kT) COMPUTED USING THE LFWF */
/* FORMALISM WITH SCALAR VERTICES; THESE MAY THEN BE CALLED TO COMPUTE, */
/* E.G., IC PDFs AND THE CHARM SIGMA TERM. */

/*     LAST EDITED: TJH, SEPT. 27th, 2016. */
/* *********************************************************************** */
doublereal f1int_(doublereal *q2, doublereal *x, doublereal *kt, doublereal *
	nq, doublereal *lq, doublereal *mq, doublereal *ms, integer *typ)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3, d__4, d__5, d__6, d__7;

    /* Builtin functions */
    double atan(doublereal), exp(doublereal), sqrt(doublereal);

    /* Local variables */
    static doublereal kpkm_sum__, pi, mn, wf, x1f, x2f, kt2, x3f, sinv, aterm,
	     bterm;


/*  SPECIFIES THE INTEGRAND OF THE x AND kT MOMENTS THAT DETERMINE F1(Q2) */
/*      THIS IS A FUNCTION OF Q2, ... AND {x, kT} ARE TO BE INTEGRATED */
/*      OVER IN THE CALLING PROGRAM */

/*  WRITTEN: T. HOBBS (OCT. 2014) */
/* *********************************************************************** */
    pi = atan(1.) * 4;
    mn = .9382720813;

/* .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . */
/* PROTON MASS; IN GeV! */
/* Computing 2nd power */
    d__1 = *kt;
    kt2 = d__1 * d__1;
/* Computing 2nd power */
    d__1 = *mq;
/* Computing 2nd power */
    d__2 = *ms;
/* Computing 2nd power */
    d__3 = 1. - *x;
    sinv = (kt2 + (1. - *x) * (d__1 * d__1) + *x * (d__2 * d__2) + d__3 * 
	    d__3 * *q2 / 4.) / (*x * (1. - *x));
/* Computing 2nd power */
    d__1 = 1. - *x;
    kpkm_sum__ = (kt2 + d__1 * d__1 * *q2 / 4.) * 2.;
/* Computing 4th power */
    d__1 = *lq, d__1 *= d__1;
/* Computing 2nd power */
    d__2 = *x;
/* Computing 2nd power */
    d__3 = 1. - *x;
    x1f = 1. / (d__1 * d__1 * 4.) * (1. / (d__2 * d__2) + 1. / (d__3 * d__3) 
	    + 2. / (*x * (1. - *x)));
/* Computing 4th power */
    d__1 = *lq, d__1 *= d__1;
/* Computing 2nd power */
    d__2 = *mq / *x;
/* Computing 2nd power */
    d__3 = *ms / (1. - *x);
/* Computing 2nd power */
    d__4 = *mq;
/* Computing 2nd power */
    d__5 = *ms;
    x2f = 1. / (d__1 * d__1 * 4.) * (d__2 * d__2 + d__3 * d__3 + (d__4 * d__4 
	    + d__5 * d__5) / (*x * (1. - *x)));
/* Computing 4th power */
    d__1 = *lq, d__1 *= d__1;
/* Computing 4th power */
    d__2 = *mq, d__2 *= d__2;
/* Computing 2nd power */
    d__3 = *x;
/* Computing 4th power */
    d__4 = *ms, d__4 *= d__4;
/* Computing 2nd power */
    d__5 = 1. - *x;
/* Computing 2nd power */
    d__6 = *mq;
/* Computing 2nd power */
    d__7 = *ms;
    x3f = 1. / (d__1 * d__1 * 4.) * (d__2 * d__2 / (d__3 * d__3) + d__4 * 
	    d__4 / (d__5 * d__5) + d__6 * d__6 * (d__7 * d__7) * 2. / (*x * (
	    1. - *x)));
/* Computing 2nd power */
    d__1 = *lq;
/* Computing 2nd power */
    d__2 = kpkm_sum__ / 2.;
    aterm = sinv / (d__1 * d__1) + 1. + x1f * (d__2 * d__2) + x2f * 
	    kpkm_sum__ + x3f;
/* Computing 2nd power */
    d__1 = 1. - *x;
    bterm = d__1 * d__1 * kt2 * *q2 * x1f;
/* .....  CHOICES FOR THE FUNCTIONAL FORM OF THE CHARM WAVEFUNCTION.... */
    if (*typ == 1) {
/* Computing 2nd power */
	d__1 = *lq;
	wf = exp(-sinv / (d__1 * d__1));

    } else if (*typ == 2) {
	wf = 0.;

    } else if (*typ == 3) {
/* Computing 3rd power */
	d__1 = aterm;
/* Computing 3rd power */
	d__2 = sqrt(aterm / (aterm - bterm));
	wf = (aterm * 2. - bterm) * .5 / (d__1 * (d__1 * d__1)) * (d__2 * (
		d__2 * d__2));

    }
/* Computing 4th power */
    d__1 = *lq, d__1 *= d__1;
/* Computing 2nd power */
    d__2 = pi;
/* Computing 2nd power */
    d__3 = *mq + *x * mn;
/* Computing 2nd power */
    d__4 = 1. - *x;
/* Computing 2nd power */
    d__5 = *x;
    ret_val = *nq / (d__1 * d__1) * (1. / (d__2 * d__2 * 16.)) * (kt2 + d__3 *
	     d__3 - d__4 * d__4 * *q2 / 4.) * (1. / (d__5 * d__5 * (1. - *x)))
	     * 2. * *kt * wf;
/* A DIMENSIONLESS */
    return ret_val;
} /* f1int_ */

/* *************************************************************************** */
/* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */
doublereal f2int_(doublereal *q2, doublereal *x, doublereal *kt, doublereal *
	nq, doublereal *lq, doublereal *mq, doublereal *ms, integer *typ)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3, d__4, d__5, d__6, d__7;

    /* Builtin functions */
    double atan(doublereal), exp(doublereal), sqrt(doublereal);

    /* Local variables */
    static doublereal kpkm_sum__, pi, mn, wf, x1f, x2f, kt2, x3f, sinv, aterm,
	     bterm;


/*  SPECIFIES THE INTEGRAND OF THE x AND kT MOMENTS THAT DETERMINE F2(Q2) */
/*      THIS IS A FUNCTION OF Q2, ... AND {x, kT} ARE TO BE INTEGRATED */
/*      OVER IN THE CALLING PROGRAM */

/*  WRITTEN: T. HOBBS (OCT. 2014) */
/* *********************************************************************** */
    pi = atan(1.) * 4;
    mn = .9382720813;

/* .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . */
/* KEEP ALL MASSES IN GeV! */
/* Computing 2nd power */
    d__1 = *kt;
    kt2 = d__1 * d__1;
/* Computing 2nd power */
    d__1 = *mq;
/* Computing 2nd power */
    d__2 = *ms;
/* Computing 2nd power */
    d__3 = 1. - *x;
    sinv = (kt2 + (1. - *x) * (d__1 * d__1) + *x * (d__2 * d__2) + d__3 * 
	    d__3 * *q2 / 4.) / (*x * (1. - *x));
/* Computing 2nd power */
    d__1 = 1. - *x;
    kpkm_sum__ = (kt2 + d__1 * d__1 * *q2 / 4.) * 2.;
/* Computing 4th power */
    d__1 = *lq, d__1 *= d__1;
/* Computing 2nd power */
    d__2 = *x;
/* Computing 2nd power */
    d__3 = 1. - *x;
    x1f = 1. / (d__1 * d__1 * 4.) * (1. / (d__2 * d__2) + 1. / (d__3 * d__3) 
	    + 2. / (*x * (1. - *x)));
/* Computing 4th power */
    d__1 = *lq, d__1 *= d__1;
/* Computing 2nd power */
    d__2 = *mq / *x;
/* Computing 2nd power */
    d__3 = *ms / (1. - *x);
/* Computing 2nd power */
    d__4 = *mq;
/* Computing 2nd power */
    d__5 = *ms;
    x2f = 1. / (d__1 * d__1 * 4.) * (d__2 * d__2 + d__3 * d__3 + (d__4 * d__4 
	    + d__5 * d__5) / (*x * (1. - *x)));
/* Computing 4th power */
    d__1 = *lq, d__1 *= d__1;
/* Computing 4th power */
    d__2 = *mq, d__2 *= d__2;
/* Computing 2nd power */
    d__3 = *x;
/* Computing 4th power */
    d__4 = *ms, d__4 *= d__4;
/* Computing 2nd power */
    d__5 = 1. - *x;
/* Computing 2nd power */
    d__6 = *mq;
/* Computing 2nd power */
    d__7 = *ms;
    x3f = 1. / (d__1 * d__1 * 4.) * (d__2 * d__2 / (d__3 * d__3) + d__4 * 
	    d__4 / (d__5 * d__5) + d__6 * d__6 * (d__7 * d__7) * 2. / (*x * (
	    1. - *x)));
/* Computing 2nd power */
    d__1 = *lq;
/* Computing 2nd power */
    d__2 = kpkm_sum__ / 2.;
    aterm = sinv / (d__1 * d__1) + 1. + x1f * (d__2 * d__2) + x2f * 
	    kpkm_sum__ + x3f;
/* Computing 2nd power */
    d__1 = 1. - *x;
    bterm = d__1 * d__1 * kt2 * *q2 * x1f;
/* .....  CHOICES FOR THE FUNCTIONAL FORM OF THE CHARM WAVEFUNCTION.... */
    if (*typ == 1) {
/* Computing 2nd power */
	d__1 = *lq;
	wf = exp(-sinv / (d__1 * d__1));

    } else if (*typ == 2) {
	wf = 0.;

    } else if (*typ == 3) {
/* Computing 3rd power */
	d__1 = aterm;
/* Computing 3rd power */
	d__2 = sqrt(aterm / (aterm - bterm));
	wf = (aterm * 2. - bterm) * .5 / (d__1 * (d__1 * d__1)) * (d__2 * (
		d__2 * d__2));

    }
/* Computing 4th power */
    d__1 = *lq, d__1 *= d__1;
/* Computing 2nd power */
    d__2 = pi;
/* Computing 2nd power */
    d__3 = *x;
    ret_val = *nq / (d__1 * d__1) * (mn / (d__2 * d__2 * 8.)) * ((*mq + *x * 
	    mn) / (d__3 * d__3)) * 2. * *kt * wf;
/* A DIMENSIONLESS RED */
    return ret_val;
} /* f2int_ */

/* *************************************************************************** */
/* *************************************************************************** */
/*    .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . */
doublereal gaint_(doublereal *q2, doublereal *x, doublereal *kt, doublereal *
	nq, doublereal *lq, doublereal *mq, doublereal *ms, integer *typ)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, 
	    d__10, d__11, d__12;

    /* Builtin functions */
    double atan(doublereal), exp(doublereal);

    /* Local variables */
    static doublereal kpkm_sum__, pi, mn, wf, kt2, mq4, sinv, kpkm2;


/*  SPECIFIES THE INTEGRAND OF THE x AND kT MOMENTS THAT DETERMINE GA(Q2) */
/*      THIS IS A FUNCTION OF Q2, ... AND {x, kT} ARE TO BE INTEGRATED */
/*      OVER IN THE CALLING PROGRAM */

/*  WRITTEN: T. HOBBS (Jan. 2015) */
/* *********************************************************************** */
    pi = atan(1.) * 4;
    mn = .9382720813;

/* .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . */
/* KEEP ALL MASSES IN GeV! */
/* Computing 2nd power */
    d__1 = *kt;
    kt2 = d__1 * d__1;
/* Computing 2nd power */
    d__1 = *mq;
/* Computing 2nd power */
    d__2 = *ms;
/* Computing 2nd power */
    d__3 = 1. - *x;
    sinv = (kt2 + (1. - *x) * (d__1 * d__1) + *x * (d__2 * d__2) + d__3 * 
	    d__3 * *q2 / 4.) / (*x * (1. - *x));
/* Computing 2nd power */
    d__1 = 1. - *x;
    kpkm_sum__ = (kt2 + d__1 * d__1 * *q2 / 4.) * 2.;
/* Computing 2nd power */
    d__2 = 1. - *x;
/* Computing 2nd power */
    d__1 = kt2 - d__2 * d__2 * *q2 / 4.;
    kpkm2 = d__1 * d__1;
/* Computing 2nd power */
    d__1 = *x;
/* Computing 2nd power */
    d__2 = 1. - *x;
/* Computing 2nd power */
    d__3 = *mq / *x;
/* Computing 2nd power */
    d__4 = *ms / (1. - *x);
/* Computing 2nd power */
    d__5 = *mq;
/* Computing 2nd power */
    d__6 = *ms;
/* Computing 4th power */
    d__7 = *mq, d__7 *= d__7;
/* Computing 2nd power */
    d__8 = *x;
/* Computing 4th power */
    d__9 = *ms, d__9 *= d__9;
/* Computing 2nd power */
    d__10 = 1. - *x;
/* Computing 2nd power */
    d__11 = *mq;
/* Computing 2nd power */
    d__12 = *ms;
    mq4 = kpkm2 * (1 / (d__1 * d__1) + 1 / (d__2 * d__2) + 2. / (*x * (1. - *
	    x))) + kpkm_sum__ * (d__3 * d__3 + d__4 * d__4 + (d__5 * d__5 + 
	    d__6 * d__6) / (*x * (1. - *x))) + d__7 * d__7 / (d__8 * d__8) + 
	    d__9 * d__9 / (d__10 * d__10) + d__11 * d__11 * (d__12 * d__12) * 
	    2. / (*x * (1. - *x));
/* .....  CHOICES FOR THE FUNCTIONAL FORM OF THE CHARM WAVEFUNCTION.... */
    if (*typ == 1) {
/* Computing 2nd power */
	d__1 = *lq;
	wf = exp(-sinv / (d__1 * d__1));

    } else if (*typ == 2) {
/* Computing 2nd power */
	d__1 = *lq;
/* Computing 4th power */
	d__2 = *lq, d__2 *= d__2;
	wf = 1. / (sinv / (d__1 * d__1) + 1. + mq4 / (d__2 * d__2 * 4.));

    } else if (*typ == 3) {
/* Computing 2nd power */
	d__2 = *lq;
/* Computing 4th power */
	d__3 = *lq, d__3 *= d__3;
/* Computing 2nd power */
	d__1 = sinv / (d__2 * d__2) + 1. + mq4 / (d__3 * d__3 * 4.);
	wf = 1. / (d__1 * d__1);

    }
/* Computing 4th power */
    d__1 = *lq, d__1 *= d__1;
/* Computing 2nd power */
    d__2 = pi;
/* Computing 2nd power */
    d__3 = *mq + *x * mn;
/* Computing 2nd power */
    d__4 = 1. - *x;
/* Computing 2nd power */
    d__5 = *x;
    ret_val = *nq / (d__1 * d__1) * (1. / (d__2 * d__2 * 16.)) * (-kt2 + d__3 
	    * d__3 + d__4 * d__4 * *q2 / 4.) * (1. / (d__5 * d__5 * (1. - *x))
	    ) * 2. * *kt * wf;
/* A DIMENSIONLESS */
    return ret_val;
} /* gaint_ */

/* *************************************************************************** */
/*    .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . */
/* *********************************************************************** */
doublereal ccint_(doublereal *q2, doublereal *x, doublereal *kt, doublereal *
	nq, doublereal *lq, doublereal *mq, doublereal *ms, integer *typ)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, 
	    d__10, d__11, d__12;

    /* Builtin functions */
    double atan(doublereal), exp(doublereal);

    /* Local variables */
    static doublereal kpkm_sum__, pi, mn, wf, kt2, mq4, sinv, kpkm2;


/*  SPECIFIES THE INTEGRAND OF THE x AND kT MOMENTS THAT DETERMINE THE */
/*      SCALAR DENSITY <cbar c> */

/*  WRITTEN: T. HOBBS (FEB. 2015) */
/*  MODIFIED:         (SEP. 2016) */
/* *********************************************************************** */
    pi = atan(1.) * 4;
    mn = .9382720813;

/* .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . */
/* KEEP ALL MASSES IN GeV! */
/* Computing 2nd power */
    d__1 = *kt;
    kt2 = d__1 * d__1;
/* Computing 2nd power */
    d__1 = *mq;
/* Computing 2nd power */
    d__2 = *ms;
/* Computing 2nd power */
    d__3 = 1. - *x;
    sinv = (kt2 + (1. - *x) * (d__1 * d__1) + *x * (d__2 * d__2) + d__3 * 
	    d__3 * *q2 / 4.) / (*x * (1. - *x));
/* Computing 2nd power */
    d__1 = 1. - *x;
    kpkm_sum__ = (kt2 + d__1 * d__1 * *q2 / 4.) * 2.;
/* Computing 2nd power */
    d__2 = 1. - *x;
/* Computing 2nd power */
    d__1 = kt2 - d__2 * d__2 * *q2 / 4.;
    kpkm2 = d__1 * d__1;
/* Computing 2nd power */
    d__1 = *x;
/* Computing 2nd power */
    d__2 = 1. - *x;
/* Computing 2nd power */
    d__3 = *mq / *x;
/* Computing 2nd power */
    d__4 = *ms / (1. - *x);
/* Computing 2nd power */
    d__5 = *mq;
/* Computing 2nd power */
    d__6 = *ms;
/* Computing 4th power */
    d__7 = *mq, d__7 *= d__7;
/* Computing 2nd power */
    d__8 = *x;
/* Computing 4th power */
    d__9 = *ms, d__9 *= d__9;
/* Computing 2nd power */
    d__10 = 1. - *x;
/* Computing 2nd power */
    d__11 = *mq;
/* Computing 2nd power */
    d__12 = *ms;
    mq4 = kpkm2 * (1 / (d__1 * d__1) + 1 / (d__2 * d__2) + 2. / (*x * (1. - *
	    x))) + kpkm_sum__ * (d__3 * d__3 + d__4 * d__4 + (d__5 * d__5 + 
	    d__6 * d__6) / (*x * (1. - *x))) + d__7 * d__7 / (d__8 * d__8) + 
	    d__9 * d__9 / (d__10 * d__10) + d__11 * d__11 * (d__12 * d__12) * 
	    2. / (*x * (1. - *x));
/* .....  CHOICES FOR THE FUNCTIONAL FORM OF THE CHARM WAVEFUNCTION.... */
    if (*typ == 1) {
/* Computing 2nd power */
	d__1 = *lq;
	wf = exp(-sinv / (d__1 * d__1));

    } else if (*typ == 2) {
/* Computing 2nd power */
	d__1 = *lq;
/* Computing 4th power */
	d__2 = *lq, d__2 *= d__2;
	wf = 1. / (sinv / (d__1 * d__1) + 1. + mq4 / (d__2 * d__2 * 4.));

    } else if (*typ == 3) {
/* Computing 2nd power */
	d__2 = *lq;
/* Computing 4th power */
	d__3 = *lq, d__3 *= d__3;
/* Computing 2nd power */
	d__1 = sinv / (d__2 * d__2) + 1. + mq4 / (d__3 * d__3 * 4.);
	wf = 1. / (d__1 * d__1);

    }
/* Computing 4th power */
    d__1 = *lq, d__1 *= d__1;
/* Computing 2nd power */
    d__2 = pi;
/* Computing 2nd power */
    d__3 = *mq + *x * mn;
/* Computing 2nd power */
    d__4 = 1. - *x;
/* Computing 3rd power */
    d__5 = *x;
    ret_val = *nq / (d__1 * d__1) * (1. / (d__2 * d__2 * 16.)) * (kt2 + d__3 *
	     d__3 - d__4 * d__4 * *q2 / 4.) * (1. / (d__5 * (d__5 * d__5) * (
	    1. - *x))) * 2. * *kt * wf;
/* A DIMENSIONLESS */
    return ret_val;
} /* ccint_ */

