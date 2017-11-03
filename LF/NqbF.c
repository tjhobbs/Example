/* NqbF.f -- translated by f2c (version 20100827).
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

/* Table of constant values */

static doublereal c_b2 = 0.;
static doublereal c_b3 = 1.;

/* ********************************************************************** */
doublereal nqbf_(doublereal *nq, doublereal *mq, doublereal *lq, doublereal *
	lqb, doublereal *ms, doublereal *msb, integer *typ)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double atan(doublereal);

    /* Local variables */
    static doublereal w, pi, mn, is;
    static integer iw;
    static doublereal isb, mqb, kti, is_i__;
    static integer ikti;
    static doublereal wmin, wmax, wint;
    extern doublereal f1int_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    ;
    static doublereal isb_i__, is_kti__, ktimin, ktimax, ktiint, isb_kti__;


/*   A SIMPLE FUNCTION TO DETERMINE THE NORMALIZATION PRE-FACTOR */
/*   ON THE ANTI-QUARK WAVEFUNCTION GIVEN A PARTICULAR SET OF */
/*   PARAMETERS (I.E., Nq, Lq, Lqb, m_S, m_Sb) */

/*    THE ANALOGUE OF `Nqb_Find.f', BUT NOW INVOLVING NUMERICAL */
/*    INTEGRATIONS OVER kT; ALSO, WAVE-FNCTS FLAGGED BY 'typ' */

/*  WRITTEN: T. Hobbs (DEC 2, 2014) */
/*  MODIFIED FOR CHARM: T. Hobbs (SEPT 27, 2016) */
/* ********************************************************************** */
/* VARIABLE DECLARATIONS */
/* *********************************************************************** */
    pi = atan(1.) * 4.;
    mn = .9382720813;
/* PROTON MASS, GeV */
    mqb = *mq;
/* --------------------------------------------------------------------- */
/* *********************************************************************** */
/* .  .  .  .  DETERMINE THE RELATION BETWEEN WF NORMALIZATIONS .  .  .  . */
/* STRUCK QUARK MASSES IDENTICAL! */
    is = 0.;
    isb = 0.;
    wmin = 0.;
    wmax = .999;
/* INTEGRATE OVER x-->w FULLY! */
    wint = (wmax - wmin) / 100.;
    for (iw = 1; iw <= 100; ++iw) {
	w = wmin + wint * (doublereal) iw;
	is_kti__ = 0.;
/* AS WELL AS INTERNALLY OVER kT */
	isb_kti__ = 0.;
	ktimin = 0.;
	ktimax = 10.;
	ktiint = (ktimax - ktimin) / 100.;
	for (ikti = 1; ikti <= 100; ++ikti) {
	    kti = ktimin + ktiint * (doublereal) ikti;
	    is_i__ = f1int_(&c_b2, &w, &kti, &c_b3, lq, mq, ms, typ);
	    isb_i__ = f1int_(&c_b2, &w, &kti, &c_b3, lqb, mq, msb, typ);
	    if (ikti == 0) {
		is_kti__ += is_i__;
		isb_kti__ += isb_i__;
	    } else if (ikti / 2 << 1 != ikti) {
		is_kti__ += is_i__ * 4.;
		isb_kti__ += isb_i__ * 4.;
	    } else if (ikti / 2 << 1 == ikti) {
		is_kti__ += is_i__ * 2.;
		isb_kti__ += isb_i__ * 2.;
	    }
	}
	is_kti__ = ktiint / 3. * is_kti__;
	isb_kti__ = ktiint / 3. * isb_kti__;
/*  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . */
	if (iw == 0) {
	    is += is_kti__;
	    isb += isb_kti__;
	} else if (iw / 2 << 1 != iw) {
	    is += is_kti__ * 4.;
	    isb += isb_kti__ * 4.;
	} else if (iw / 2 << 1 == iw) {
	    is += is_kti__ * 2.;
	    isb += isb_kti__ * 2.;
	}
    }
    is = wint / 3. * is;
    isb = wint / 3. * isb;
/* ....THESE ENABLE US TO RELATE THE NORMALIZATION CONSTANTS; THEY ARE THEN: */
    ret_val = *nq * (is / isb);
/* ________________________________________________________________________________________ */
    return ret_val;
} /* nqbf_ */

