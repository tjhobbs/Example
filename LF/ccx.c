/* ccx.f -- translated by f2c (version 20100827).
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

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__5 = 5;

/* ********************************************************************** */
/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    doublereal d__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    double atan(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void), f_open(olist *), f_clos(cllist *);

    /* Local variables */
    static doublereal c__[1000], x, c0, cb[1000], cc[1000], dc[1000], q20, pi,
	     lq, mq, kt, ms;
    static integer ix;
    static doublereal nq, cc0, dc0, f2c[1000], eq2, c_i__, lqb, mqb, nqb, msb,
	     xdc[1000];
    static integer ikt;
    static doublereal xpc[1000], xpt[1000];
    static integer typ;
    static doublereal dxc0, pxc0, cb_i__, cc_i__;
    extern doublereal nqbf_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    static doublereal xmin, xmax, xint;
    extern doublereal f1int_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    , ccint_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    static doublereal ktmin, ktmax, ktint;

    /* Fortran I/O blocks */
    static cilist io___39 = { 0, 6, 0, 0, 0 };
    static cilist io___40 = { 0, 6, 0, 0, 0 };
    static cilist io___41 = { 0, 6, 0, 0, 0 };
    static cilist io___42 = { 0, 6, 0, 0, 0 };
    static cilist io___43 = { 0, 6, 0, 0, 0 };
    static cilist io___44 = { 0, 11, 0, 0, 0 };
    static cilist io___45 = { 0, 12, 0, 0, 0 };
    static cilist io___46 = { 0, 13, 0, 0, 0 };



/*  A SIMPLE CODE TO COMPUTE/PLOT THE x DEPENDENCE OF THE FORM FACTOR */
/*  INTEGRANDS; ESSENTIALLY, AT Q2=0, THE UNINTEGRATED F1 = s(x) - sbar(x), */
/*  ETC. */

/*  WE INTEGRATE NUMERICALLY OVER kT2; MOREOVER, */
/*  THE CODE CALLS THE EXTERNAL FUNCTION F1int CONTAINED in `GEGM_5P.f' */

/*  WRITTEN: T. Hobbs (SEPT 27, 2016) */
/* ********************************************************************** */
/* VARIABLE DECLARATIONS */
/* *********************************************************************** */
    pi = atan(1.) * 4.;
    eq2 = .44444444444444442;
/* C QUARK CHARGE SQUARED */
    q20 = 0.;
/* COMPUTE AT Q2=0 */
    typ = 1;
/* USE THE STAND. GAUSSIAN WF */
    mq = 1.3;
/* CHARM MASS */
    mqb = mq;
/* __________________________________________________________________ */
/* _ _ _ _ VALUES OF THE FITTED PARAMETERS WE PLACE HERE _ _ _ _ _ _ _ _ */
/* STRUCK QUARK MASSES IDENTICAL! */
    nq = 2.9312962778399867;
    lq = 3.;
/* PARAMETER GUESSES; INITIALLY, FROM 2014 PRD */
    lqb = 3.;
    ms = 2.865;
/* APPROX. D MESON MASS */
    msb = 3.;
/* DETERMINE THE ANTI-QUARK NORM. CONSTANT FROM THE GPDs CONDITION. */
/* FROM EBERT ET AL., (2011) */
    nqb = nqbf_(&nq, &mq, &lq, &lqb, &ms, &msb, &typ);
/* __________________________________________________________________ */
/* --------------------------------------------------------------------- */
/* *********************************************************************** */
/* .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . */
/* .................PROCEED WITH THE PDF CALCULATION............. */
    c0 = 0.;
    dc0 = 0.;
    dxc0 = 0.;
    pxc0 = 0.;
    cc0 = 0.;
/* DEFINE THE BOUNDS IN x: */
    xmin = 0.;
    xmax = .999;
/* INTEGRATE OVER x FULLY! */
    xint = (xmax - xmin) / 1e3;
    for (ix = 1; ix <= 1000; ++ix) {
	x = xmin + xint * (doublereal) ix;
	xpt[ix - 1] = x;
/* DEFINE THE BOUNDS OF THE kT INTEGRAL --- i.e., kT \in [0, \infty): */
	c__[ix - 1] = 0.;
	cb[ix - 1] = 0.;
	cc[ix - 1] = 0.;
	ktmin = 0.;
	ktmax = 10.;
	ktint = (ktmax - ktmin) / 1e3;
	for (ikt = 1; ikt <= 1000; ++ikt) {
	    kt = ktmin + ktint * (doublereal) ikt;
/* ________________________________________________________________________________________ */
/* .... FOR c(x): */
	    c_i__ = f1int_(&q20, &x, &kt, &nq, &lq, &mq, &ms, &typ);
/* .... FOR cb(x): */
	    cb_i__ = f1int_(&q20, &x, &kt, &nqb, &lqb, &mqb, &msb, &typ);
/* .... FOR <cc>: */
	    cc_i__ = ccint_(&q20, &x, &kt, &nq, &lq, &mq, &ms, &typ) + ccint_(
		    &q20, &x, &kt, &nqb, &lqb, &mqb, &msb, &typ);
/* ________________________________________________________________________________________ */
/*  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . */
/* WE ADD THE EVALUATED INTEGRAND TO THE CUMULANT IN A STANDARD RUNGE-KUTTA METHOD */
/* ----- FIRST FOR THE kT INTEGRAL: */
	    if (ikt == 0) {
		c__[ix - 1] += c_i__;
		cb[ix - 1] += cb_i__;
		cc[ix - 1] += cc_i__;
	    } else if (ikt / 2 << 1 != ikt) {
		c__[ix - 1] += c_i__ * 4.;
		cb[ix - 1] += cb_i__ * 4.;
		cc[ix - 1] += cc_i__ * 4.;
	    } else if (ikt / 2 << 1 == ikt) {
		c__[ix - 1] += c_i__ * 2.;
		cb[ix - 1] += cb_i__ * 2.;
		cc[ix - 1] += cc_i__ * 2.;
	    }
	}
	c__[ix - 1] = ktint / 3. * c__[ix - 1];
	cb[ix - 1] = ktint / 3. * cb[ix - 1];
	cc[ix - 1] = ktint / 3. * cc[ix - 1];
/*  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . */
	dc[ix - 1] = c__[ix - 1] - cb[ix - 1];
	xdc[ix - 1] = x * dc[ix - 1];
	xpc[ix - 1] = x * (c__[ix - 1] + cb[ix - 1]);
	f2c[ix - 1] = eq2 * x * (c__[ix - 1] + cb[ix - 1]);
/*  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . */
/* ----- AND THEN DETERMINE NUMERICAL VALUES FOR THE x MOMENTS ALSO: */
	if (ix == 0) {
	    c0 += c__[ix - 1];
	    dc0 += dc[ix - 1];
	    dxc0 += xdc[ix - 1];
	    pxc0 += xpc[ix - 1];
	    cc0 += cc[ix - 1];
	} else if (ix / 2 << 1 != ix) {
	    c0 += c__[ix - 1] * 4.;
	    dc0 += dc[ix - 1] * 4.;
	    dxc0 += xdc[ix - 1] * 4.;
	    pxc0 += xpc[ix - 1] * 4.;
	    cc0 += cc[ix - 1] * 4.;
	} else if (ix / 2 << 1 == ix) {
	    c0 += c__[ix - 1] * 2.;
	    dc0 += dc[ix - 1] * 2.;
	    dxc0 += xdc[ix - 1] * 2.;
	    pxc0 += xpc[ix - 1] * 2.;
	    cc0 += cc[ix - 1] * 2.;
	}
    }
    c0 = xint / 3. * c0;
    dc0 = xint / 3. * dc0;
    dxc0 = xint / 3. * dxc0;
    pxc0 = xint / 3. * pxc0;
    cc0 = xint / 3. * cc0;
    s_wsle(&io___39);
    do_lio(&c__9, &c__1, "THE TOTAL CHARM PROB.:", (ftnlen)22);
    do_lio(&c__5, &c__1, (char *)&c0, (ftnlen)sizeof(doublereal));
    e_wsle();
    s_wsle(&io___40);
    do_lio(&c__9, &c__1, "THE FIRST MOMENT, x*{c+cbar}:", (ftnlen)29);
    do_lio(&c__5, &c__1, (char *)&pxc0, (ftnlen)sizeof(doublereal));
    e_wsle();
    s_wsle(&io___41);
    do_lio(&c__9, &c__1, "THE ZEROTH MOMENT, c-cbar:", (ftnlen)26);
    do_lio(&c__5, &c__1, (char *)&dc0, (ftnlen)sizeof(doublereal));
    e_wsle();
    s_wsle(&io___42);
    do_lio(&c__9, &c__1, "THE FIRST MOMENT, x*{c-cbar}:", (ftnlen)29);
    do_lio(&c__5, &c__1, (char *)&dxc0, (ftnlen)sizeof(doublereal));
    e_wsle();
    s_wsle(&io___43);
    do_lio(&c__9, &c__1, "THE CHARM SIGMA TERM, m_c <cc>, MeV:", (ftnlen)36);
    d__1 = mq * cc0 * 1e3;
    do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsle();
/* _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _ */
/* ...WRITE DATA TO FILE */
/*  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . */
    o__1.oerr = 0;
    o__1.ounit = 11;
    o__1.ofnmlen = 22;
    o__1.ofnm = "DATA/alt_c_cb_0.5%.dat";
    o__1.orl = 0;
    o__1.osta = "UNKNOWN";
    o__1.oacc = 0;
    o__1.ofm = "FORMATTED";
    o__1.oblnk = 0;
    f_open(&o__1);
    for (ix = 1; ix <= 1000; ++ix) {
	s_wsle(&io___44);
	do_lio(&c__5, &c__1, (char *)&xpt[ix - 1], (ftnlen)sizeof(doublereal))
		;
	do_lio(&c__5, &c__1, (char *)&c__[ix - 1], (ftnlen)sizeof(doublereal))
		;
	do_lio(&c__5, &c__1, (char *)&cb[ix - 1], (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&f2c[ix - 1], (ftnlen)sizeof(doublereal))
		;
	e_wsle();
    }
    cl__1.cerr = 0;
    cl__1.cunit = 11;
    cl__1.csta = 0;
    f_clos(&cl__1);
    o__1.oerr = 0;
    o__1.ounit = 12;
    o__1.ofnmlen = 22;
    o__1.ofnm = "DATA/alt_c-cb_0.5%.dat";
    o__1.orl = 0;
    o__1.osta = "UNKNOWN";
    o__1.oacc = 0;
    o__1.ofm = "FORMATTED";
    o__1.oblnk = 0;
    f_open(&o__1);
    for (ix = 1; ix <= 1000; ++ix) {
	s_wsle(&io___45);
	do_lio(&c__5, &c__1, (char *)&xpt[ix - 1], (ftnlen)sizeof(doublereal))
		;
	do_lio(&c__5, &c__1, (char *)&xdc[ix - 1], (ftnlen)sizeof(doublereal))
		;
	e_wsle();
    }
    cl__1.cerr = 0;
    cl__1.cunit = 12;
    cl__1.csta = 0;
    f_clos(&cl__1);
    o__1.oerr = 0;
    o__1.ounit = 13;
    o__1.ofnmlen = 23;
    o__1.ofnm = "DATA/alt_sigma_0.5%.dat";
    o__1.orl = 0;
    o__1.osta = "UNKNOWN";
    o__1.oacc = 0;
    o__1.ofm = "FORMATTED";
    o__1.oblnk = 0;
    f_open(&o__1);
    for (ix = 1; ix <= 1000; ++ix) {
	s_wsle(&io___46);
	do_lio(&c__5, &c__1, (char *)&xpt[ix - 1], (ftnlen)sizeof(doublereal))
		;
	do_lio(&c__5, &c__1, (char *)&cc[ix - 1], (ftnlen)sizeof(doublereal));
	e_wsle();
    }
    cl__1.cerr = 0;
    cl__1.cunit = 13;
    cl__1.csta = 0;
    f_clos(&cl__1);
/* ________________________________________________________________________________________ */
    return 0;
} /* MAIN__ */

/* Main program alias */ int ccx_ () { MAIN__ (); return 0; }
