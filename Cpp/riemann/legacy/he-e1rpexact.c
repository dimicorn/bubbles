/* he-e1rpexact.f -- translated by f2c (version 20200916).
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

/* Common Block Declarations */

struct {
    real gamma, g1, g2, g3, g4, g5, g6, g7, g8;
} gammas_;

#define gammas_1 gammas_

struct {
    real dl, ul, pl, cl, dr, ur, pr, cr;
} states_;

#define states_1 states_

/* Table of constant values */

static integer c__4 = 4;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__9 = 9;

/* ----------------------------------------------------------* */
/*                                                          * */
/*         EXACT RIEMANN SOLVER                             * */
/*       FOR THE EULER EQUATIONS                            * */
/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_20[] = "(5(f14.6,2x))";

    /* System generated locals */
    integer i__1;
    real r__1, r__2;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(olist *), s_rsle(cilist *), do_lio(integer *, integer *, 
	    char *, ftnlen), e_rsle(void), f_clos(cllist *);
    double sqrt(doublereal);
    integer s_wsle(cilist *), e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__;
    static real s, ds, dx, pm, ps, um, us, mpa, xpos, diaph;
    static integer cells;
    static real domlen;
    extern /* Subroutine */ int sample_(real *, real *, real *, real *, real *
	    , real *), starpu_(real *, real *, real *);
    static real timeout;

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 1, 0, 0, 0 };
    static cilist io___3 = { 0, 1, 0, 0, 0 };
    static cilist io___5 = { 0, 1, 0, 0, 0 };
    static cilist io___7 = { 0, 1, 0, 0, 0 };
    static cilist io___8 = { 0, 1, 0, 0, 0 };
    static cilist io___10 = { 0, 1, 0, 0, 0 };
    static cilist io___11 = { 0, 1, 0, 0, 0 };
    static cilist io___12 = { 0, 1, 0, 0, 0 };
    static cilist io___13 = { 0, 1, 0, 0, 0 };
    static cilist io___14 = { 0, 1, 0, 0, 0 };
    static cilist io___15 = { 0, 1, 0, 0, 0 };
    static cilist io___16 = { 0, 1, 0, 0, 0 };
    static cilist io___18 = { 0, 6, 0, 0, 0 };
    static cilist io___19 = { 0, 6, 0, 0, 0 };
    static cilist io___20 = { 0, 6, 0, 0, 0 };
    static cilist io___21 = { 0, 6, 0, 0, 0 };
    static cilist io___31 = { 0, 2, 0, fmt_20, 0 };


/* Declaration of variables: */
    o__1.oerr = 0;
    o__1.ounit = 1;
    o__1.ofnmlen = 9;
    o__1.ofnm = "exact.ini";
    o__1.orl = 0;
    o__1.osta = "UNKNOWN";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
/* Initial data and parameters are read in */
    s_rsle(&io___1);
    do_lio(&c__4, &c__1, (char *)&domlen, (ftnlen)sizeof(real));
    e_rsle();
/* Domain length */
    s_rsle(&io___3);
    do_lio(&c__4, &c__1, (char *)&diaph, (ftnlen)sizeof(real));
    e_rsle();
/* Initial discontinuity position */
    s_rsle(&io___5);
    do_lio(&c__3, &c__1, (char *)&cells, (ftnlen)sizeof(integer));
    e_rsle();
/* Number of computing cells */
    s_rsle(&io___7);
    do_lio(&c__4, &c__1, (char *)&gammas_1.gamma, (ftnlen)sizeof(real));
    e_rsle();
/* Ratio of specific heats */
    s_rsle(&io___8);
    do_lio(&c__4, &c__1, (char *)&timeout, (ftnlen)sizeof(real));
    e_rsle();
/* Output time */
    s_rsle(&io___10);
    do_lio(&c__4, &c__1, (char *)&states_1.dl, (ftnlen)sizeof(real));
    e_rsle();
/* Initial density on left state */
    s_rsle(&io___11);
    do_lio(&c__4, &c__1, (char *)&states_1.ul, (ftnlen)sizeof(real));
    e_rsle();
/* Initial velocity on left state */
    s_rsle(&io___12);
    do_lio(&c__4, &c__1, (char *)&states_1.pl, (ftnlen)sizeof(real));
    e_rsle();
/* Initial pressure on left state */
    s_rsle(&io___13);
    do_lio(&c__4, &c__1, (char *)&states_1.dr, (ftnlen)sizeof(real));
    e_rsle();
/* Initial density on right state */
    s_rsle(&io___14);
    do_lio(&c__4, &c__1, (char *)&states_1.ur, (ftnlen)sizeof(real));
    e_rsle();
/* Initial velocity on right state */
    s_rsle(&io___15);
    do_lio(&c__4, &c__1, (char *)&states_1.pr, (ftnlen)sizeof(real));
    e_rsle();
/* Initial pressure on right state */
    s_rsle(&io___16);
    do_lio(&c__4, &c__1, (char *)&mpa, (ftnlen)sizeof(real));
    e_rsle();
/* Normalising constant */
    cl__1.cerr = 0;
    cl__1.cunit = 1;
    cl__1.csta = 0;
    f_clos(&cl__1);
/* Compute gamma related constants */
    gammas_1.g1 = (gammas_1.gamma - 1.f) / (gammas_1.gamma * 2.f);
    gammas_1.g2 = (gammas_1.gamma + 1.f) / (gammas_1.gamma * 2.f);
    gammas_1.g3 = gammas_1.gamma * 2.f / (gammas_1.gamma - 1.f);
    gammas_1.g4 = 2.f / (gammas_1.gamma - 1.f);
    gammas_1.g5 = 2.f / (gammas_1.gamma + 1.f);
    gammas_1.g6 = (gammas_1.gamma - 1.f) / (gammas_1.gamma + 1.f);
    gammas_1.g7 = (gammas_1.gamma - 1.f) / 2.f;
    gammas_1.g8 = gammas_1.gamma - 1.f;
/* Compute sound speeds */
    states_1.cl = sqrt(gammas_1.gamma * states_1.pl / states_1.dl);
    states_1.cr = sqrt(gammas_1.gamma * states_1.pr / states_1.dr);
/* The pressure positivity condition is tested for */
    if (gammas_1.g4 * (states_1.cl + states_1.cr) <= states_1.ur - 
	    states_1.ul) {
/* The initial data is such that vaccum is generated. */
/* Program stopped */
	s_wsle(&io___18);
	e_wsle();
	s_wsle(&io___19);
	do_lio(&c__9, &c__1, "***Vacuum is generated by data***", (ftnlen)33);
	e_wsle();
	s_wsle(&io___20);
	do_lio(&c__9, &c__1, "***Program stopped***", (ftnlen)21);
	e_wsle();
	s_wsle(&io___21);
	e_wsle();
	s_stop("", (ftnlen)0);
    }
/* Exact solution for pressure and velocity in star region is found */
    starpu_(&pm, &um, &mpa);
    dx = domlen / (real) cells;
/* Complete solution at time TIMOUT is found */
    o__1.oerr = 0;
    o__1.ounit = 2;
    o__1.ofnmlen = 9;
    o__1.ofnm = "test.out";
    o__1.orl = 0;
    o__1.osta = "UNKNOWN";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    i__1 = cells;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xpos = ((real) i__ - .5f) * dx;
	s = (xpos - diaph) / timeout;
/* Solution at point (X, T) = (XPOS - DIAPH, TIMEOUT) is found */
	sample_(&pm, &um, &s, &ds, &us, &ps);
/* Exact solution profiles are written to exact.out */
	s_wsfe(&io___31);
	do_fio(&c__1, (char *)&xpos, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ds, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&us, (ftnlen)sizeof(real));
	r__1 = ps / mpa;
	do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	r__2 = ps / ds / gammas_1.g8 / mpa;
	do_fio(&c__1, (char *)&r__2, (ftnlen)sizeof(real));
	e_wsfe();
/* L10: */
    }
    cl__1.cerr = 0;
    cl__1.cunit = 2;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
} /* MAIN__ */

/* Subroutine */ int starpu_(real *p, real *u, real *mpa)
{
    /* Initialized data */

    static real tolpre = 1e-6f;
    static integer nriter = 20;

    /* Format strings */
    static char fmt_30[] = "(5x,i5,15x,f12.7)";
    static char fmt_40[] = "(2(f14.6,5x))";

    /* System generated locals */
    integer i__1;
    real r__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void), s_wsfe(cilist *), do_fio(integer *, char *, ftnlen),
	     e_wsfe(void);

    /* Local variables */
    static integer i__;
    static real fl, fr, fld, frd, pold, udiff, change;
    extern /* Subroutine */ int prefun_(real *, real *, real *, real *, real *
	    , real *), guessp_(real *);
    static real pstart;

    /* Fortran I/O blocks */
    static cilist io___37 = { 0, 6, 0, 0, 0 };
    static cilist io___38 = { 0, 6, 0, 0, 0 };
    static cilist io___39 = { 0, 6, 0, 0, 0 };
    static cilist io___46 = { 0, 6, 0, fmt_30, 0 };
    static cilist io___47 = { 0, 6, 0, 0, 0 };
    static cilist io___48 = { 0, 6, 0, 0, 0 };
    static cilist io___49 = { 0, 6, 0, 0, 0 };
    static cilist io___50 = { 0, 6, 0, 0, 0 };
    static cilist io___51 = { 0, 6, 0, fmt_40, 0 };
    static cilist io___52 = { 0, 6, 0, 0, 0 };


/* Declaration of variables */
    guessp_(&pstart);
    pold = pstart;
    udiff = states_1.ur - states_1.ul;
    s_wsle(&io___37);
    do_lio(&c__9, &c__1, "----------------------------------------", (ftnlen)
	    40);
    e_wsle();
    s_wsle(&io___38);
    do_lio(&c__9, &c__1, "   Iteration number      Change  ", (ftnlen)33);
    e_wsle();
    s_wsle(&io___39);
    do_lio(&c__9, &c__1, "----------------------------------------", (ftnlen)
	    40);
    e_wsle();
    i__1 = nriter;
    for (i__ = 1; i__ <= i__1; ++i__) {
	prefun_(&fl, &fld, &pold, &states_1.dl, &states_1.pl, &states_1.cl);
	prefun_(&fr, &frd, &pold, &states_1.dr, &states_1.pr, &states_1.cr);
	*p = pold - (fl + fr + udiff) / (fld + frd);
	change = (r__1 = (*p - pold) / (*p + pold), dabs(r__1)) * 2.f;
	s_wsfe(&io___46);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&change, (ftnlen)sizeof(real));
	e_wsfe();
	if (change <= tolpre) {
	    goto L20;
	}
	if (*p < 0.f) {
	    *p = tolpre;
	}
	pold = *p;
/* L10: */
    }
    s_wsle(&io___47);
    do_lio(&c__9, &c__1, "Divergence in Newton-Raphson iteration", (ftnlen)38)
	    ;
    e_wsle();
L20:
/* Compute velocity in Star Region */
    *u = (states_1.ul + states_1.ur + fr - fl) * .5f;
    s_wsle(&io___48);
    do_lio(&c__9, &c__1, "---------------------------------------", (ftnlen)
	    39);
    e_wsle();
    s_wsle(&io___49);
    do_lio(&c__9, &c__1, "   Pressure        Velocity", (ftnlen)27);
    e_wsle();
    s_wsle(&io___50);
    do_lio(&c__9, &c__1, "---------------------------------------", (ftnlen)
	    39);
    e_wsle();
    s_wsfe(&io___51);
    r__1 = *p / *mpa;
    do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&(*u), (ftnlen)sizeof(real));
    e_wsfe();
    s_wsle(&io___52);
    do_lio(&c__9, &c__1, "---------------------------------------", (ftnlen)
	    39);
    e_wsle();
    return 0;
} /* starpu_ */

/* Subroutine */ int guessp_(real *pm)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static real pq, um, gel, ger, cup, ptl, ppv, ptr, pmin, pmax, qmax, quser;

/* Declaration of variables */
    quser = 2.f;
/* Compute guess pressure from PVRS Riemann solver */
    cup = (states_1.dl + states_1.dr) * .25f * (states_1.cl + states_1.cr);
    ppv = (states_1.pl + states_1.pr) * .5f + (states_1.ul - states_1.ur) * 
	    .5f * cup;
    ppv = dmax(0.f,ppv);
    pmin = dmin(states_1.pl,states_1.pr);
    pmax = dmax(states_1.pl,states_1.pr);
    qmax = pmax / pmin;
    if (qmax <= quser && (pmin <= ppv && ppv <= pmax)) {
/* Select PVRS Riemann solver */
	*pm = ppv;
    } else {
	if (ppv < pmin) {
/* Select Two-Rarefaction Riemann solver */
	    d__1 = (doublereal) (states_1.pl / states_1.pr);
	    d__2 = (doublereal) gammas_1.g1;
	    pq = pow_dd(&d__1, &d__2);
	    um = (pq * states_1.ul / states_1.cl + states_1.ur / states_1.cr 
		    + gammas_1.g4 * (pq - 1.f)) / (pq / states_1.cl + 1.f / 
		    states_1.cr);
	    ptl = gammas_1.g7 * (states_1.ul - um) / states_1.cl + 1.f;
	    ptr = gammas_1.g7 * (um - states_1.ur) / states_1.cr + 1.f;
	    d__1 = (doublereal) ptl;
	    d__2 = (doublereal) gammas_1.g3;
	    d__3 = (doublereal) ptr;
	    d__4 = (doublereal) gammas_1.g3;
	    *pm = (states_1.pl * pow_dd(&d__1, &d__2) + states_1.pr * pow_dd(&
		    d__3, &d__4)) * .5f;
	} else {
/* Select Two-Shock Riemann solver with PVRS as estimate */
	    gel = sqrt(gammas_1.g5 / states_1.dl / (gammas_1.g6 * states_1.pl 
		    + ppv));
	    ger = sqrt(gammas_1.g5 / states_1.dr / (gammas_1.g6 * states_1.pr 
		    + ppv));
	    *pm = (gel * states_1.pl + ger * states_1.pr - (states_1.ur - 
		    states_1.ul)) / (gel + ger);
	}
    }
    return 0;
} /* guessp_ */

/* Subroutine */ int prefun_(real *f, real *fd, real *p, real *dk, real *pk, 
	real *ck)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static real ak, bk, qrt, prat;

/* Declaration of variables */
    if (*p <= *pk) {
/* Rarefaction wave */
	prat = *p / *pk;
	d__1 = (doublereal) prat;
	d__2 = (doublereal) gammas_1.g1;
	*f = gammas_1.g4 * *ck * (pow_dd(&d__1, &d__2) - 1.f);
	d__1 = (doublereal) prat;
	d__2 = (doublereal) (-gammas_1.g2);
	*fd = 1.f / (*dk * *ck) * pow_dd(&d__1, &d__2);
    } else {
/* Shock wave */
	ak = gammas_1.g5 / *dk;
	bk = gammas_1.g6 * *pk;
	qrt = sqrt(ak / (bk + *p));
	*f = (*p - *pk) * qrt;
	*fd = (1.f - (*p - *pk) * .5f / (bk + *p)) * qrt;
    }
    return 0;
} /* prefun_ */

/* Subroutine */ int sample_(real *pm, real *um, real *s, real *d__, real *u, 
	real *p)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static real c__, sl, sr, cml, cmr, shl, pml, shr, pmr, stl, str;

/* Declaration of variables */
    if (*s <= *um) {
	if (*pm <= states_1.pl) {
	    shl = states_1.ul - states_1.cl;
	    if (*s <= shl) {
		*d__ = states_1.dl;
		*u = states_1.ul;
		*p = states_1.pl;
	    } else {
		d__1 = (doublereal) (*pm / states_1.pl);
		d__2 = (doublereal) gammas_1.g1;
		cml = states_1.cl * pow_dd(&d__1, &d__2);
		stl = *um - cml;
		if (*s > stl) {
		    d__1 = (doublereal) (*pm / states_1.pl);
		    d__2 = (doublereal) (1.f / gammas_1.gamma);
		    *d__ = states_1.dl * pow_dd(&d__1, &d__2);
		    *u = *um;
		    *p = *pm;
		} else {
		    *u = gammas_1.g5 * (states_1.cl + gammas_1.g7 * 
			    states_1.ul + *s);
		    c__ = gammas_1.g5 * (states_1.cl + gammas_1.g7 * (
			    states_1.ul - *s));
		    d__1 = (doublereal) (c__ / states_1.cl);
		    d__2 = (doublereal) gammas_1.g4;
		    *d__ = states_1.dl * pow_dd(&d__1, &d__2);
		    d__1 = (doublereal) (c__ / states_1.cl);
		    d__2 = (doublereal) gammas_1.g3;
		    *p = states_1.pl * pow_dd(&d__1, &d__2);
		}
	    }
	} else {
/* Left shock */
	    pml = *pm / states_1.pl;
	    sl = states_1.ul - states_1.cl * sqrt(gammas_1.g2 * pml + 
		    gammas_1.g1);
	    if (*s <= sl) {
		*d__ = states_1.dl;
		*u = states_1.ul;
		*p = states_1.pl;
	    } else {
		*d__ = states_1.dl * (pml + gammas_1.g6) / (pml * gammas_1.g6 
			+ 1.f);
		*u = *um;
		*p = *pm;
	    }
	}
    } else {
	if (*pm > states_1.pr) {
/* Right shock */
	    pmr = *pm / states_1.pr;
	    sr = states_1.ur + states_1.cr * sqrt(gammas_1.g2 * pmr + 
		    gammas_1.g1);
	    if (*s >= sr) {
		*d__ = states_1.dr;
		*u = states_1.ur;
		*p = states_1.pr;
	    } else {
		*d__ = states_1.dr * (pmr + gammas_1.g6) / (pmr * gammas_1.g6 
			+ 1.f);
		*u = *um;
		*p = *pm;
	    }
	} else {
/* Right rarefaction */
	    shr = states_1.ur + states_1.cr;
	    if (*s >= shr) {
		*d__ = states_1.dr;
		*u = states_1.ur;
		*p = states_1.pr;
	    } else {
		d__1 = (doublereal) (*pm / states_1.pr);
		d__2 = (doublereal) gammas_1.g1;
		cmr = states_1.cr * pow_dd(&d__1, &d__2);
		str = *um + cmr;
		if (*s <= str) {
		    d__1 = (doublereal) (*pm / states_1.pr);
		    d__2 = (doublereal) (1.f / gammas_1.gamma);
		    *d__ = states_1.dr * pow_dd(&d__1, &d__2);
		    *u = *um;
		    *p = *pm;
		} else {
		    *u = gammas_1.g5 * (-states_1.cr + gammas_1.g7 * 
			    states_1.ur + *s);
		    c__ = gammas_1.g5 * (states_1.cr - gammas_1.g7 * (
			    states_1.ur - *s));
		    d__1 = (doublereal) (c__ / states_1.cr);
		    d__2 = (doublereal) gammas_1.g4;
		    *d__ = states_1.dr * pow_dd(&d__1, &d__2);
		    d__1 = (doublereal) (c__ / states_1.cr);
		    d__2 = (doublereal) gammas_1.g3;
		    *p = states_1.pr * pow_dd(&d__1, &d__2);
		}
	    }
	}
    }
    return 0;
} /* sample_ */

