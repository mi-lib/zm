/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_pex - interpolation: polynomial function.
 */

#ifndef __ZM_IP_PEX_H__
#define __ZM_IP_PEX_H__

/* NOTE: never include this header file in user programs. */

#include <zm/zm_pex.h>

__BEGIN_DECLS

/*! \struct zPexIP
 * \brief polynomial function class.
 *
 * zPexIP is a class for the polynomial function, which is mathematically expressed as
 *   f(x) = c0 + c1*t/T + c2*(t/T)^2 + c3*(t/T)^3 + ...,
 * where \a x is the output variable and \a t is the input variable. It assumes that
 * \a x is for displacement and \a t is for time. \a T represents the total duration
 * during which \a t is defined.
 * zPexIP stores the coefficients \a c0, \a c1, ..., up to the dimension of the function.
 * \sa
 * zPex
 */
ZDEF_STRUCT( __ZM_CLASS_EXPORT, zPexIP ){
  double term;
  zPex c;
};

#define zPexIPDim(p)          zPexDim( (p)->c )

#define zPexIPTerm(p)         (p)->term
#define zPexIPSetTerm(p,t)    ( (p)->term = (t) )

#define zPexIPCoeff(p,i)      zPexCoeff((p)->c,i)
#define zPexIPSetCoeff(p,i,v) zPexSetCoeff((p)->c,i,v)

/*! \brief allocate and free a polynomial function.
 *
 * zPexIPAlloc() allocates a new polynomial function \a pc with a dimension \a dim.
 * \a term is the total duration during which the function is defined.
 *
 * zPexIPFree() frees inner parameters of the polynomial function \a pc.
 * \return
 * zPexIPAlloc() returns the true value if it succeeds to allocate internal variables.
 * Otherwise, it returns the false value.
 *
 * zPexIPFree() does not return any value.
 */
__ZM_EXPORT bool zPexIPAlloc(zPexIP *pc, double term, int dim);
__ZM_EXPORT void zPexIPFree(zPexIP *pc);

/*! \brief creates a polynomial function from the boundary conditions and a sequence of points to fit.
 *
 * zPexIPCreateBoundary() creates a polynomial function \a pc that satisfies the specified
 * boundary conditions. The conditions are defined by the total duration \a term, the initial
 * displacement \a x1, velocity \a v1, and acceleration \a a1, and the terminal displacement
 * \a x2, velocity \a v2, and acceleration \a a2. Namely,
 *                f(0) = \a x1
 *            df(0)/dt = \a v1
 *        d^2f(0)/dt^2 = \a a1
 *          f(\a term) = \a x2
 *      df(\a term)/dt = \a v2
 *  d^2f(\a term)/dt^2 = \a a2
 * \a v is for the coefficients of the third-order to n-3 order terms, where n is the
 * dimension of \a pc, which is the size of \a v +5. The null pointer is acceptable as \a v,
 * in which case the dimension will be 5. The remaining coefficients are calculated from
 * the specified boundary condition.
 *
 * zPexIPCreateLSM() creates a polynomial function \a pc so that the error from the sequence of
 * points {(t_i, x_i)} specified by vectors \a t and \a x is minimized based on the least square
 * method. \a term is the total duration of the function, and \a dim is the dimension of the
 * polynomial.
 *
 * zPexIPCreateBounderyLSM() creates a polynomial function \a pc that satisfies the boundary
 * conditions specified by \a term, \a x1, \a v1, \a a1, \a x2, \a v2, and \a a2, and also
 * fits the sequence of points {(t_i, x_i)} specified by vectors \a t and \a x. The meanings
 * of those arguments are the same with those for zPexIPCreateBoundary() and zPexIPCreateLSM().
 * \a dim is the dimension of the polynomial.
 * \return
 * zPexIPCreateBoundary(), zPexIPCreateLSM() and zPexIPCreateBounderyLSM() return the true value
 * if they succeed to create the polynomial function, or the false value otherwise.
 */
__ZM_EXPORT bool zPexIPCreateBoundary(zPexIP *pc, double term, double x1, double v1, double a1, double x2, double v2, double a2, zVec v);
__ZM_EXPORT bool zPexIPCreateLSM(zPexIP *pc, double term, int dim, zVec t, zVec x);
/*! \brief create a polynomial function from the boundary condition and a sequence of points to fit. */
__ZM_EXPORT bool zPexIPCreateBounderyLSM(zPexIP *pc, double term, double x1, double v1, double a1, double x2, double v2, double a2, int dim, zVec t, zVec x);

/*! \brief value of polynomial function.
 *
 * \return
 * zPexIPVal() returns the value of a polynomial function \a pc at the time \a t.
 */
__ZM_EXPORT double zPexIPVal(zPexIP *pc, double t);

/*! \brief velocity of a polynomial function. */
__ZM_EXPORT double zPexIPVel(zPexIP *pc, double t);

/*! \brief acceleration of a polynomial function. */
__ZM_EXPORT double zPexIPAcc(zPexIP *pc, double t);

/*! \brief print a polynomial function.
 *
 * zPexIPFPrint() prints out a polynomial function \a pc in the expression form with polynomial terms.
 *
 * zPexIPPrint() prints \a pc to the standerd output.
 * \return
 * zPexIPFPrint() and zPexIPPrint() return no value.
 */
__ZM_EXPORT void zPexIPFPrint(FILE *fp, zPexIP *pc);
#define zPexIPPrint(p) zPexIPFPrint( stdout, (p) )

__END_DECLS

#endif /* __ZM_IP_PEX_H__ */
