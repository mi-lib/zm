/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_pex - interpolation: polynomial curve.
 */

#ifndef __ZM_IP_PEX_H__
#define __ZM_IP_PEX_H__

/* NOTE: never include this header file in user programs. */

#include <zm/zm_pex.h>

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zPexIP
 * polynomial curve class with designable coefficients
 * ********************************************************** */

typedef struct{
  double term;
  zPex c;
} zPexIP;

/*! \brief allocate and free polynomial curve with designable coefficients.
 *
 * zPexIP is a class for the polynomial curve, which is
 * mathematically expressed as:
 *   f(x) = c0 + c1*t/T + c2*(t/T)^2 + c3*(t/T)^3 + ...
 *
 * zPexIPAlloc() allocates a new polynomial curve \a pc as \a dim
 * dimension polynomial expression. \a term is the terminal time.
 *
 * zPexIPAllocBoundary() allocates \a pc from the boundary
 * condition. \a x1, \a v1 and \a a1 are for the initial condition
 * of position, velocity and acceleration, respectively, while
 * \a x2, \a v2 and \a a2 are for the terminal condition at the time
 * \a term.
 * \a v is for the coefficients of the third-order to n-3 order
 * terms, where n is the dimension of \a pc, which is the size of
 * \a v +5. When \a v is the null vector, it is ignored.
 * The remainder of coefficients are calculated from the boundary
 * condition.
 *
 * zPexIPCreateLSM() creates \a pc based on the least square
 * method, making it fit to a point set given as (t_i, x_i), where
 * t_i and x_i are the i th components of \a t and \a x, respectively.
 *
 * zPexIPFree() frees the inner parameters of a polynomial curve
 * \a pc.
 * \return
 * zPexIPAlloc(), zPexIPAllocBoundary() and zPexIPCreateLSM()
 * return the true value if they succeed to allocate the internal
 * vector, or the false value otherwize.
 *
 * zPexIPFree() returns no value.
 */
__EXPORT bool zPexIPAlloc(zPexIP *pc, double term, int dim);
__EXPORT void zPexIPFree(zPexIP *pc);

/*! \brief set the boundary condition of a polynomial curve. */
__EXPORT bool zPexIPBoundary(zPexIP *pc, double x1, double v1, double a1, double x2, double v2, double a2, zVec v);

/*! \brief allocate a polynomial curve from the boundary condition. */
__EXPORT bool zPexIPCreateBoundary(zPexIP *pc, double term, double x1, double v1, double a1, double x2, double v2, double a2, zVec v);

/*! \brief fit a polynomial curve to a sequence of points based on the least square method. */
__EXPORT bool zPexIPLSM(zPexIP *pc, zVec t, zVec x);

/*! \brief create a polynomial curve that fits a sequence of points based on the least square method. */
__EXPORT bool zPexIPCreateLSM(zPexIP *pc, double term, int dim, zVec t, zVec x);

/*! \brief create a polynomial curve from the boundary condition and a sequence of points to fit. */
__EXPORT bool zPexIPCreateBounderyLSM(zPexIP *pc, double term, double x1, double v1, double a1, double x2, double v2, double a2, int dim, zVec t, zVec x);

#define zPexIPDim(p)          zPexDim( (p)->c )

#define zPexIPTerm(p)         (p)->term
#define zPexIPSetTerm(p,t)    ( (p)->term = (t) )

#define zPexIPCoeff(p,i)      zPexCoeff((p)->c,i)
#define zPexIPSetCoeff(p,i,v) zPexSetCoeff((p)->c,i,v)

/*! \brief value of polynomial curve.
 *
 * zPexIPVal() calculates the value of a polynomial curve
 * \a pc at the time \a t.
 * \return
 * zPexIPVal() returns the value calculated.
 */
__EXPORT double zPexIPVal(zPexIP *pc, double t);

/*! \brief velocity of a polynomial curve. */
__EXPORT double zPexIPVel(zPexIP *pc, double t);

/*! \brief acceleration of a polynomial curve. */
__EXPORT double zPexIPAcc(zPexIP *pc, double t);

/*! \brief print a polynomial curve.
 *
 * zPexIPFPrint() prints out a polynomial curve \a pc in
 * the expression form with polynomial terms.
 *
 * zPexIPPrint() prints \a pc to the standerd output.
 * \return
 * zPexIPFPrint() and zPexIPPrint() return no value.
 */
__EXPORT void zPexIPFPrint(FILE *fp, zPexIP *pc);
#define zPexIPPrint(p) zPexIPFPrint( stdout, (p) )

__END_DECLS

#endif /* __ZM_IP_PEX_H__ */
