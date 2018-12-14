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
 * 'zPexIP' is a class for the polynomial curve, which is
 * mathematically expressed as:
 *   f(x) = c0 + c1*t/T + c2*(t/T)^2 + c3*(t/T)^3 + ...
 *
 * 'zPexIPAlloc()' creates the polynomial curve 'pc' as 'dim'
 * dimension polynomial expression. 'term' is the terminal time.
 *
 * 'zPexIPAllocBoundary()' creates 'pc' from the boundary
 * condition. 'x1', 'v1' and 'a1' are the initial condition for
 * position, velocity and acceleration, respectively, while
 * 'x2', 'v2' and 'a2' are for the terminal condition at the time
 * 'term'.
 * 'v' is for the coefficients of the third-order to n-3 order
 * terms, where n is the dimension of 'pc', which is the size of
 * 'v' plus 5. When 'v' is the null vector, it is ignored.
 * The remainder of coefficients are calculated from the boundary
 * condition.
 *
 * 'zPexIPCreateLSM()' creates 'pc' based on the least square
 * method, making it fit to a point set given as (t_i, x_i), where
 * t_i and x_i are the i th components of 't' and 'x', respectively.
 *
 * 'zPexIPFree()' frees the inner parameters the polynomial curve
 * \a pc.
 * \return
 * 'zPexIPAlloc()', 'zPexIPAllocBoundary()' and
 * 'zPexIPCreateLSM()' return the true value if they succeed
 * to allocate the internal vector, or the false value, otherwize.
 *
 * zPexIPFree() returns no value.
 */
__EXPORT bool zPexIPAlloc(zPexIP *pc, double term, int dim);
__EXPORT bool zPexIPAllocBoundary(zPexIP *pc, double term, double x1, double v1, double a1, double x2, double v2, double a2, zVec v);
__EXPORT bool zPexIPCreateLSM(zPexIP *pc, double term, int dim, zVec t, zVec x);
__EXPORT void zPexIPFree(zPexIP *pc);

#define zPexIPDim(p)          zPexDim( (p)->c )

#define zPexIPTerm(p)         (p)->term
#define zPexIPSetTerm(p,t)    ( (p)->term = (t) )

#define zPexIPCoeff(p,i)      zPexCoeff((p)->c,i)
#define zPexIPSetCoeff(p,i,v) zPexSetCoeff((p)->c,i,v)

/*! \brief value of polynomial curve.
 *
 * zPexIPVal() calculates the value of the polynomial curve
 * \a pc at the time \a t.
 * \return
 * zPexIPVal() returns the value calculated.
 */
__EXPORT double zPexIPVal(zPexIP *pc, double t);

/*! \brief output of polynomial curve.
 *
 * zPexIPFWrite() prints out the polynomial curve \a pc in
 * the expression form with polynomial terms. zPexIPWrite()
 * outputs \a pc to the standerd output.
 * \return
 * zPexIPFWrite() and zPexIPWrite() return no value.
 */
__EXPORT void zPexIPFWrite(FILE *fp, zPexIP *pc);
#define zPexIPWrite(p) zPexIPFWrite( stdout, (p) )

__END_DECLS

#endif /* __ZM_IP_PEX_H__ */
