/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_spline - interpolation: spline interpolation.
 */

#ifndef __ZM_IP_SPLINE_H__
#define __ZM_IP_SPLINE_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief value on Ferguson curve.
 *
 * zFergusonVal() returns a value on a Ferguson curve, which connects two boundary values \a x0 and \a x1
 * with specified gradients \a v0 and \a v1, respectively, by a cubic Hermite functions. \a term is the
 * interval between the boundaries. \a t is the current time, which should be between 0 and \a term.
 * \return
 * zFergusonVal() returns the computed value on a Ferguson curve corresponding to \a t/\a term.
 */
__ZM_EXPORT double zFergusonVal(double t, double term, double x0, double v0, double x1, double v1);

/* spline interpolation, which expresses each segment divided by the sections for a third-order polynomial
 * curve with keeping the continuity of velocity at every sections.
 *
 * For the description of \a etype1, \a etype2, \a v1 and \a v2, see zm_ip.h.
 */
enum{ ZSPLINE_INVALID, ZSPLINE_FIX_EDGE, ZSPLINE_FREE_EDGE };

/* \brief create spline interpolator.
 *
 * zIPCreateSpline() creates a 3rd-order spline interpolator \a ip. \a etype1 and \a etype2 are the types
 * of edges at the beginning point and the termination point, respectively. They should be chosen from
 * ZSPLINE_FIX_EDGE and ZSPLINE_FREE_EDGE, which are for fixed-edge and free-edge, respectively.
 * One can set the velocity at each endpoint for \a v1 (or \a v2) by choosing ZSPLINE_FIX_EDGE. When
 * ZSPLINE_FREE_EDGE is specified, \a v1 (or \a v2) is ignored.
 * \return
 * zIPCreateSpline() returns a pointer \a ip when it succeeds to create the interpolator. Otherwise, it
 * returns the null pointer.
 */
__ZM_EXPORT bool zIPCreateSpline(zIP *ip, const zSeq *seq, int etype1, const zVec v1, int etype2, const zVec v2);

/*! \brief find coefficients of a cubic polynomial of a segment of a spline interpolation.
 *
 * zIPSplineCoeff() finds coefficients of a cubic polynomial of the \a i th segment of a spline interpolator
 * \a ip. Namely, \a a, \a b, \a c, and \a d corresponds to the following expression:
 * x = \a a (t-t_\a i)^3 + \a b (t-t_\a i)^2 + \a c (t-t_\a i) + \a d
 * \return
 * zIPSplineCoeff() returns the false value if the sizes of \a a, \a b, \a c, \a d, and internal vectors
 * of \a ip are not equal. Otherwise, it returns the true value.
 */
__ZM_EXPORT bool zIPSplineCoeff(const zIP *ip, int i, zVec a, zVec b, zVec c, zVec d);

/*! \brief create a Piecewise Cubic Hermite Interporating Polynomial interpolator.
 *
 * zIPCreatePCHIP() creates a Piecewise Cubic Hermite Interporating Polynomial
 * interpolator. It connects a series of monotonously increasing/decreasing points
 * by monotonously increasing/decreasing continuous curve.
 * The algorithm to create the interpolator was proposed by
 * F. N. Fritsch and R. E. Carlson in 1980.
 *
 * \a ip is a pointer to an instance interpolator.
 * \a seq is a sequence of points to be interpolated.
 */
__ZM_EXPORT bool zIPCreatePCHIP(zIP *ip, const zSeq *seq);

__END_DECLS

#endif /* __ZM_IP_SPLINE_H__ */
