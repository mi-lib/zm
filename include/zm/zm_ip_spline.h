/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_spline - interpolation: spline interpolation.
 */

#ifndef __ZM_IP_SPLINE_H__
#define __ZM_IP_SPLINE_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

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

__END_DECLS

#endif /* __ZM_IP_SPLINE_H__ */
