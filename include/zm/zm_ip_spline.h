/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_spline - interpolation: spline interpolation.
 */

#ifndef __ZM_IP_SPLINE_H__
#define __ZM_IP_SPLINE_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* spline interpolation, which expresses each segment
 * divided by the sections for a third-order polynomial
 * curve with keeping the continuity of velocity at every
 * sections.
 * #
 * For the description of 'etype1', 'etype2', 'v1' and
 * 'v2, see 'zm_ip.h'.
 */
enum{ ZSPLINE_INVALID, ZSPLINE_FIX_EDGE, ZSPLINE_FREE_EDGE };

/* METHOD:
 * zIPCreateSpline - create spline interpolator interpolator.
 *
 * zIPCreateSpring() creates an interpolator \a ip.
 * \a etype1 and \a etype2 are the types of edges at the
 * beginning point and the termination point, respectively.
 * They should be chosen from ZSPLINE_FIX_EDGE and
 * ZSPLINE_FREE_EDGE, which are for fixed-edge and
 * free-edge, respectively.
 * One can set the velocity at each endpoint for \a v1
 * (or \a v2) with choosing ZSPLINE_FIX_EDGE.
 * When ZSPLINE_FREE_EDGE is chosen, \a v1 (or \a v2) is ignored.
 * [RETURN VALUE]
 * 'zIPCreate()' returns a pointer 'ip' when it succeeds
 * to create interpolator, or the null pointer otherwise.
 */
__EXPORT bool zIPCreateSpline(zIP *ip, zSeq *seq, int etype1, zVec v1, int etype2, zVec v2);

__END_DECLS

#endif /* __ZM_IP_SPLINE_H__ */
