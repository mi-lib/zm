/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_lagrange - interpolation: Lagrange's interpolation.
 */

#ifndef __ZM_IP_LAGRANGE_H__
#define __ZM_IP_LAGRANGE_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* \brief create Lagrange interpolator.
 *
 * zIPCreateLagrange() creates a Lagrange interpolator \a ip.
 * It connects given n points by n-1 th order polynomial curve.
 * \a seq is a sequence of points to be interpolated.
 * \return
 * zIPCreateLagrange() returns a pointer \a ip when it succeeds to create the interpolator. Otherwise, it
 * returns the null pointer.
 */
__ZM_EXPORT bool zIPCreateLagrange(zIP *ip, const zSeq *seq);

__END_DECLS

#endif /* __ZM_IP_LAGRANGE_H__ */
