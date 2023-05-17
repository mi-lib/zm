/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_pchip - interpolation: Piecewise Cubic Hermite Interporating Polynomial.
 */

#ifndef __ZM_IP_PCHIP_H__
#define __ZM_IP_PCHIP_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief value on Ferguson curve.
 *
 * zFergusonVal() returns a value on a Ferguson curve, which connects two
 * boundary values \a x0 and \a x1 with specified gradients \a v0 and \a v1,
 * respectively, by a cubic Hermite functions. \a term is the interval
 * between the boundaries. \a t is the current time, which should be between
 * 0 and \a term
 */
__ZM_EXPORT double zFergusonVal(double t, double term, double x0, double v0, double x1, double v1);

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
__ZM_EXPORT bool zIPCreatePCHIP(zIP *ip, zSeq *seq);

__END_DECLS

#endif /* __ZM_IP_PCHIP_H__ */
