/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_pchip - interpolation: Piecewise Cubic Hermite Interporating Polynomial.
 */

#ifndef __ZM_IP_PCHIP_H__
#define __ZM_IP_PCHIP_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* Piecewise Cubic Hermite Interporating Polynomial, which connects
 * a series of monotonously increasing/decreasing points by monotonously
 * increasing/decreasing continuous curve.
 * The algorithm to create the interpolator was proposed by
 * F. N. Fritsch and R. E. Carlson in 1980.
 */
__EXPORT bool zIPCreatePCHIP(zIP *ip, zSeq *seq);

__END_DECLS

#endif /* __ZM_IP_PCHIP_H__ */
