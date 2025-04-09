/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_lin - interpolation: linear interpolation.
 */

#ifndef __ZM_IP_LIN_H__
#define __ZM_IP_LIN_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* \brief linear interpolator.
 *
 * zIPCreateLinear() creates a linear interpolator \a ip.
 * It connects given n points by segmented lines.
 * \a seq is a sequence of points to be interpolated.
 * \return
 * zIPCreateLinear() returns a pointer \a ip when it succeeds to create the interpolator. Otherwise, it
 * returns the null pointer.
 */
__ZM_EXPORT bool zIPCreateLinear(zIP *ip, const zSeq *seq);

__END_DECLS

#endif /* __ZM_IP_LIN_H__ */
