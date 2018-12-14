/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_lagrange - interpolation: Lagrange's interpolation.
 */

#ifndef __ZM_IP_LAGRANGE_H__
#define __ZM_IP_LAGRANGE_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* Lagrange interpolation, which connects the 'n' sections
 * by 'n' th order polynomial curve. */
__EXPORT bool zIPCreateLagrange(zIP *ip, zSeq *seq);

__END_DECLS

#endif /* __ZM_IP_LAGRANGE_H__ */
