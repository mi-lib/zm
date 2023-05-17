/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_chebyshev - interpolation: Chebyshev's interpolation.
 */

#ifndef __ZM_IP_CHEBYSHEV_H__
#define __ZM_IP_CHEBYSHEV_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* Chebyshev's interpolation, which interpolates n points by
 * a combination of Chebyshev base functions. */
__ZM_EXPORT bool zIPCreateChebyshev(zIP *ip, zSeq *seq);

__END_DECLS

#endif /* __ZM_IP_CHEBYSHEV_H__ */
