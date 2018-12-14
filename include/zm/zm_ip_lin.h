/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_lin - interpolation: linear interpolation.
 */

#ifndef __ZM_IP_LIN_H__
#define __ZM_IP_LIN_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* linear interpolation, which connects the sections by lines. */
__EXPORT bool zIPCreateLinear(zIP *ip, zSeq *seq);

__END_DECLS

#endif /* __ZM_IP_LIN_H__ */
