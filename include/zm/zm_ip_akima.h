/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_akima - interpolation: Akima's interpolation (1970).
 */

#ifndef __ZM_IP_AKIMA_H__
#define __ZM_IP_AKIMA_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* Akima's interpolation, which connects the sections in a similar
   way to a handwriting. It was proposed by H. Akima in 1970. */
__EXPORT bool zIPCreateAkima(zIP *ip, zSeq *seq);

__END_DECLS

#endif /* __ZM_IP_AKIMA_H__ */
