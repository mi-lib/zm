/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_akima - interpolation: Akima's interpolation (1970).
 */

#ifndef __ZM_IP_AKIMA_H__
#define __ZM_IP_AKIMA_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief create Akima interpolator.
 *
 * zIPCreateAkima() creates an Akima interpolator \a ip proposed by H. Akima in 1970.
 * \a seq is a sequence of points to be interpolated.
 * It connects the points in a similar way to a handwriting.
 *
 * zIPCreateModifiedAkima() creates a modified Akima interporator \a ip from a sequence of points \a seq
 * to be interpolated.
 * \return
 * zIPCreateAkima() and zIPCreateModifiedAkima() return a pointer \a ip when they succeed to create the
 * interpolator. Otherwise, the null pointer is returned.
 */
__ZM_EXPORT bool zIPCreateAkima(zIP *ip, const zSeq *seq);
__ZM_EXPORT bool zIPCreateModifiedAkima(zIP *ip, const zSeq *seq);

__END_DECLS

#endif /* __ZM_IP_AKIMA_H__ */
