/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_sf_fresnel.h
 * \brief Fresnel integral.
 * \author Zhidao
 */

#ifndef __ZM_SF_FRESNEL_H__
#define __ZM_SF_FRESNEL_H__

#include <zm/zm_misc.h>

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief Fresnel integral
 *
 * zFresnelIntgPI_2() computes the Fresnel integral, which is defined as
 *  S(x) = int_0^x sin(pi/2 u^2) du
 *  C(x) = int_0^x cos(pi/2 u^2) du
 * zFresnelIntg() finds a scaled version of the Fresnel integral as
 *  S(x) = int_0^x sin(u^2) du
 *  C(x) = int_0^x cos(u^2) du
 * In either case, S(x) and C(x) are stored in \a s and \a c, respectively.
 * \return
 * zFresnelIntgPI_2() and zFresnelIntg() return the false value if an
 * internal iteration fails. Otherwise, the true value is returned.
 */
__EXPORT bool zFresnelIntgPI_2(double x, double *s, double *c);
__EXPORT bool zFresnelIntg(double x, double *s, double *c);

__END_DECLS

#endif /* __ZM_SF_FRESNEL_H__ */
