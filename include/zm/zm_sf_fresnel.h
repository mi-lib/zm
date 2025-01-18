/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_sf_fresnel.h
 * \brief Fresnel integral.
 * \author Zhidao
 */

#ifndef __ZM_SF_FRESNEL_H__
#define __ZM_SF_FRESNEL_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief Fresnel integral
 *
 * zFresnelIntgPI_2() computes the Fresnel integral, which is defined as
 *  S(\a x) = int_0^\a x sin(pi/2 u^2) du
 *  C(\a x) = int_0^\a x cos(pi/2 u^2) du
 * zFresnelIntgScale() finds a scaled version of the Fresnel integral as
 *  S(\a x) = int_0^\a x sin(\a f u^2) du
 *  C(\a x) = int_0^\a x cos(\a f u^2) du
 * where \a f is a scale factor.
 * zFresnelIntg() computes the normalized Fresnel integral as
 *  S(\a x) = int_0^x sin(u^2) du
 *  C(\a x) = int_0^x cos(u^2) du
 * In whichever case, S(x) and C(x) are stored in \a s and \a c, respectively.
 * \return
 * zFresnelIntgPI_2(), zFresnelIntgScale() and zFresnelIntg() return the
 * false value if an internal iteration fails. Otherwise, the true value
 * is returned.
 */
__ZM_EXPORT bool zFresnelIntgPI_2(double x, double *s, double *c);
__ZM_EXPORT bool zFresnelIntgScale(double x, double f, double *s, double *c);
__ZM_EXPORT bool zFresnelIntg(double x, double *s, double *c);

/*! \brief generalized Fresnel integral
 *
 * zFresnelIntgGen() computes a generalized Fresnel integral as
 *  S(\a x) = int_0^\a x sin(\a f0 + \a f1 u + \a f2 u^2) du
 *  C(\a x) = int_0^\a x cos(\a f0 + \a f1 u + \a f2 u^2) du
 * where \a f0, \a f1 and \a f2 can be any real numbers.
 * \return
 * zFresnelIntgGen() internally calls zFresnelIntgScale() The returned
 * value conforms to it.
 */
__ZM_EXPORT bool zFresnelIntgGen(double x, double f0, double f1, double f2, double *s, double *c);

__END_DECLS

#endif /* __ZM_SF_FRESNEL_H__ */
