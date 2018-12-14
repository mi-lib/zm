/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_sf_gamma.h
 * \brief gamma and beta function.
 * \author Zhidao
 */

#ifndef __ZM_SF_GAMMA_H__
#define __ZM_SF_GAMMA_H__

#include <zm/zm_misc.h>

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* ********************************************************** */
/*! \defgroup sf special functions.
 * \{ *//* ************************************************** */

/*! \brief gamma function.
 *
 * zGamma() is the gamma function defined by the integral
 * of t^(\a x-1)*exp(-t) about t from zero to infinity for
 * an arbitrary value \a x.
 * \sa
 * zLnGamma
 */
__EXPORT double zGamma(double x);

/*! \brief natural logarithmic gamma function.
 *
 * zLnGamma() is a natural logarithm of zGamma().
 * Due to the property of logarithm functions, it is undefined
 * for negative gamma function values. It is useful to avoid
 * an overflow even for a modest \a x value.
 * \note
 * As is noted in the above description, exp(zLnGamma())
 * is not necessarily equivalent to zGamma(). Programmers
 * should pay attention to the definition range of zLnGamma()
 * when using.
 * \sa
 * zGamma
 */
__EXPORT double zLnGamma(double x);

/*! \brief beta function.
 *
 * zBeta() is the beta function defined by
 * zGamma(\a z)*zGamma(\a w)/zGamma(\a z+\a w) for a set of
 * any positive values \a z and \a w.
 */
__EXPORT double zBeta(double z, double w);

/*! \} */

__END_DECLS

#endif /* __ZM_SF_GAMMA_H__ */
