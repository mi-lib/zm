/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_sf_erf.h
 * \brief Gauss's error function.
 * \author Zhidao
 */

#ifndef __ZM_SF_ERF_H__
#define __ZM_SF_ERF_H__

#include <zm/zm_misc.h>

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* ********************************************************** */
/*! \defgroup sf special functions.
 * \{ *//* ************************************************** */

/*! \brief error function approximated by Taylor series.
 *
 * zErfN() is an approximation of Gauss's error function, which
 * is defined by the integral of 2/sqrt(pi) exp(-t^2) about t
 * from minus infinity to \a x. The approximation is done by
 * Taylor series expantion up to \a n 'th term.
 */
__EXPORT double zErfN(double x, int n);

/*! \brief error function.
 *
 * zErf() is Gauss's error function, which is defined by the
 * integral of 2/sqrt(pi) exp(-t^2) about t from minus infinity
 * to \a x. It is equivalent with zErfN( x, Z_ERF_DEFAULT_ORDER ).
 */
#define Z_ERF_DEFAULT_ORDER 35
#define zErf(x) zErfN( x, Z_ERF_DEFAULT_ORDER )

/*! \} */

__END_DECLS

#endif /* __ZM_SF_ERF_H__ */
