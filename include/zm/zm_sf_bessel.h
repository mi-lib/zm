/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_sf_bessel.h
 * \brief Bessel functions of integer order.
 *
 * This implementation of the Bessel function family is
 * a rearrangement of the code written by Mr. Takuya Ooura,
 * conforming to his redistribution policy.
 * The original sources are available at:
 *  http://www.kurims.kyoto-u.ac.jp/~ooura/bessel-j.html
 * \author Zhidao
 */

#ifndef __ZM_SF_BESSEL_H__
#define __ZM_SF_BESSEL_H__

#include <zm/zm_misc.h>

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* ********************************************************** */
/*! \ingroup sf
 * \{ *//* ************************************************** */

/*! \brief the first kind of the zero-order Bessel function. */
__EXPORT double zBesselJ0(double x);

/*! \brief the first kind of the first-order Bessel function. */
__EXPORT double zBesselJ1(double x);

/*! \brief the first kind of n-order Bessel function.
 * \a n is an integer index. */
__EXPORT double zBesselJ(int n, double x);

/*! \brief the second kind of the zero-order Bessel function.
 * \note
 * From the definition, it is undefined for non-positive \a x.
 */
__EXPORT double zBesselY0(double x);

/*! \brief the second kind of the first-order Bessel function.
 * \note
 * From the definition, it is undefined for non-positive \a x.
 */
__EXPORT double zBesselY1(double x);

/*! \brief the second kind of n-order Bessel function.
 * \a n is an integer index.
 * \note
 * From the definition, it is undefined for non-positive \a x.
 */
__EXPORT double zBesselY(int n, double x);

/*! \brief the first kind of the zero-order modified Bessel function. */
__EXPORT double zBesselI0(double x);

/*! \brief the first kind of the first-order modified Bessel function. */
__EXPORT double zBesselI1(double x);

/*! \brief the first kind of n-order modified Bessel function.
 * \a n is an integer index. */
__EXPORT double zBesselI(int n, double x);

/*! \brief the second kind of the zero-order modified Bessel function.
 * \note
 * From the definition, it is undefined for non-positive \a x.
 */
__EXPORT double zBesselK0(double x);

/*! \brief the second kind of the first-order modified Bessel function.
 * \note
 * From the definition, it is undefined for non-positive \a x.
 */
__EXPORT double zBesselK1(double x);

/*! \brief the second kind of n-order modified Bessel function.
 * \a n is an integer index.
 * \note
 * From the definition, it is undefined for non-positive \a x.
 */
__EXPORT double zBesselK(int n, double x);

/*! \} */

__END_DECLS

#endif /* __ZM_SF_BESSEL_H__ */
