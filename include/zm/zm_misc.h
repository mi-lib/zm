/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_misc.h
 * \brief miscellanies.
 * \author Zhidao
 */

#ifndef __ZM_MISC_H__
#define __ZM_MISC_H__

#include <zeda/zeda.h>
#include <math.h>

#include <zm/zm_pi.h>
#include <zm/zm_ieee.h>
#include <zm/zm_errmsg.h>

__BEGIN_DECLS

/*! \brief tolerance on computation; tiny quantity. */
/*! \note machine epsilon conforming to IEEE754 is about 2.22e-16 */
#define zTOL ( 1.0e-12 )

/*! \brief check if value is under the tolerance level.
 *
 * zIsTol() checks if the absolute value of a given
 * \a x is less than or equal to \a tol.
 * \note \a tol must be positive.
 */
#define zIsTol(x,tol) ( fabs(x) <= (tol) )

/*! \brief check if value is under the defined tolerance zTOL. */
#define zIsTiny(x)    zIsTol( x, zTOL )

/*! \brief checks if an integer value \a d is an even number. */
#define zIsEven(d)    (~(d) & 0x1 )
/*! \brief checks if an integer value \a d is an odd number. */
#define zIsOdd(d)     ( (d) & 0x1 )

/*! \brief check if signs of two values are opposite.
 *
 * zIsSgnOpp() checks if the signs of two values \a a and \a b
 * are opposite or not.
 * \retval true signs of \a a and \a b are opposite.
 * \retval false signs of \a a and \a b are the same.
 */
__EXPORT bool zIsSgnOpp(double a, double b);

/*! \brief a sine and cosine set.
 *
 * zSinCos() computes a set of sine and cosine values
 * of a given radian value \a angle. sin(\a angle) is
 * set for \a s, while cos(\a angle) for \a c.
 */
#define _zSinCos(a,s,c) do{\
  *(s) = sin(a);\
  *(c) = cos(a);\
} while(0)
__EXPORT void zSinCos(double angle, double *s, double *c);

/*! \brief normalize phase within the range from minus pi to pi.
 *
 * zPhaseNormalize() normalizes a radian value \a angle
 * within the range from minus pi to pi. It returns a value
 * \a x that satisfies \a angle = \a x + 2*pi*n, where n is
 * an integer number which satisfies -pi<\a x<=pi.
 */
__EXPORT double zPhaseNormalize(double angle);

/*! \brief the base of natural logarithms. */
#define zE 2.7182818284590452354

/*! \brief maximum number of iteration times. */
#define Z_MAX_ITER_NUM 10000
/*! \brief set the muximum iteration number. */
#define ZITERINIT(iter) do{\
  if( (iter) <= 0 ) (iter) = Z_MAX_ITER_NUM; } while(0)
/*! \brief warn when an iteration counts over the maximum. */
#define ZITERWARN(iter) \
  ZRUNWARN( "iteration unconverged within %d steps", (iter) )

/*! \brief sign function. */
#define zSgn(x)     ( (x)>=0 ? 1 : -1 )
/*! \brief the absolute-ceiling value.
 *  ( ex. zCeil(0.5) = 1.0 , zCeil(-0.5) = -1.0 ) */
#define zCeil(x)    ( (x)>=0 ? ceil(x) : floor(x) )
/*! \brief rounded value.
 * \note This is defined for a compiler which is not conforming C99. */
#define zRound(x)   (double)( (x)>=0 ? (int)((x)+0.5) : (int)((x)-0.5) )
/*! \brief a fructuation of a value \a x devided by \a y.
 * \note if \a y is zero, anything might happen. */
#define zFruct(x,y) ( (x) - floor((x)/(y))*(y) )
/*! \brief the squared value of \a x. */
__EXPORT double zSqr(double x);
/*! \brief the cubed value of \a x. */
__EXPORT double zCube(double x);

/*! \brief a non-negative integer \a n power of \a x. */
__EXPORT double zPowN(double x, unsigned int n);

/*! \brief a line inter-extrapolation, connecting two
 * points on a plane, (\a x0, \a y0) and (\a x1, \a y1). */
__EXPORT double zLine(double x, double x0, double y0, double x1, double y1);

/*! \brief the x value of a cycloid curve which goes from
 * (\a from, \a base) to (\a to, \a base) with a height \a h
 * with respect to \a phase, mostly expected in between 0 and 1.
 * \sa zCycloidY
 */
__EXPORT double zCycloidX(double from, double to, double phase);
/*! \brief the y value of a cycloid curve which goes from
 * (\a from, \a base) to (\a to, \a base) with a height \a h
 * with respect to \a phase, mostly expected in between 0 and 1.
 * \sa zCycloidX
 */
__EXPORT double zCycloidY(double base, double h, double phase);

/*! \brief a sigmoid function. \a s is the time constant. */
__EXPORT double zSigmoid(double x, double s);

/*! \brief the cubic root of \a x. */
__EXPORT double zCbrt(double x);

__END_DECLS

#endif /* __ZM_MISC_H__ */
