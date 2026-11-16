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

#include <zm/zm_export.h>

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

#define _zTestExpr(testval,baseval,expr) ( zIsNan(testval) ? false : ( zIsNan(baseval) ? true : expr(testval,baseval) ) )
#define _zIsMax(testval,baseval)    ( (testval) > (baseval) )
#define _zIsMin(testval,baseval)    ( (testval) < (baseval) )
#define _zIsAbsMax(testval,baseval) ( fabs(testval) > fabs(baseval) )
#define _zIsAbsMin(testval,baseval) ( fabs(testval) < fabs(baseval) )

/*! \brief check if a value is larger than a base value. */
#define zIsMax(testval,baseval)    _zTestExpr( testval, baseval, _zIsMax )
/*! \brief check if a value is smaller than a base value. */
#define zIsMin(testval,baseval)    _zTestExpr( testval, baseval, _zIsMin )
/*! \brief check if the absolute of a value is larger than a base value. */
#define zIsAbsMax(testval,baseval) _zTestExpr( testval, baseval, _zIsAbsMax )
/*! \brief check if the absolute of a value is smaller than a base value. */
#define zIsAbsMin(testval,baseval) _zTestExpr( testval, baseval, _zIsAbsMin )

/*! \brief check if a test value is less than a base value. */
__ZM_EXPORT bool zIsLess(double testval, double baseval);

/*! \brief check if a test value is greater than a base value. */
__ZM_EXPORT bool zIsGreater(double testval, double baseval);

/*! \brief check if two values are equal. */
__ZM_EXPORT bool zEqual(double a, double b, double tol);

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
__ZM_EXPORT void zSinCos(double angle, double *s, double *c);

/*! \brief normalize phase within the range from minus pi to pi.
 *
 * zPhaseNormalize() normalizes a radian value \a angle
 * within the range from minus pi to pi. It returns a value
 * \a x that satisfies \a angle = \a x + 2*pi*n, where n is
 * an integer number which satisfies -pi<\a x<=pi.
 */
__ZM_EXPORT double zPhaseNormalize(double angle);

/*! \brief logarithmic function with an arbitrary base. */
#define zLog(base,x) ( log(x) / log(base) )

/*! \brief Napier's constant (the base of natural logarithms). */
#define zE 2.7182818284590452354

/*! \brief maximum number of iteration times. */
#define Z_MAX_ITER_NUM 10000
/*! \brief set the muximum iteration number. */
#define ZITERINIT(iter) do{ if( (iter) <= 0 ) (iter) = Z_MAX_ITER_NUM; } while(0)
/*! \brief warn when an iteration counts over the maximum. */
#define ZITERWARN(iter) ZRUNWARN( ZM_WARN_ITERATION, (iter) )

/*! \brief sign function. */
#define zSgn(x)     ( (x)>=0 ? 1 : -1 )
/*! \brief the absolute-ceiling value.
 *  ( ex. zCeil(0.5) = 1.0 , zCeil(-0.5) = -1.0 ) */
#define zCeil(x)    ( (x)>=0 ? ceil(x) : floor(x) )
/*! \brief rounded value.
 * \note This is defined for a compiler which is not conforming to C99. */
#define zRound(x)   (double)( (x)>=0 ? (int)((x)+0.5) : (int)((x)-0.5) )
/*! \brief a fructuation of a value \a x devided by \a y.
 * \note if \a y is zero, anything might happen. */
#define zFruct(x,y) ( (x) - floor((x)/(y))*(y) )
/*! \brief the squared value of \a x. */
#define _zSqr(x)    ( (x) * (x) )
__ZM_EXPORT double zSqr(double x);
/*! \brief the cubed value of \a x. */
#define _zCube(x)   ( (x) * (x) * (x) )
__ZM_EXPORT double zCube(double x);

/*! \brief a non-negative integer \a n power of \a x. */
__ZM_EXPORT double zPowN(double x, unsigned int n);

/*! \brief a line inter-extrapolation, connecting two
 * points on a plane, (\a x0, \a y0) and (\a x1, \a y1). */
__ZM_EXPORT double zLine(double x, double x0, double y0, double x1, double y1);

/*! \brief the x value of a cycloid curve with respect to \a phase,
 * which is mostly expected to be between 0 and 1.
 * \sa zCycloidY
 */
__ZM_EXPORT double zCycloidX(double phase);

/*! \brief the y value of a cycloid curve with respect to \a phase,
 * which is mostly expected to be between 0 and 1.
 * \sa zCycloidX
 */
__ZM_EXPORT double zCycloidY(double phase);

/*! \brief sigmoid function. */
__ZM_EXPORT double zSigmoid(double x);

/*! \brief the cubic root of \a x. */
__ZM_EXPORT double zCbrt(double x);

/*! \brief the power of two not greater than the absolute value of a value.
 *
 * zPow2Floor() returns the power of two 2^(e-1) such that the absolute value of \a x divided by it is
 * in [1,2). Since multiplications and divisions by it are exact unless they overflow or underflow,
 * it is useful to scale values to avoid overflow and underflow.
 * \return
 * zPow2Floor() returns the power of two. It returns zero for zero, infinity for infinity, and NaN for
 * NaN, respectively.
 */
__ZM_EXPORT double zPow2Floor(double x);

/*! \brief error-free transformations of the sum and the product of two values.
 *
 * zAdd() computes the sum of \a a and \a b in floating-point arithmetic, and puts the rounding error of
 * it into \a error, so that \a a + \a b = (the returned value) + \a error exactly (Knuth's TwoSum).
 *
 * zMul() computes the product of \a a and \a b in floating-point arithmetic, and puts the rounding error
 * of it into \a error, so that \a a \a b = (the returned value) + \a error exactly.
 *
 * zSplitMant() splits the mantissa of \a a into the upper and the lower halves, namely, two values with
 * at most DBL_MANT_DIG/2 (rounded down; 26 for IEEE 754 double-precision floating-point numbers)
 * significant bits, of which the lower is put into \a lower, so that \a a = (the returned value) +
 * \a lower exactly (Veltkamp's splitting). The product of two such halves is exactly representable.
 * It is used in zMul() if fma() is not available.
 *
 * They are useful to compute sums and inner products accurately, e.g., by Dot2 algorithm by T. Ogita,
 * S. M. Rump and S. Oishi:
 *  T. Ogita, S. M. Rump and S. Oishi, Accurate Sum and Dot Product, SIAM Journal of Science and Computing,
 *  26(6):1955-1988, 2005.
 *
 * They assume binary floating-point arithmetic with rounding to nearest. They are not affected by
 * contractions of multiplications and additions into fused multiply-add operations, but are broken by
 * optimizations which reorder floating-point operations, e.g., -ffast-math of GCC.
 * The transformations are exact unless overflow or underflow happens. Note that zSplitMant() overflows
 * if the absolute value of the argument is greater than about 1e+300, and so does zMul() if fma() is not
 * available.
 * \notes
 * In zMul(), the error is computed by fma() if it is available (C99 or C++11 and later), or by Dekker's
 * TwoProduct otherwise.
 * \return
 * zAdd() returns the sum computed in floating-point arithmetic.
 * zMul() returns the product computed in floating-point arithmetic.
 * zSplitMant() returns the upper half.
 */
__ZM_EXPORT double zAdd(double a, double b, double *error);
__ZM_EXPORT double zMul(double a, double b, double *error);
__ZM_EXPORT double zSplitMant(double a, double *lower);

/*! \brief permutation.
 *
 * zPermutation() computes \a n_P_\a i.
 * The result is returned as a double-precision floating-point value.
 */
__ZM_EXPORT double zPermutation(int n, int i);

/*! \brief factorial.
 *
 * zFactorial() computes \a n!.
 * The result is returned as a double-precision floating-point value.
 */
__ZM_EXPORT double zFactorial(int n);

/*! \brief combination.
 *
 * zCombination() computes \a n_C_i.
 * zCombinationRecursive() also computes \a n_C_i but in a recursive way.
 * The result is returned as a double-precision floating-point value.
 */
__ZM_EXPORT double zCombination(int n, int i);
__ZM_EXPORT double zCombinationRecursive(int n, int i);

/*! \brief series of combination.
 *
 * zCombinationSeries() computes a series of combination \a n_C_i where \a i is from 0 to \a n.
 * The result is stored in the array of double-precision floating-point values \a c. \a size is
 * the size of \a c.
 * If \a n is more than \a size, equal to \a size or less than 0, it fails to compute the series
 * and the null pointer is returned. Otherwise, a pointer to \a c is returned.
 */
__ZM_EXPORT double *zCombinationSeries(uint n, size_t size, double c[]);

/*! \brief smoothstep function. */
__ZM_EXPORT double zSmoothStep(double x, uint order);
/*! \brief derivative of smoothstep function. */
__ZM_EXPORT double zSmoothStepDif(double x, uint order);

__END_DECLS

#endif /* __ZM_MISC_H__ */
