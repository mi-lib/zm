/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_complex.h
 * \brief complex number class.
 * \author Zhidao
 */

#ifndef __ZM_COMPLEX_H__
#define __ZM_COMPLEX_H__

#include <zm/zm_misc.h>

__BEGIN_DECLS

/* ********************************************************** */
/*! \defgroup complex complex number.
 * \{ *//* ************************************************** */

/* ********************************************************** */
/*! \brief complex number class.
 *//* ******************************************************* */
typedef struct{
  double re; /*!< \brief real part */
  double im; /*!< \brief imaginary part */
} zComplex;

/*! \brief create a complex number.
 *
 * zComplexCreate() creates a complex number whose real part
 * is \a r and imaginary part is \a i. The result is put
 * where \a c points.
 */
#define _zComplexCreate(c,r,i) do{\
  (c)->re = (r);\
  (c)->im = (i);\
} while(0)
__ZM_EXPORT zComplex *zComplexCreate(zComplex *c, double r, double i);

/*! \brief create a complex number from the polar expression.
 *
 * zComplexCreatePolar() creates a complex number from the polar
 * expression with radius \a r and argument angle \a t
 * in Gaussian plane, where \a t is given in radians.
 * The result is put where \a c points.
 */
#define _zComplexCreatePolar(c,r,t) _zComplexCreate( c, (r)*cos(t), (r)*sin(t) )
__ZM_EXPORT zComplex *zComplexCreatePolar(zComplex *c, double r, double t);

/*! \brief copy a complex number to another.
 *
 * zComplexCopy() copies a complex number \a src to another \a dest. */
#define _zComplexCopy(src,dest) ( *(dest) = *(src) )
__ZM_EXPORT zComplex *zComplexCopy(const zComplex *src, zComplex *dest);

/*! \brief create a zero complex number.
 *
 * zComplexZero() creates a zero complex number \a c by
 * setting both real and imaginary parts for zeros. */
#define _zComplexZero(c) _zComplexCreate(c,0,0)
#define zComplexZero(c) zComplexCreate(c,0,0)

/*! \brief touchup a complex number.
 *
 * zComplexTouchup() replaces real part or imaginary part of a complex number \a c for zero if either
 * value relative to the other part is smaller than \a tol.
 * \return
 * zComplexTouchup() returns a pointer \a c.
 */
__ZM_EXPORT zComplex *zComplexTouchup(zComplex *c, double tol);

/*! \brief test if a complex number is under the tolerance.
 *
 * zComplexIsTol() tests if the absolute values of
 * both real and imaginary parts of a given complex
 * number \a c is less than the tolerance \a tol,
 * where \a tol must be positive.
 * The result is returned as a boolean value.
 */
#define zComplexIsTol(c,tol) ( zIsTol((c)->re,tol) && zIsTol((c)->im,tol) )

/*! \brief test if a complex nubmer is under the default tolerance. */
#define zComplexIsTiny(c)    zComplexIsTol( c, zTOL )

/*! \brief check if a complex number is a real number. */
#define zComplexIsReal(c,tol) zIsTol( (c)->im, tol )

/*! \brief check if two complex numbers are equal. */
#define zComplexEqual(c1,c2,tol) ( zEqual( (c1)->re, (c2)->re, tol ) && zEqual( (c1)->im, (c2)->im, tol ) )

/*! \brief check if two complex numbers are exactly the same. */
#define zComplexMatch(c1,c2) ( (c1)->re == (c2)->re && (c1)->im == (c2)->im )

/*! \brief check if two complex numbers are co-conjugate. */
#define zComplexCoconj(c1,c2,tol) ( zEqual( (c1)->re, (c2)->re, tol ) && zEqual( (c1)->im, -(c2)->im, tol ) )

/*! \brief read a complex number from a string. */
__ZM_EXPORT zComplex *zComplexFromStr(zComplex *c, const char *str);

/*! \brief read a complex number from a ZTK format processor. */
__ZM_EXPORT zComplex *zComplexFromZTK(zComplex *c, ZTK *ztk);

/*! \brief primt a complex number.
 *
 * zComplexFPrint() prints a complex number \a c to the
 * current position of a file \a fp in the following style:
 *  x + y i
 */
__ZM_EXPORT void zComplexFPrint(FILE *fp, const zComplex *c);
/*! \brief print a complex number to the standard output. */
#define zComplexPrint(c) zComplexFPrint( stdout, (c) )

/*! \brief print the coordinates of a complex number.
 *
 * zComplexCoordFPrint() prints the coordinates of a
 * complex number \a c on Gaussian plane to the current
 * position of \a file \a fp in the following style.
 *  x y
 */
__ZM_EXPORT void zComplexCoordFPrint(FILE *fp, const zComplex *c);
/*! \brief prints the coordinates of a complex number
 * to the standard output. */
#define zComplexCoordPrint(c) zComplexCoordFPrint( stdout, (c) )

/*! \brief check if a complex number is a member of an array.
 *
 * zComplexValIsIncluded() checks if a complex number \a c is
 * included in an array \a array.
 * zComplexValConjIsIncluded() checks if conjugate of \a c is
 * included in an array \a array.
 * For both functions, \a size is the size of the array, and
 * \a tol is the tolerance to regard two values are the same.
 * \return
 * zComplexValIsIncluded() returns the true value if \a c is
 * included in \a \array. Otherwise, the false value is returned.
 * zComplexValConjIsIncluded() returns the true value if conjugate
 * of \a c is included in \a \array. Otherwise, the false value
 * is returned.
 */
__ZM_EXPORT bool zComplexValIsIncluded(const zComplex *array, int size, const zComplex *c, double tol);
__ZM_EXPORT bool zComplexValConjIsIncluded(const zComplex *array, int size, const zComplex *c, double tol);

/*! \} */

__END_DECLS

#include <zm/zm_complex_arith.h> /* arithmatics */
#include <zm/zm_complex_pe.h> /* polynomial equation solver */

#endif /* __ZM_COMPLEX_H__ */
