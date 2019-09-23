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
__EXPORT zComplex *zComplexCreate(zComplex *c, double r, double i);

/*! \brief create a complex number from the polar expression.
 *
 * zComplexCreatePolar() creates a complex number from the polar
 * expression with radius \a r and argument angle \a t
 * in Gaussian plane, where \a t is given in radians.
 * The result is put where \a c points.
 */
__EXPORT zComplex *zComplexCreatePolar(zComplex *c, double r, double t);

/*! \brief copy a complex number to another.
 *
 * zComplexCopy() copies a complex number \a src to another \a dest. */
#define zComplexCopy(src,dest) ( *(dest) = *(src) );

/*! \brief create a zero complex number.
 *
 * zComplexZero() creates a zero complex number \a c by
 * setting both real and imaginary parts for zeros. */
#define zComplexZero(c) zComplexCreate(c,0,0)

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

/*! \brief primt a complex number.
 *
 * zComplexFPrint() prints a complex number \a c to the
 * current position of a file \a fp in the following style:
 *  x + y i
 */
__EXPORT void zComplexFPrint(FILE *fp, zComplex *c);
/*! \brief print a complex number to the standard output. */
#define zComplexPrint(c) zComplexFPrint( stdout, (c) )

/*! \brief print the coordinates of a complex number.
 *
 * zComplexCoordFPrint() prints the coordinates of a
 * complex number \a c on Gaussian plane to the current
 * position of \a file \a fp in the following style.
 *  x y
 */
__EXPORT void zComplexCoordFPrint(FILE *fp, zComplex *c);
/*! \brief prints the coordinates of a complex number
 * to the standard output. */
#define zComplexCoordPrint(c) zComplexCoordFPrint( stdout, (c) )

/*! \} */

__END_DECLS

#include <zm/zm_complex_arith.h> /* arithmatics */
#include <zm/zm_complex_pe.h> /* polynomial equation solver */

#endif /* __ZM_COMPLEX_H__ */
