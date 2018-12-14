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

/*! \brief zero of the complex number */
/*! \cond */
extern const zComplex zcomplexzero;
/*! \endcond */
#define ZCOMPLEXZERO ( (zComplex *)&zcomplexzero )

/*! \brief create a complex number.
 *
 * zComplexCreate() creates a complex number whose real part
 * is \a r and imaginary part is \a i. The result is put
 * where \a c points.
 */
__EXPORT zComplex *zComplexCreate(zComplex *c, double r, double i);

/*! \brief create a complex number from a polar expression.
 *
 * zComplexPolar() creates a complex number from a polar
 * expression with radius \a r and argument angle \a t
 * in Gaussian plane, where \a t is given in radians.
 * The result is put where \a c points.
 */
__EXPORT zComplex *zComplexPolar(zComplex *c, double r, double t);

/*! \brief copy a complex number to another.
 *
 * zComplexCopy() copies a complex number \a src to another \a dest. */
#define zComplexCopy(src,dest) ( *(dest) = *(src) );

/*! \brief create a zero complex number.
 *
 * zComplexClear() creates a zero complex number \a c,
 * setting both real and imaginary parts for zero. */
#define zComplexClear(c) zComplexCreate(c,0,0)

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

/*! \brief write a complex number.
 *
 * zComplexFWrite() writes a complex number \a c to
 * the current position of a file \a fp in the following
 * style:
 *  x + y i
 */
__EXPORT void zComplexFWrite(FILE *fp, zComplex *c);
/*! \brief write a complex number to the standard output. */
#define zComplexWrite(c) zComplexFWrite( stdout, (c) )

/*! \brief write the coordinates of a complex number.
 *
 * zComplexCoordFWrite() writes the coordinates of a
 * complex number \a c on Gaussian plane to the current
 * position of \a file \a fp in the following style.
 *  x y
 */
__EXPORT void zComplexCoordFWrite(FILE *fp, zComplex *c);
/*! \brief writes the coordinates of a complex number
 * to the standard output. */
#define zComplexCoordWrite(c) zComplexCoordFWrite( stdout, (c) )

/*! \} */

__END_DECLS

#include <zm/zm_complex_arith.h> /* arithmatics */
#include <zm/zm_complex_pe.h> /* polynomial equation solver */

#endif /* __ZM_COMPLEX_H__ */
