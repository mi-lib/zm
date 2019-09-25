/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_complex - complex number class.
 */

#include <zm/zm_complex.h>

/* ********************************************************** */
/* CLASS: zComplex
 * complex number class
 * ********************************************************** */

/* create a complex number. */
zComplex *zComplexCreate(zComplex *c, double r, double i)
{
  c->re = r;
  c->im = i;
  return c;
}

/* create a complex number based on the polar expression. */
zComplex *zComplexCreatePolar(zComplex *c, double r, double t)
{
  return zComplexCreate( c, r*cos(t), r*sin(t) );
}

/* touchup a complex number. */
zComplex *zComplexTouchup(zComplex *c)
{
  double ri, ir;

  ri = c->re / c->im;
  ir = c->im / c->re;
  if( zIsTiny(ri) ) c->re = 0;
  if( zIsTiny(ir) ) c->im = 0;
  return c;
}

/* print a complex number to a file. */
void zComplexFPrint(FILE *fp, zComplex *c)
{
  fprintf( fp, "%.10g", c->re );
  if( c->im > 0 )
    fprintf( fp, " + %.10g i", c->im );
  else if( c->im < 0 )
    fprintf( fp, " - %.10g i",-c->im );
}

/* print the coordinates of a complex number on the
 * Gaussian plane to a file. */
void zComplexCoordFPrint(FILE *fp, zComplex *c)
{
  fprintf( fp, "%.10g %.10g", c->re, c->im );
}

/* check if a complex number is a member of an array. */
bool zComplexValIsIncluded(zComplex *array, int size, zComplex *c)
{
  register int i;

  for( i=0; i<size; i++ )
    if( zComplexIsEqual( &array[i], c ) ) return true;
  return false;
}

/* check if conjugate of a complex number is a member of an array. */
bool zComplexValConjIsIncluded(zComplex *array, int size, zComplex *c)
{
  zComplex cc;

  zComplexConj( c, &cc );
  return zComplexValIsIncluded( array, size, &cc );
}
