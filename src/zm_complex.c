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

/* zero of the complex number */
const zComplex zcomplexzero = { 0, 0 };

/* create a complex number. */
zComplex *zComplexCreate(zComplex *c, double r, double i)
{
  c->re = r;
  c->im = i;
  return c;
}

/* create a complex number based on the polar expression. */
zComplex *zComplexPolar(zComplex *c, double r, double t)
{
  return zComplexCreate( c, r*cos(t), r*sin(t) );
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
