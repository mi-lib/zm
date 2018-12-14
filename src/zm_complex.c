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

/* zComplexCreate
 * - creation of a complex number.
 */
zComplex *zComplexCreate(zComplex *c, double r, double i)
{
  c->re = r;
  c->im = i;
  return c;
}

/* zComplexPolar
 * - creation of a complex number from a polar expression.
 */
zComplex *zComplexPolar(zComplex *c, double r, double t)
{
  return zComplexCreate( c, r*cos(t), r*sin(t) );
}

/* zComplexFWrite
 * - output of complex number to file.
 */
void zComplexFWrite(FILE *fp, zComplex *c)
{
  fprintf( fp, "%.10g", c->re );
  if( c->im > 0 )
    fprintf( fp, " + %.10g i", c->im );
  else if( c->im < 0 )
    fprintf( fp, " - %.10g i",-c->im );
}

/* zComplexCoordFWrite
 * - output of the coordinates of complex number on the
 * Gaussian plane to file.
 */
void zComplexCoordFWrite(FILE *fp, zComplex *c)
{
  fprintf( fp, "%.10g %.10g", c->re, c->im );
}
