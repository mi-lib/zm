/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_complex - complex number class.
 */

#include <zm/zm_complex.h>

/* create a complex number. */
zComplex *zComplexCreate(zComplex *c, double r, double i)
{
  _zComplexCreate( c, r, i );
  return c;
}

/* create a complex number based on the polar expression. */
zComplex *zComplexCreatePolar(zComplex *c, double r, double t)
{
  _zComplexCreatePolar( c, r, t );
  return c;
}

/* copy a complex number to another. */
zComplex *zComplexCopy(const zComplex *src, zComplex *dest)
{
  _zComplexCopy( src, dest );
  return dest;
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

/* read a complex number from a string. */
zComplex *zComplexFromStr(zComplex *c, const char *str)
{
  double re, im;
  char buf[BUFSIZ], tkn[BUFSIZ];

  strcpy( buf, str );
  if( *buf == 'i' ) return zComplexCreate( c, 0, 1 );
  if( *buf == '+' && *(buf+1) == 'i' ) return zComplexCreate( c, 0, 1 );
  if( *buf == '-' && *(buf+1) == 'i' ) return zComplexCreate( c, 0, -1 );
  re = atof( zSNumToken( buf, tkn, BUFSIZ ) );
  if( *buf == 'i' )
    return zComplexCreate( c, 0, re );

  if( *buf == 'i' ) return zComplexCreate( c, re, 1 );
  if( *buf == '+' && *(buf+1) == 'i' ) return zComplexCreate( c, re, 1 );
  if( *buf == '-' && *(buf+1) == 'i' ) return zComplexCreate( c, re, -1 );
  zSNumToken( buf, tkn, BUFSIZ );
  im = *buf == 'i' ? atof( tkn ) : 0;
  return zComplexCreate( c, re, im );
}

/* read a complex number from a ZTK format processor. */
zComplex *zComplexFromZTK(zComplex *c, ZTK *ztk)
{
  return zComplexFromStr( c, ZTKVal(ztk) );
}

/* print imaginary part of a complex number to a file. */
static void _zComplexImFPrint(FILE *fp, const zComplex *c, char ps)
{
  double im;

  fprintf( fp, "%c", c->im > 0 ? ps : '-' );
  if( ( im = fabs( c->im ) ) != 1 )
    fprintf( fp, "%.10g", im );
  fprintf( fp, "i" );
}

/* print a complex number to a file. */
void zComplexFPrint(FILE *fp, const zComplex *c)
{
  if( c->re == 0 ){
    if( c->im == 0 )
      fprintf( fp, "0" );
    else
      _zComplexImFPrint( fp, c, '\0' );
  } else{
    fprintf( fp, "%.10g", c->re );
    if( c->im != 0 )
      _zComplexImFPrint( fp, c, '+' );
  }
}

/* print the coordinates of a complex number on the
 * Gaussian plane to a file. */
void zComplexCoordFPrint(FILE *fp, const zComplex *c)
{
  fprintf( fp, "%.10g %.10g", c->re, c->im );
}

/* check if a complex number is a member of an array. */
bool zComplexValIsIncluded(const zComplex *array, int size, const zComplex *c, double tol)
{
  int i;

  for( i=0; i<size; i++ )
    if( zComplexEqual( &array[i], c, tol ) ) return true;
  return false;
}

/* check if conjugate of a complex number is a member of an array. */
bool zComplexValConjIsIncluded(const zComplex *array, int size, const zComplex *c, double tol)
{
  int i;

  for( i=0; i<size; i++ )
    if( zComplexCoconj( &array[i], c, tol ) ) return true;
  return false;
}
