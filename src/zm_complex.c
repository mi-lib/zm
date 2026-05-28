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
zComplex *zComplexTouchup(zComplex *c, double tol)
{
  double ri, ir;

  ri = c->re / c->im;
  ir = c->im / c->re;
  if( zIsTol(ri,tol) || zIsTol(c->re,tol) ) c->re = 0;
  if( zIsTol(ir,tol) || zIsTol(c->im,tol) ) c->im = 0;
  return c;
}

/* arithmatics */

/* squared absolute value of a complex number. */
double zComplexSqrAbs(const zComplex *c)
{
  return _zComplexSqrAbs( c );
}

/* absolute value of a complex number. */
double zComplexAbs(const zComplex *c)
{
  return _zComplexAbs( c );
}

/* argument angle of a complex number. */
double zComplexArg(const zComplex *c)
{
  return _zComplexArg( c );
}

/* complex conjugate. */
zComplex *zComplexConj(const zComplex *c, zComplex *cc)
{
  _zComplexConj( c, cc );
  return cc;
}

/* add two complex numbers. */
zComplex *zComplexAdd(const zComplex *c1, const zComplex *c2, zComplex *c)
{
  _zComplexAdd( c1, c2, c );
  return c;
}

/* subtract complex number from the other. */
zComplex *zComplexSub(const zComplex *c1, const zComplex *c2, zComplex *c)
{
  _zComplexSub( c1, c2, c );
  return c;
}

/* reverse complex number. */
zComplex *zComplexRev(const zComplex *c, zComplex *rc)
{
  _zComplexRev( c, rc );
  return rc;
}

/* multiply a complex number by a real number. */
zComplex *zComplexMul(const zComplex *c, double k, zComplex *ec)
{
  _zComplexMul( c, k, ec );
  return ec;
}

/* divide a complex number by a real number. */
zComplex *zComplexDiv(const zComplex *c, double k, zComplex *rc)
{
  if( zIsTiny( k ) ){
    ZRUNERROR( ZM_ERR_ZERODIV );
    return NULL;
  }
  k = 1 / k;
  _zComplexMul( c, k, rc );
  return rc;
}

/* multiply two complex numbers. */
zComplex *zComplexCMul(const zComplex *c1, const zComplex *c2, zComplex *c)
{
  _zComplexCMul( c1, c2, c );
  return c;
}

/* multiply a complex numbers by the conjugate of another complex number. */
zComplex *zComplexCMulConj(const zComplex *c1, const zComplex *c2, zComplex *c)
{
  _zComplexCMulConj( c1, c2, c );
  return c;
}

/* invert a complex number. */
zComplex *zComplexInv(const zComplex *c, zComplex *ic)
{
  double r2;
  zComplex cc;

  r2 = _zComplexSqrAbs( c );
  _zComplexConj( c, &cc );
  return zComplexDiv( &cc, r2, ic );
}

/* divide a complex number by another. */
zComplex *zComplexCDiv(const zComplex *c1, const zComplex *c2, zComplex *c)
{
  zComplex cc;

  zComplexInv( c2, &cc );
  _zComplexCMul( c1, &cc, c );
  return c;
}

/* multiply two complex numbers. */
zComplex *zComplexCMulDRC(zComplex *c1, const zComplex *c2)
{
  zComplex tmp;

  _zComplexCMul( c1, c2, &tmp );
  _zComplexCopy( &tmp, c1 );
  return c1;
}

/* multiply a complex numbers by the conjugate of another complex number. */
zComplex *zComplexCMulConjDRC(zComplex *c1, const zComplex *c2)
{
  zComplex tmp;

  _zComplexCMulConj( c1, c2, &tmp );
  _zComplexCopy( &tmp, c1 );
  return c1;
}

/* divide a complex number by another. */
zComplex *zComplexCDivDRC(zComplex *c1, const zComplex *c2)
{
  zComplex tmp;

  zComplexCDiv( c1, c2, &tmp );
  _zComplexCopy( &tmp, c1 );
  return c1;
}

/* power a complex number by a real number. */
zComplex *zComplexPow(const zComplex *c, double z, zComplex *pc)
{
  double r, theta;

  r = _zComplexAbs( c );
  theta = _zComplexArg( c );
  _zComplexCreatePolar( pc, pow( r, z ), theta * z );
  return pc;
}

/* power complex number by a real number. */
zComplex *zComplexPowRef(const zComplex *c, double z, const zComplex *ref, zComplex *pc)
{
  double r, theta, w;

  r = pow( _zComplexAbs(c), z );
  theta = _zComplexArg(c) * z;
  w = 2 * zPI * z;
  if( ref && !zIsTiny(z) )
    theta += w * ceil( ( _zComplexArg(ref) - theta ) / w );
  _zComplexCreatePolar( pc, r, theta );
  return pc;
}

/* power complex number by another complex number. */
zComplex *zComplexCPow(const zComplex *c, const zComplex *z, zComplex *pc)
{
  double r, theta;

  r = _zComplexAbs( c );
  theta = _zComplexArg( c );
  return zComplexCreatePolar( pc,
    pow( r, z->re ) * exp( -z->im * theta ), z->re*theta + z->im*log(r) );
}

/* the real number base logarithm of complex number. */
zComplex *zComplexLog(const zComplex *c, double base, zComplex *lc)
{
  double r, theta, d;

  if( base <= 0 || base == 1.0 ){
    ZRUNERROR( ZM_ERR_INVALID_LOGBASE );
    return NULL;
  }
  r = _zComplexAbs( c );
  theta = _zComplexArg( c );
  d = log( base );
  _zComplexCreate( lc, log(r)/d, theta/d );
  return lc;
}

/* the complex number base logarithm of complex number. */
zComplex *zComplexCLog(const zComplex *c, const zComplex *base, zComplex *lc)
{
  double rb, thetab, lrb, d, rc, thetac, lrc;

  rb = _zComplexAbs( base );
  thetab = _zComplexArg( base );
  if( rb == 0 ){
    ZRUNERROR( ZM_ERR_INVALID_LOGBASE );
    return NULL;
  }
  lrb = log( rb );
  if( ( d = lrb*lrb + thetab*thetab ) == 0 ){
    ZRUNERROR( ZM_ERR_INVALID_LOGBASE );
    return NULL;
  }
  rc = _zComplexAbs( c );
  thetac = _zComplexArg( c );
  lrc = log( rc );
  _zComplexCreate( lc,
    (lrb*lrc+thetab*thetac)/d, (thetac*lrb-thetab*lrc)/d );
  return lc;
}

/* normalization of complex number. */
zComplex *zComplexNormalize(const zComplex *c, zComplex *nc)
{
  return zComplexDiv( c, _zComplexAbs( c ), nc );
}

/* I/O */

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
