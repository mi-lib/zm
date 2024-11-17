/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_complex_arith - complex number class: arithmetics.
 */

#include <zm/zm_complex.h>

/* squared absolute value of a complex number. */
double zComplexSqrAbs(zComplex *c)
{
  return _zComplexSqrAbs( c );
}

/* absolute value of a complex number. */
double zComplexAbs(zComplex *c)
{
  return _zComplexAbs( c );
}

/* argument angle of a complex number. */
double zComplexArg(zComplex *c)
{
  return _zComplexArg( c );
}

/* complex conjugate. */
zComplex *zComplexConj(zComplex *c, zComplex *cc)
{
  _zComplexConj( c, cc );
  return cc;
}

/* add two complex numbers. */
zComplex *zComplexAdd(zComplex *c1, zComplex *c2, zComplex *c)
{
  _zComplexAdd( c1, c2, c );
  return c;
}

/* subtract complex number from the other. */
zComplex *zComplexSub(zComplex *c1, zComplex *c2, zComplex *c)
{
  _zComplexSub( c1, c2, c );
  return c;
}

/* reverse complex number. */
zComplex *zComplexRev(zComplex *c, zComplex *rc)
{
  _zComplexRev( c, rc );
  return rc;
}

/* multiply a complex number by a real number. */
zComplex *zComplexMul(zComplex *c, double k, zComplex *ec)
{
  _zComplexMul( c, k, ec );
  return ec;
}

/* divide a complex number by a real number. */
zComplex *zComplexDiv(zComplex *c, double k, zComplex *rc)
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
zComplex *zComplexCMul(zComplex *c1, zComplex *c2, zComplex *c)
{
  _zComplexCMul( c1, c2, c );
  return c;
}

/* multiply a complex numbers by the conjugate of another complex number. */
zComplex *zComplexCMulConj(zComplex *c1, zComplex *c2, zComplex *c)
{
  _zComplexCMulConj( c1, c2, c );
  return c;
}

/* invert a complex number. */
zComplex *zComplexInv(zComplex *c, zComplex *ic)
{
  double r2;
  zComplex cc;

  r2 = _zComplexSqrAbs( c );
  _zComplexConj( c, &cc );
  return zComplexDiv( &cc, r2, ic );
}

/* divide a complex number by another. */
zComplex *zComplexCDiv(zComplex *c1, zComplex *c2, zComplex *c)
{
  zComplex cc;

  zComplexInv( c2, &cc );
  _zComplexCMul( c1, &cc, c );
  return c;
}

/* multiply two complex numbers. */
zComplex *zComplexCMulDRC(zComplex *c1, zComplex *c2)
{
  zComplex tmp;

  _zComplexCMul( c1, c2, &tmp );
  _zComplexCopy( &tmp, c1 );
  return c1;
}

/* multiply a complex numbers by the conjugate of another complex number. */
zComplex *zComplexCMulConjDRC(zComplex *c1, zComplex *c2)
{
  zComplex tmp;

  _zComplexCMulConj( c1, c2, &tmp );
  _zComplexCopy( &tmp, c1 );
  return c1;
}

/* divide a complex number by another. */
zComplex *zComplexCDivDRC(zComplex *c1, zComplex *c2)
{
  zComplex tmp;

  zComplexCDiv( c1, c2, &tmp );
  _zComplexCopy( &tmp, c1 );
  return c1;
}

/* power a complex number by a real number. */
zComplex *zComplexPow(zComplex *c, double z, zComplex *pc)
{
  double r, theta;

  r = _zComplexAbs( c );
  theta = _zComplexArg( c );
  _zComplexCreatePolar( pc, pow( r, z ), theta * z );
  return pc;
}

/* power complex number by a real number. */
zComplex *zComplexPowRef(zComplex *c, double z, zComplex *ref, zComplex *pc)
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
zComplex *zComplexCPow(zComplex *c, zComplex *z, zComplex *pc)
{
  double r, theta;

  r = _zComplexAbs( c );
  theta = _zComplexArg( c );
  return zComplexCreatePolar( pc,
    pow( r, z->re ) * exp( -z->im * theta ), z->re*theta + z->im*log(r) );
}

/* the real number base logarithm of complex number. */
zComplex *zComplexLog(zComplex *c, double base, zComplex *lc)
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
zComplex *zComplexCLog(zComplex *c, zComplex *base, zComplex *lc)
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
zComplex *zComplexNormalize(zComplex *c, zComplex *nc)
{
  return zComplexDiv( c, _zComplexAbs( c ), nc );
}
