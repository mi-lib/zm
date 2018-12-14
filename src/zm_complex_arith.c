/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_complex_arith - complex number class: arithmetics.
 */

#include <zm/zm_complex.h>

/* zComplexSqrAbs
 * - absolute value of complex.
 */
double zComplexSqrAbs(zComplex *c)
{
  return zSqr(c->re) + zSqr(c->im);
}

/* zComplexArg
 * - argument angle of complex.
 */
double zComplexArg(zComplex *c)
{
  return atan2( c->im, c->re );
}

/* zComplexConj
 * - complex conjugate.
 */
zComplex *zComplexConj(zComplex *c, zComplex *cc)
{
  return zComplexCreate( cc, c->re,-c->im );
}

/* zComplexAdd
 * - add two complex numbers.
 */
zComplex *zComplexAdd(zComplex *c1, zComplex *c2, zComplex *c)
{
  return zComplexCreate( c, c1->re+c2->re, c1->im+c2->im );
}

/* zComplexSub
 * - subtract complex number from the other.
 */
zComplex *zComplexSub(zComplex *c1, zComplex *c2, zComplex *c)
{
  return zComplexCreate( c, c1->re-c2->re, c1->im-c2->im );
}

/* zComplexRev
 * - reverse complex number.
 */
zComplex *zComplexRev(zComplex *c, zComplex *rc)
{
  return zComplexCreate( rc, -c->re, -c->im );
}

/* zComplexMul
 * - multiply a complex number by a real number.
 */
zComplex *zComplexMul(zComplex *c, double k, zComplex *ec)
{
  return zComplexCreate( ec, c->re*k, c->im*k );
}

/* zComplexDiv
 * - divide a complex number by a real number.
 */
zComplex *zComplexDiv(zComplex *c, double k, zComplex *rc)
{
  if( zIsTiny( k ) ){
    ZRUNERROR( ZM_ERR_ZERODIV );
    return NULL;
  }
  return zComplexCreate( rc, c->re/k, c->im/k );
}

/* zComplexCMul
 * - multiply two complex numbers.
 */
zComplex *zComplexCMul(zComplex *c1, zComplex *c2, zComplex *c)
{
  return zComplexCreate( c,
    c1->re*c2->re - c1->im*c2->im, c1->re*c2->im + c1->im*c2->re );
}

/* zComplexCDiv
 * - divide complex number by the other.
 */
zComplex *zComplexCDiv(zComplex *c1, zComplex *c2, zComplex *c)
{
  double r;
  zComplex cc;

  r = zComplexAbs( c2 );
  zComplexConj( c2, &cc );
  zComplexDiv( &cc, r, &cc );
  zComplexCMul( c1, &cc, c );
  return zComplexDiv( c, r, c );
}

/* zComplexPow
 * - power complex number by a real number.
 */
zComplex *zComplexPow(zComplex *c, double z, zComplex *pc)
{
  double r, theta;

  r = zComplexAbs( c );
  theta = zComplexArg( c );
  return zComplexPolar( pc, pow( r, z ), theta * z );
}

/* zComplexPowRef
 * - power complex number by a real number.
 */
zComplex *zComplexPowRef(zComplex *c, double z, zComplex *ref, zComplex *pc)
{
  double r, theta, w;

  r = pow( zComplexAbs(c), z );
  theta = zComplexArg(c) * z;
  w = 2 * zPI * z;
  if( ref && !zIsTiny(z) )
    theta += w * ceil( ( zComplexArg(ref) - theta ) / w );
  return zComplexPolar( pc, r, theta );
}

/* zComplexCPow
 * - power complex number by another complex number.
 */
zComplex *zComplexCPow(zComplex *c, zComplex *z, zComplex *pc)
{
  double r, theta;

  r = zComplexAbs( c );
  theta = zComplexArg( c );
  return zComplexPolar( pc,
    pow( r, z->re ) * exp( -z->im * theta ),
    z->re*theta + z->im*log(r) );
}

/* zComplexLog
 * - the real number base logarithm of complex number.
 */
zComplex *zComplexLog(zComplex *c, double base, zComplex *lc)
{
  double r, theta, d;

  if( base <= 0 || base == 1.0 ){
    ZRUNERROR( ZM_ERR_LOG_INVALID );
    return NULL;
  }
  r = zComplexAbs( c );
  theta = zComplexArg( c );
  d = log(base);
  return zComplexCreate( lc, log(r)/d, theta/d );
}

/* zComplexCLog
 * - the complex number base logarithm of complex number.
 */
zComplex *zComplexCLog(zComplex *c, zComplex *base, zComplex *lc)
{
  double rb, thetab, lrb, d, rc, thetac, lrc;

  rb = zComplexAbs( base );
  thetab = zComplexArg( base );
  if( rb == 0 ){
    ZRUNERROR( ZM_ERR_LOG_INVALID );
    return NULL;
  }
  lrb = log( rb );
  d = lrb*lrb + thetab*thetab;
  if( d == 0 ){
    ZRUNERROR( ZM_ERR_LOG_INVALID );
    return NULL;
  }
  rc = zComplexAbs( c );
  thetac = zComplexArg( c );
  lrc = log( rc );
  return zComplexCreate( lc,
    (lrb*lrc+thetab*thetac)/d, (thetac*lrb-thetab*lrc)/d );
}

/* zComplexNormalize
 * - normalization of complex number.
 */
zComplex *zComplexNormalize(zComplex *c, zComplex *nc)
{
  return zComplexDiv( c, zComplexAbs( c ), nc );
}
