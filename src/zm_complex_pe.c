/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_complex_pe - complex number class: polynomial equation solver.
 */

#include <zm/zm_complex.h>

/* solve quadratic equation. the coefficients are divided by the power of two not greater than
 * max(|b|, sqrt(|a||c|)), so that b^2 - 4ac neither overflows nor underflows. since the scaling is
 * exact and the answers q/a and c/q, where q = -( b + sgn(b) sqrt(b^2 - 4ac) )/2, are invariant to it,
 * the answers are not affected by the scaling unless they overflow or underflow by themselves. */
zComplex *zQESolve(double a, double b, double c, zComplex ans[])
{
  double u, d, q;

  if( a == 0 ){
    ZRUNERROR( ZM_ERR_PE_DEFLATED );
    return NULL;
  }
  if( b == 0 && c == 0 ){ /* a double zero answer */
    zComplexCreate( &ans[0], 0, 0 );
    zComplexCreate( &ans[1], 0, 0 );
    return ans;
  }
  u = zPow2Floor( zMax( fabs( b ), sqrt( fabs( a ) ) * sqrt( fabs( c ) ) ) );
  a /= u;
  b /= u;
  c /= u;
  d = b*b - 4*a*c;
  if( d >= 0 ){
    d = sqrt( d );
    q = -0.5 * ( b + ( b >= 0 ? d : -d ) );
    zComplexCreate( &ans[0], q/a, 0 );
    zComplexCreate( &ans[1], c/q, 0 );
  } else{
    a *= 2;
    d = sqrt( -d ) / a;
    b /= -a;
    zComplexCreate( &ans[0], b,  d );
    zComplexCreate( &ans[1], b, -d );
  }
  return ans;
}

/* scale a cubic equation coeff[3] x^3 + ... + coeff[0] = 0 with coeff[3] != 0, so that the
 * intermediate values of Cardano's formula, which are of the order of the sixth power of the answers,
 * neither overflow nor underflow. the variable is scaled as x = s y, where s is the power of two not
 * greater than an upper bound of the magnitudes of the answers max_k |coeff[3-k]/coeff[3]|^(1/k), and
 * coeff[k] is replaced by coeff[k] / ( 2^e s^(3-k) ), where 2^e is the power of two not greater than
 * |coeff[3]|. since it is computed by one multiplication by a power of two, it is exact unless the
 * result underflows. the scaled constant term is the product of the answers divided by s^3, which
 * underflows if the ratio of the magnitudes of the largest and the smallest answers exceeds about
 * 1e+150; then the smaller answers become inaccurate.
 * it returns s, by which the answers of the scaled equation are to be multiplied. */
static double _zCESolveScale(double coeff[])
{
  double r = 0, rk, s;
  int k, e, es;

  for( k=1; k<=3; k++ ) /* the k-th roots are taken before the division so as not to overflow */
    if( ( rk = pow( fabs( coeff[3-k] ), 1.0/k ) / pow( fabs( coeff[3] ), 1.0/k ) ) > r ) r = rk;
  s = r > 0 && zIsFinite( r ) ? zPow2Floor( r ) : 1;
  frexp( s, &es ); /* s = 2^(es-1) */
  frexp( coeff[3], &e ); /* 2^(e-1) <= |coeff[3]| < 2^e */
  for( k=0; k<=3; k++ )
    coeff[k] = ldexp( coeff[k], -( e - 1 ) - ( 3 - k ) * ( es - 1 ) );
  return s;
}

/* deflate a cubic equation x^3 + 3 b x^2 + c x + d = 0 by a real answer r, and solve the remaining
 * quadratic equation x^2 + e x + f = 0, where the coefficients are computed backward from the constant
 * term, i.e., f = -d/r and e = ( f - c )/r. the backward deflation is stable if r is the answer of
 * the largest magnitude. */
static void _zCESolveDeflate(double c, double d, double r, zComplex ans[])
{
  double f;

  if( r == 0 ){ /* all answers are zero */
    zComplexCreate( &ans[0], 0, 0 );
    zComplexCreate( &ans[1], 0, 0 );
    return;
  }
  f = -d / r;
  zQESolve( 1, ( f - c ) / r, f, ans );
}

/* solve a cubic equation by Cardano's formula. only the answers of the largest magnitude are directly
 * computed by the formula, since the others might suffer from cancellation, and the remaining answers
 * are computed by deflating the cubic equation. */
zComplex *zCESolve(double a, double b, double c, double d, zComplex ans[])
{
  double coeff[4], s, g, p, q, w, z1, z2, r, m2, x[3];
  int i, imax;

  if( a == 0 ){
    ZRUNERROR( ZM_ERR_PE_DEFLATED );
    return NULL;
  }
  coeff[0] = d; coeff[1] = c; coeff[2] = b; coeff[3] = a;
  s = _zCESolveScale( coeff );
  d = coeff[0]; c = coeff[1]; b = coeff[2]; a = coeff[3];
  b /= 3 * a;
  c /= a;
  d /= a;
  p = c/3 - b*b;
  q = -b*(2*b*b-c) - d;
  g = q*q + 4*p*p*p;
  if( g > 0 ){ /* a real answer and a pair of complex conjugate answers */
    w = sqrt( g );
    /* z1 z2 = -p, where z1 is computed without cancellation, which happens in q - w when p ~ 0. */
    z1 = zCbrt( 0.5 * ( q >= 0 ? q + w : q - w ) );
    z2 = -p / z1;
    r = z1 + z2 - b;
    zComplexCreate( &ans[1], -0.5*(z1+z2)-b, 0.5*sqrt(3)*fabs(z1-z2) );
    m2 = zSqr( ans[1].re ) + zSqr( ans[1].im );
    if( r*r >= m2 ){ /* the real answer is the largest */
      _zCESolveDeflate( c, d, r, &ans[1] );
    } else{ /* the complex answers are the largest, and r is from the product of the answers */
      r = -d / m2;
      zComplexConj( &ans[1], &ans[2] );
    }
  } else{ /* three real answers, where p is necessarily negative */
    z1 = ( q != 0 ) ? atan2( sqrt(-g), q ) : zPI_2;
    w = 2 * sqrt( -p );
    x[0] =  w * cos(     z1 /3) - b;
    x[1] = -w * cos((zPI-z1)/3) - b;
    x[2] = -w * cos((zPI+z1)/3) - b;
    for( imax=0, i=1; i<3; i++ )
      if( fabs( x[i] ) > fabs( x[imax] ) ) imax = i;
    _zCESolveDeflate( c, d, ( r = x[imax] ), &ans[1] );
  }
  zComplexCreate( &ans[0], r, 0 );
  for( i=0; i<3; i++ ){
    ans[i].re *= s;
    ans[i].im *= s;
  }
  return ans;
}
