/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_complex_pe - complex number class: polynomial equation solver.
 */

#include <zm/zm_complex.h>

/* zQESolve
 * - solve quadratic equation.
 */
zComplex *zQESolve(double a, double b, double c, zComplex ans[])
{
  double d, q;

  if( a == 0 ){
    ZRUNERROR( ZM_ERR_PE_DEFL );
    return NULL;
  }
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

/* zCardano
 * - solve cubic equation by Cardano's formula: algebric roots
 *   of cubic equation, originally derived by Fontana=Tartaglia.
 */
zComplex *zCardano(double a, double b, double c, double d, zComplex ans[])
{
  double g, p, q, w, z1, z2;

  if( a == 0 ){
    ZRUNERROR( ZM_ERR_PE_DEFL );
    return NULL;
  }
  b /= 3 * a;
  c /= a;
  d /= a;
  p = c/3 - b*b;
  q = -b*(2*b*b-c) - d;
  g = q*q + 4*p*p*p;
  if( g > 0 ){
    w = sqrt( g );
    z1 = zCbrt( 0.5 * ( q + w ) );
    z2 = zCbrt( 0.5 * ( q - w ) );
    zComplexCreate( &ans[0], z1 + z2 - b, 0 );
    zComplexCreate( &ans[1], -0.5*(z1+z2)-b, 0.5*sqrt(3)*(z1-z2) );
    zComplexConj( &ans[1], &ans[2] );
  } else{
    z1 = ( q != 0 ) ? atan2( sqrt(-g), q ) : zPI_2;
    w = 2 * sqrt( -p ); /* 'p' necessarily be negative. */
    zComplexCreate( &ans[0],  w * cos(     z1 /3) - b, 0 );
    zComplexCreate( &ans[1], -w * cos((zPI-z1)/3) - b, 0 );
    zComplexCreate( &ans[2], -w * cos((zPI+z1)/3) - b, 0 );
  }
  return ans;
}
