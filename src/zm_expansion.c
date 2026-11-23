/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_expansion - arithmetics of expansions.
 */

#include <zm/zm_expansion.h>

/* add a value to an expansion, eliminating zero components (Shewchuk's Grow-Expansion).
 * since h[k] is written after e[i] is read with k <= i, h can be the same with e. */
int zGrowExpansion(int n, const double e[], double b, double h[])
{
  double q, hh;
  int i, k = 0;

  for( q=b, i=0; i<n; i++ ){
    q = zAdd( q, e[i], &hh );
    if( hh != 0 ) h[k++] = hh;
  }
  if( q != 0 || k == 0 ) h[k++] = q;
  return k;
}

/* multiply an expansion by a value, eliminating zero components (Shewchuk's Scale-Expansion).
 * |p1| >= |sum| holds in the loop (see Theorem 19 of the paper). */
int zScaleExpansion(int n, const double e[], double b, double h[])
{
  double q, sum, hh, p1, p0;
  int i, k = 0;

  if( n <= 0 ){
    h[0] = 0;
    return 1;
  }
  q = zMul( e[0], b, &hh );
  if( hh != 0 ) h[k++] = hh;
  for( i=1; i<n; i++ ){
    p1 = zMul( e[i], b, &p0 );
    sum = zAdd( q, p0, &hh );
    if( hh != 0 ) h[k++] = hh;
    q = zAddFast( p1, sum, &hh );
    if( hh != 0 ) h[k++] = hh;
  }
  if( q != 0 || k == 0 ) h[k++] = q;
  return k;
}

/* add a product of three values to an expansion. the product is an expansion with at most four
 * components, each of which is added to the expansion by zGrowExpansion() in place. */
int zAddExpansionTripleProd(int n, double e[], double a, double b, double c)
{
  double ab[2], abc[4];
  int i, k;

  ab[1] = zMul( a, b, &ab[0] );
  k = zScaleExpansion( 2, ab, c, abc );
  for( i=0; i<k; i++ )
    n = zGrowExpansion( n, e, abc[i], e );
  return n;
}
