/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_misc - miscellanies.
 */

#include <zm/zm_misc.h>

/* check if two values are equal. */
bool zEqual(double a, double b, double tol)
{
  return zIsTol( a-b, tol*zMax( 1, zMax( fabs(a), fabs(b) ) ) );
}

/* a pair of sine and cosine values of an angle */
void zSinCos(double angle, double *s, double *c)
{
  _zSinCos( angle, s, c );
}

/* normalize phase within the range from minus pi to pi. */
double zPhaseNormalize(double angle)
{
  double s, c;

  _zSinCos( angle, &s, &c );
  return atan2( s, c );
}

/* square a value. */
double zSqr(double x){ return _zSqr(x); }

/* cube a value. */
double zCube(double x){ return _zCube(x); }

/* a non-negative integer power. */
double zPowN(double x, unsigned int n)
{
  double ret = 1;

  while( n > 0 ){
    if( n % 2 == 0 ){
      x *= x;
      n >>= 1;
    } else{
      ret *= x;
      n--;
    }
  }
  return ret;
}

/* planer line connecting two points. */
double zLine(double x, double x0, double y0, double x1, double y1)
{
  return ( (y1-y0) * x - x0*y1-x1*y0 ) / ( x1 - x0 );
}

/* x-displacement of cycloid. */
double zCycloidX(double phase)
{
  return phase - sin( zPIx2 * phase ) / zPIx2;
}

/* y-displacement of cycloid. */
double zCycloidY(double phase)
{
  return 0.5 * ( 1 - cos( zPIx2 * phase ) );
}

/* sigmoid function. */
double zSigmoid(double x)
{
  return 1.0 / ( 1 + exp(-x) );
}

/* cubic root. */
double zCbrt(double x)
{
  double w;

  if( x == 0 ) return 0;
  w = ( x > 0 ) ? pow( x, 1.0/3.0 ) : -pow( -x, 1.0/3.0 );
  return 0.5 * ( w + 3.0*x / (2.0*w*w + x/w) );
}

/* permutation. */
double zPermutation(int n, int i)
{
  double result = 1.0;

  if( i > n || i < 0 ) return 0;
  for( ; i>0; i--, n-- ) result *= n;
  return result;
}

/* factorial. */
double zFactorial(int n)
{
  return zPermutation( n, n );
}

/* combination. */
double zCombination(int n, int i)
{
  double val;
  int j;

  if( n - i < i ) i = n - i;
  if( i == 0 ) return 1;
  for( val=1.0, j=0; j<i; j++ ){
    val *= (double)( n - j ) / ( i - j );
  }
  return val;
}

/* combination (recursive version). */
double zCombinationRecursive(int n, int i)
{
  if( n <= 0 ) return 1;
  if( i == 0 || i == n ) return 1;
  if( i < 0 || i > n ) return NAN;
  return zCombinationRecursive( n-1, i-1 ) + zCombinationRecursive( n-1, i );
}

/* series of combination. */
double *zCombinationSeries(uint n, size_t size, double c[])
{
  uint i;
  int j;

  if( n >= size ){
    ZRUNERROR( ZM_ERR_OUTOFRANGE );
    return NULL;
  }
  for( i=0; i<=n; i++ ){
    c[i] = 1.0;
    for( j=i-1; j>0; j-- )
      c[j] += c[j-1];
  }
  return c;
}

/* smoothstep function. */
double zSmoothStep(double x, uint order)
{
  uint k;
  double c, val;

  if( x <= 0 ) return 0;
  if( x >= 1 ) return 1;
  for( val=0, k=0; k<=order; k++ ){
    c = zCombination( order+k, k ) * zCombination( 2*order+1, order-k ) * pow( x, k );
    val += k % 2 ? -c : c;
  }
  return val * pow( x, order+1 );
}

/* derivative of smoothstep function. */
double zSmoothStepDif(double x, uint order)
{
  uint k;
  double c, val;

  if( x <= 0 || x >= 1 ) return 0;
  for( val=0, k=0; k<=order; k++ ){
    c = ( order + k + 1 ) * zCombination( order+k, k ) * zCombination( 2*order+1, order-k ) * pow( x, k );
    val += k % 2 ? -c : c;
  }
  return val * pow( x, order );
}
