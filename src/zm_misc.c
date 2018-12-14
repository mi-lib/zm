/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_misc - miscellanies.
 */

#include <zm/zm_misc.h>

/* zIsSgnOpp
 * - check if signs of two values are opposite.
 */
bool zIsSgnOpp(double a, double b)
{
  return ( a < 0 && b > 0 ) || ( b < 0 && a > 0 );
}

/* zSinCos
 * - compute sine and cosine values of an angle.
 */
void zSinCos(double angle, double *s, double *c)
{
  _zSinCos( angle, s, c );
}

/* zPhaseNormalize
 * - normalize phase within the range from minus pi to pi.
 */
double zPhaseNormalize(double angle)
{
  double s, c;

  _zSinCos( angle, &s, &c );
  return atan2( s, c );
}

/* zSqr
 * - squared value.
 */
double zSqr(double x){ return x * x; }

/* zCube
 * - cubic value.
 */
double zCube(double x){ return x * x * x; }

/* zPowN
 * - a non-negative integer power.
 */
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

/* zLine
 * - planer line connecting two points.
 */
double zLine(double x, double x0, double y0, double x1, double y1)
{
  return ( (y1-y0) * x - x0*y1-x1*y0 ) / ( x1 - x0 );
}

/* zCycloidX
 * - x-displacement of cycloid.
 */
double zCycloidX(double from, double to, double phase)
{
  return ( to - from ) * ( phase -sin( zPIx2 * phase ) / zPIx2 ) + from;
}

/* zCycloidY
 * - y-displacement of cycloid.
 */
double zCycloidY(double base, double h, double phase)
{
  return 0.5 * h * ( 1 - cos( zPIx2 * phase ) ) + base;
}

/* zSigmoid
 * - sigmoid function.
 */
double zSigmoid(double x, double s)
{
  return 1.0 / ( 1 + exp(-x/s) );
}

/* zCbrt
 * - cubic root.
 */
double zCbrt(double x)
{
  double w;

  if( x == 0 ) return 0;
  w = ( x > 0 ) ? pow( x, 1.0/3.0 ) : -pow( -x, 1.0/3.0 );
  return 0.5 * ( w + 3.0*x / (2.0*w*w + x/w) );
}
