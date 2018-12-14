/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_sf_gamma - special functions : gamma and beta function.
 */

#include <zm/zm_sf.h>

static double _zGammaSeed(double x);
static double _zGammaBias(double x);

/* (static)
 * _zGammaSeed
 * - a seed value of gamma function. see 'zGamma()' and 'zLnGamma()'.
 */
double _zGammaSeed(double x)
{
  double s = 1.000000000190015;

  s -=  0.000005395239384953 / (x+6);
  s +=  0.001208650973866179 / (x+5);
  s -=  1.231739572450155    / (x+4);
  s += 24.01409824083091     / (x+3);
  s -= 86.50532032941677     / (x+2);
  s += 76.18009172947146     / (x+1);
  return 2.5066282746310005 * s / x;
}

/* (static)
 * _zGammaBias
 * - a vias value of gamma function. see 'zGamma()' and 'zLnGamma()'.
 */
double _zGammaBias(double x)
{
  return ( x + 0.5 ) * log( x + 5.5 ) - x - 5.5;
}

/* zGamma
 * - gamma function.
 */
double zGamma(double x)
{
  return exp( _zGammaBias(x) ) * _zGammaSeed(x);
}

/* zLnGamma
 * - logarithm of gamma function.
 */
double zLnGamma(double x)
{
  return _zGammaBias(x) + log( _zGammaSeed(x) );
}

/* zBeta
 * - beta function.
 */
double zBeta(double z, double w)
{
  return exp( zLnGamma(z) + zLnGamma(w) - zLnGamma(z+w) );
}
