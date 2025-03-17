/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_stat - statistics.
 */

#include <zm/zm_stat.h>
#include <zm/zm_sf.h>

/* basic distribution functions */

/* normal distribution. */
double zNormalDistrib(double x, double mu, double sigma)
{
  return zND / sigma * exp( -0.5 * zSqr( ( x - mu ) / sigma ) );
}

/* normal cumulative distribution. */
#define Z_ND_CUM_MAX 200
double zNormalCumDistrib(double x, double mu, double sigma)
{
  int i;
  double x2, p, pp, t;

  x -= mu;
  x2 = zSqr( x / sigma );
  t = p = pp = x * exp(-0.5*x2) / ( sqrt(zPIx2) * sigma );
  for( i=3; i<Z_ND_CUM_MAX; pp=p, i+=2 )
    if( ( p += ( t *= x2 / i ) ) == pp ) return 0.5 + p;
  return x > 0 ? 1.0 : 0.0;
}

/* Poisson distribution. */
double zPoissonDistrib(int x, double lambda)
{
  return exp( -lambda ) * pow( lambda, x ) / zFactorial( x );
}

/* binomial distribution. */
double zBinDistrib(int x, int n, double p)
{
  return zCombination( n, x ) * pow( p, x ) * pow( 1-p, n-x );
}

/* chi-squared distribution. */
double zChi2Distrib(double x, int k)
{
  double kd;

  kd = 0.5 * k;
  return pow(0.5,kd)*pow(x,kd-1)*exp(-0.5*x)/zGamma(kd);
}

/* chi-squared cumulative distribution. */
double zChi2CumDistrib(double x, int k)
{
  int i;
  double s, t, chi;

  if( k == 1 )
    return 2 * zNormalCumDistrib(sqrt(x),0,1) - 1;
  if( k % 2 == 1 ){
    s = t = ( chi = sqrt(x) ) * exp(-0.5*x) / sqrt(zPIx2);
    for( i=3; i<k; i+=2 )
      s += ( t *= x / i );
    return 2 * zNormalCumDistrib(chi,0,1) - 1 - 2 * s;
  }
  s = t = exp(-0.5*x);
  for( i=2; i<k; i+=2 )
    s += ( t *= x / i );
  return 1 - s;
}

/* basic statistics computation */

/* maximum component of data. */
double zDataMax(const double *data, int size, int *im)
{
  int i, __im;
  double max;

  if( !im ) im = &__im;
  for( *im=0, max=data[0], i=1; i<size; i++ )
    if( data[i] > max ){
      max = data[i];
      *im = i;
    }
  return max;
}

/* minimum component of data. */
double zDataMin(const double *data, int size, int *im)
{
  int i, __im;
  double min;

  if( !im ) im = &__im;
  for( *im=0, min=data[0], i=1; i<size; i++ )
    if( data[i] < min ){
      min = data[i];
      *im = i;
    }
  return min;
}

/* minimum and maximum elements of data. */
void zDataMinMax(const double *data, int size, double *min, int *imin, double *max, int *imax)
{
  double __min, __max;
  int __imin, __imax;
  int i;

  if( !min ) min = &__min;
  if( !max ) max = &__max;
  if( !imin ) imin = &__imin;
  if( !imax ) imax = &__imax;
  *min = *max = data[( *imin = *imax = 0 )];
  for( i=0; i<size; i++ ){
    if( data[i] < *min ){
      *imin = i;
      *min = data[i];
    }
    if( data[i] > *max ){
      *imax = i;
      *max = data[i];
    }
  }
}

/* maximum absolute component of data. */
double zDataAbsMax(const double *data, int size, int *im)
{
  int i, __im;
  double val, max;

  if( !im ) im = &__im;
  for( *im=0, max=fabs( data[0] ), i=1; i<size; i++ )
    if( ( val = fabs( data[i] ) ) > max ){
      max = val;
      *im = i;
    }
  return max;
}

/* minimum absolute component of data. */
double zDataAbsMin(const double *data, int size, int *im)
{
  int i, __im;
  double val, min;

  if( !im ) im = &__im;
  for( *im=0, min=fabs( data[0] ), i=1; i<size; i++ )
    if( ( val = fabs( data[i] ) ) < min ){
      min = val;
      *im = i;
    }
  return min;
}

/* sum up all values of data. */
double zDataSum(const double *data, int size)
{
  double s=0, s_prev=0, q=0, r;
  int i;

  for( i=0; i<size; i++ ){
    s = s_prev + data[i];
    r = s - s_prev;
    q += data[i] - r;
    s_prev = s;
  }
  return s + q;
}

/* mean of all values of data. */
double zDataMean(const double *data, int size)
{
  return zDataSum( data, size ) / size;
}

/* calculate the variance of data. */
double zDataVar(const double *data, int size)
{
  int i;
  double mean, result;

  mean = zDataMean( data, size );
  for( result=0, i=0; i<size; i++ )
    result += zSqr( data[i] - mean );
  return result / size;
}

/* calculate the standard deviation of data. */
double zDataStandardDeviation(const double *data, int size)
{
  return sqrt( zDataVar( data, size ) );
}

/* check if a value is a member of data. */
bool zDataIsIncluded(const double *data, int size, double val, double tol)
{
  int i;

  for( i=0; i<size; i++ )
    if( zEqual( data[i], val, tol ) ) return true;
  return false;
}
