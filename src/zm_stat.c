/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_stat - statistics.
 */

#include <zm/zm_stat.h>
#include <zm/zm_sf.h>

/* zPermut
 * - permutation.
 */
double zPermut(int n, int i)
{
  double result = 1.0;

  if( i > n || i < 0 ) return 0;
  for( ; i>0; i--, n-- ) result *= n;
  return result;
}

/* zFacto
 * - factorial.
 */
double zFacto(int n)
{
  return zPermut( n, n );
}

/* zCombi
 * - combination.
 */
double zCombi(int n, int i)
{
  if( n-i < i ) i = n-i;
  return zPermut( n, i ) / zFacto( i );
}

/* zCombiSeries
 * - series of combination.
 */
double *zCombiSeries(int n, size_t size, double c[])
{
  register int i, j;

  if( n < 0 || n >= size ){
    ZRUNERROR( ZM_ERR_STAT_ILLS );
    return NULL;
  }
  for( i=0; i<=n; i++ ){
    c[0] = c[i] = 1.0;
    for( j=i-1; j>0; j-- )
      c[j] += c[j-1];
  }
  return c;
}

/* basic distribution functions */

/* zNormalDistrib
 * - normal distribution.
 */
double zNormalDistrib(double x, double mu, double sigma)
{
  return zND / sigma * exp( -0.5 * zSqr( ( x - mu ) / sigma ) );
}

/* zNormalCumDistrib
 * - normal cumulative distribution.
 */
#define Z_ND_CUM_MAX 200
double zNormalCumDistrib(double x, double mu, double sigma)
{
  register int i;
  double x2, p, pp, t;

  x -= mu;
  x2 = zSqr( x / sigma );
  t = p = pp = x * exp(-0.5*x2) / ( sqrt(zPIx2) * sigma );
  for( i=3; i<Z_ND_CUM_MAX; pp=p, i+=2 )
    if( ( p += ( t *= x2 / i ) ) == pp ) return 0.5 + p;
  return x > 0 ? 1.0 : 0.0;
}

/* zPoissonDistrib
 * - Poisson distribution.
 */
double zPoissonDistrib(int x, double lambda)
{
  return exp( -lambda ) * pow( lambda, x ) / zFacto( x );
}

/* zBinDistrib
 * - binomial distribution.
 */
double zBinDistrib(int x, int n, double p)
{
  return zCombi( n, x ) * pow( p, x ) * pow( 1-p, n-x );
}

/* zChi2Distrib
 * - chi-squared distribution.
 */
double zChi2Distrib(double x, int k)
{
  double kd;

  kd = 0.5 * k;
  return pow(0.5,kd)*pow(x,kd-1)*exp(-0.5*x)/zGamma(kd);
}

/* zChi2CumDistrib
 * - chi-squared cumulative distribution.
 */
double zChi2CumDistrib(double x, int k)
{
  register int i;
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

/* zDataMax
 * - maximum component of data.
 */
double zDataMax(double *data, int num, int *im)
{
  register int i;
  double max;

  if( im ) *im = 0;
  for( max = data[0], i=1; i<num; i++ )
    if( data[i] > max ){
      max = data[i];
      if( im ) *im = i;
    }
  return max;
}

/* zDataMin
 * - minimum component of data.
 */
double zDataMin(double *data, int num, int *im)
{
  register int i;
  double min;

  if( im ) *im = 0;
  for( min = data[0], i=1; i<num; i++ )
    if( data[i] < min ){
      min = data[i];
      if( im ) *im = i;
    }
  return min;
}

/* zDataAbsMax
 * - maximum absolute component of data.
 */
double zDataAbsMax(double *data, int num, int *im)
{
  register int i;
  double val, max;

  if( im ) *im = 0;
  for( max = fabs( data[0] ), i=1; i<num; i++ )
    if( ( val = fabs( data[i] ) ) > max ){
      max = val;
      if( im ) *im = i;
    }
  return max;
}

/* zDataAbsMin
 * - minimum component of data.
 */
double zDataAbsMin(double *data, int num, int *im)
{
  register int i;
  double val, min;

  if( im ) *im = 0;
  for( min = fabs( data[0] ), i=1; i<num; i++ )
    if( ( val = fabs( data[i] ) ) < min ){
      min = val;
      if( im ) *im = i;
    }
  return min;
}

/* zDataSum
 * - sum up all data values.
 */
double zDataSum(double *data, int num)
{
  double s=0, s_prev=0, q=0, r;
  register int i;

  for( i=0; i<num; i++ ){
    s = s_prev + data[i];
    r = s - s_prev;
    q += data[i] - r;
    s_prev = s;
  }
  return s + q;
}

/* zDataAve
 * - calculate the average of data.
 */
double zDataAve(double *data, int num)
{
  return zDataSum(data,num) / num;
}

/* zDataVar
 * - calculate the variance of data.
 */
double zDataVar(double *data, int num)
{
  register int i;
  double ave, result;

  ave = zDataAve( data, num );
  for( result=0, i=0; i<num; i++ )
    result += zSqr( data[i] - ave );
  return result / num;
}

/* zDataSD
 * - calculate the standard deviation of data.
 */
double zDataSD(double *data, int num)
{
  return sqrt( zDataVar( data, num ) );
}
