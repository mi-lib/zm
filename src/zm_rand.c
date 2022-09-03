/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_rand - random number generator.
 */

#include <zm/zm_rand.h>

/* ********************************************************** */
/* a variety of random distributions
 * ********************************************************** */

/* gamma distribution family */

/* gamma distribution. */
double zRandGamma(zRandMT *mt, double a)
{
  double t, u, x, y;

  if( a > 1 ){
    a -= 1.0;
    t = sqrt(2*a+1);
    do{
      do{
        do{
          x = 1 - zRandMTN(mt);
          y = 2 * zRandMTN(mt) - 1;
        } while( x*x+y*y > 1 );
        y /= x;
        x = t*y + a;
      } while( x <= 0 );
      u = a * log(x/a) - t*y;
    } while( u<-50 || zRandMTN(mt) > (1+y*y)*exp(u) );
  } else{
    t = zE / ( a + zE );
    do{
      if( zRandMTN(mt) < t ){
        x = pow( zRandMTN(mt), 1/a );
        y = exp(-x);
      } else {
        x = 1 - log( 1 - zRandMTN(mt) );
        y = pow( x, a-1 );
      }
    } while( zRandMTN(mt) >= y );
  }
  return x;
}

/* Chi squared distribution. */
double zRandChiSqr(zRandMT *mt, double nu)
{
  return 2 * zRandGamma( mt, 0.5*nu );
}

/* beta distribution. */
double zRandBeta(zRandMT *mt, double a, double b)
{
  double x;

  x = zRandGamma(mt,a);
  return x / ( x + zRandGamma(mt,b) );
}

/* F distribution. */
double zRandFD(zRandMT *mt, double a, double b)
{
  double x1, x2;

  x1 = zRandChiSqr(mt,a);
  x2 = zRandChiSqr(mt,b);
  return (x1*b) / (x2*a);
}

/* normal distribution based on Box-Muller's method. */
double zRandND0(zRandMT *mt)
{
  double t, u;

  if( !mt ) mt = zRandMTDefault();
  if( !mt->nd_sw ){
    t = sqrt( -2*log( 1-zRandMTNU(mt) ) );
    u = zPIx2*zRandMTNU(mt);
    mt->nd_last = t * sin(u);
    mt->nd_sw = true;
    return t * cos(u);
  }
  mt->nd_sw = false;
  return mt->nd_last;
}

/* a random value generator based on the normal distribution based on Box-Muller's method. */
double zRandND(zRandMT *mt, double mu, double sigma)
{
  return mu + sigma * zRandND0( mt );
}

/* t distribution. */
double zRandT(zRandMT *mt, double n)
{
  double a, b, c;

  if( n <= 2 ){
    do{
      a = zRandChiSqr(mt,n);
    } while( zIsTiny(a) );
    return zRandND0(mt) / sqrt(a/n);
  }
  do{
    a = zRandND0(mt);
    b = a*a / (n-2);
    c = log( 1 - zRandMTNU(mt) ) / ( 1 - 0.5*n );
  } while( exp(-b-c) > 1-b );
  return a / sqrt( (1-2.0/n) * (1-b) );
}

/* continuous distribution */

/* exponential distribution with unit meanvalue. */
double zRandExp(zRandMT *mt)
{
  return -log( 1 - zRandMTNU(mt) );
}

/* logistic distribution. */
double zRandLog(zRandMT *mt)
{
  double r;

  do{
    r = zRandMTNU(mt);
  } while( zIsTiny(r) );
  return log( r / (1-r) );
}

/* power distribution. */
double zRandPower(zRandMT *mt, double n)
{
  return pow( zRandMTNU(mt), 1.0/(n+1) );
}

/* triangle distribution. */
double zRandTri(zRandMT *mt)
{
  return zRandMTN(mt) - zRandMTN(mt);
}

/* Cauchy's distribution. */
double zRandCauchy(zRandMT *mt)
{
  double x, y;

  do{
    x = 1 - zRandMTF(mt,0,1);
    y = 2 * zRandMTF(mt,0,1)-1;
  } while( x*x + y*y > 1 );
  return y/x;
}

/* Weibull distribution. */
double zRandWeibull(zRandMT *mt, double alpha)
{
  return pow( -log( 1 - zRandMTNU(mt) ), 1/alpha );
}

/* discrete distribution */

/* binomial distribution. */
int zRandBinom(zRandMT *mt, int n, double p)
{
  int i, r;

  for( r=0, i=0; i<n; i++ )
    if( zRandMTN(mt) < p ) r++;
  return r;
}

/* two-variable binomial distribution with relation. */
void zRandBinom2(zRandMT *mt, double r, double *x, double *y)
{
  double r1, r2, s;

  do{
    r1 = 2 * zRandMTNU(mt) - 1;
    r2 = 2 * zRandMTNU(mt) - 1;
    s = r1*r1 + r2*r2;
  } while( s > 1 || zIsTiny(s) );
  s = -log(s) / s;
  r1 = sqrt( (1+r)*s ) * r1;
  r2 = sqrt( (1-r)*s ) * r2;
  *x = r1+r2;
  *y = r1-r2;
}

/* Poisson's distribution. */
int zRandPoisson(zRandMT *mt, double lambda)
{
  int i;

  lambda = exp(lambda) * zRandMTF(mt,0,1);
  for( i=0; lambda>1; i++ )
    lambda *= zRandMTF(mt,0,1);
  return i;
}

/* geometric distribution. */
int zRandGeo(zRandMT *mt, double p)
{
  return (int)ceil( log(1-zRandMTNU(mt)) / log(1-p) );
}
