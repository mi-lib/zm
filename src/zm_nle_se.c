/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_nle_se - nonlinear equation: single nonlinear equation solver.
 */

#include <zm/zm_nle.h>

/* scale of the variable for the judgement of the step, which is the larger of |x1| and |x2|, or 1 if
 * both are zero. */
static double _zNLE1Scale(double x1, double x2)
{
  double s0;

  return ( s0 = _zMax( fabs( x1 ), fabs( x2 ) ) ) == 0 ? 1.0 : s0;
}

/* check if a step dx is small relatively to x, namely, |dx| <= tol_step ( |x| + s0 ). */
static bool _zNLE1StepIsSmall(double dx, double x, double s0, double tol_step)
{
  return fabs( dx ) <= tol_step * ( fabs( x ) + s0 );
}

/* check if a value of the function is valid, i.e., neither not-a-number nor infinity. */
static bool _zNLE1IsValid(double fx)
{
  return !zIsNan( fx ) && fabs( fx ) < HUGE_VAL;
}

/* evaluate the function at x. if the value is invalid, x is moved toward xo by halving the distance
 * until the value becomes valid. it returns false if the distance becomes small relatively to xo or
 * the number of evaluations count exceeds iter. */
static bool _zNLE1Eval(double (* f)(double,void*), void *util, double xo, double *x, double *fx, double s0, double tol_step, int iter, int *count)
{
  while( *count < iter ){
    (*count)++;
    if( _zNLE1IsValid( ( *fx = f( *x, util ) ) ) ) return true;
    *x = xo + 0.5 * ( *x - xo );
    if( _zNLE1StepIsSmall( *x - xo, xo, s0, tol_step ) ) return false;
  }
  return false;
}

/* search a range in which the function changes its sign, given a range [x1,x2] where the function has
 * the same sign at the both end points. the range is divided into 2^k equal subranges for k=1,2,...,
 * and the new dividing points are evaluated from the left. since the function has the same sign as
 * f1 at all the points evaluated before, the first point x with the opposite sign brackets a root with
 * the previous point evaluated at the same level (or x1 for the first point). points where the value
 * of the function is invalid are skipped.
 * it returns 1 if a root is found, which is stored to answer, 0 if a bracketing range is found, which
 * is stored to x1, x2, f1 and f2, or -1 if the subranges become small or the number of evaluations
 * exceeds iter, where the point with the smallest absolute value of the function is stored to answer.
 * count is the number of evaluations. */
static int _zNLE1Search(double (* f)(double,void*), double *answer, double *x1, double *x2, void *util, double *f1, double *f2, double s0, double tol_step, double tol_err, int iter, int *count)
{
  double d, h, x, fx, xp, fp, fmin;
  int n, j;

  *answer = fabs( *f1 ) < fabs( *f2 ) ? *x1 : *x2;
  fmin = _zMin( fabs( *f1 ), fabs( *f2 ) );
  d = *x2 - *x1;
  for( n=2; ; n*=2 ){
    if( _zNLE1StepIsSmall( ( h = d / n ), _zMax( fabs( *x1 ), fabs( *x2 ) ), s0, tol_step ) ) break;
    xp = *x1; fp = *f1;
    for( j=1; j<n; j+=2 ){
      if( *count >= iter ) goto FAILURE;
      x = *x1 + h * j;
      fx = f( x, util );
      (*count)++;
      if( !_zNLE1IsValid( fx ) ) continue; /* skip an invalid point */
      if( zIsTol( fx, tol_err ) ){
        *answer = x;
        return 1;
      }
      if( fx * *f1 < 0 ){
        *x1 = xp; *f1 = fp;
        *x2 = x;  *f2 = fx;
        return 0;
      }
      if( fabs( fx ) < fmin ){
        *answer = x;
        fmin = fabs( fx );
      }
      xp = x; fp = fx;
    }
  }
 FAILURE:
  ZRUNWARN( ZM_WARN_NLE_NOROOT );
  return -1;
}

/* check if a range [x1,x2] brackets a root, and search a bracketing range if not.
 * it returns 1 if a root is found, which is stored to answer, 0 if a bracketing range is found, which
 * is stored to x1, x2, f1 and f2, or -1 otherwise, including the case where the value of the function
 * is invalid at x1 or x2. count is the number of evaluations except at x1 and x2. */
static int _zNLE1Bracket(double (* f)(double,void*), double *answer, double *x1, double *x2, void *util, double *f1, double *f2, double s0, double tol_step, double tol_err, int iter, int *count)
{
  *count = 0;
  *answer = *x1;
  if( zIsTol( ( *f1 = f( *x1, util ) ), tol_err ) ) return 1;
  *answer = *x2;
  if( zIsTol( ( *f2 = f( *x2, util ) ), tol_err ) ) return 1;
  if( !_zNLE1IsValid( *f1 ) || !_zNLE1IsValid( *f2 ) ) return -1; /* invalid end points */
  if( *f1 * *f2 > 0 )
    return _zNLE1Search( f, answer, x1, x2, util, f1, f2, s0, tol_step, tol_err, iter, count );
  return 0;
}

/* the iteration terminates if the number of evaluations count exceeds iter. */
static int _zNLE1Fail(int count, int iter)
{
  if( count >= iter ) ZITERWARN( iter );
  return -1;
}

/* bisection method. */
int zNLE1Bisec(double (* f)(double,void*), double *answer, double x1, double x2, void *util, double tol_step, double tol_err, int iter)
{
  double xm, f1, f2, fm, s0;
  int count, ret;

  ZITERINIT( iter );
  s0 = _zNLE1Scale( x1, x2 );
  if( ( ret = _zNLE1Bracket( f, answer, &x1, &x2, util, &f1, &f2, s0, tol_step, tol_err, iter, &count ) ) != 0 )
    return ret > 0 ? count : -1;
  while( count < iter ){
    xm = 0.5 * ( x1 + x2 );
    if( !_zNLE1Eval( f, util, fabs(f1) < fabs(f2) ? x1 : x2, &xm, &fm, s0, tol_step, iter, &count ) ) break;
    *answer = xm;
    if( zIsTol( fm, tol_err ) ) return count; /* succeed */
    if( f1*fm < 0 ){
      x2 = xm; f2 = fm;
    } else{
      x1 = xm; f1 = fm;
    }
    if( _zNLE1StepIsSmall( x2 - x1, xm, s0, tol_step ) ) return -1; /* stagnate */
  }
  return _zNLE1Fail( count, iter );
}

/* secant method. */
int zNLE1Secant(double (* f)(double,void*), double *answer, double x1, double x2, void *util, double tol_step, double tol_err, int iter)
{
  double f1, f2, fc, xc, s0;
  int count = 0, ret;

  ZITERINIT( iter );
  s0 = _zNLE1Scale( x1, x2 );
  *answer = x1;
  if( zIsTol( ( f1 = f( x1, util ) ), tol_err ) ) return 0;
  *answer = x2;
  if( zIsTol( ( f2 = f( x2, util ) ), tol_err ) ) return 0;
  if( !_zNLE1IsValid( f1 ) || !_zNLE1IsValid( f2 ) ) return -1; /* invalid initial points */
  /* search a range in which the function changes its sign if the secant is horizontal */
  if( f1 == f2 &&
      ( ret = _zNLE1Search( f, answer, &x1, &x2, util, &f1, &f2, s0, tol_step, tol_err, iter, &count ) ) != 0 )
    return ret > 0 ? count : -1;
  while( count < iter ){
    if( f1 == f2 ) return -1; /* horizontal secant */
    xc = ( f2*x1 - f1*x2 ) / ( f2 - f1 );
    if( !_zNLE1Eval( f, util, x2, &xc, &fc, s0, tol_step, iter, &count ) ) break;
    *answer = xc;
    if( zIsTol( fc, tol_err ) ) return count; /* succeed */
    if( _zNLE1StepIsSmall( xc - x2, xc, s0, tol_step ) ) return -1; /* stagnate */
    x1 = x2; f1 = f2;
    x2 = xc; f2 = fc;
  }
  return _zNLE1Fail( count, iter );
}

/* Regula-Falsi method with the Illinois modification (M. Dowell and P. Jarratt, 1971). if the new
 * estimate replaces the latest estimate, i.e., the other end point is retained, the value of the
 * function at the retained end point is halved, so that it is replaced in the following iterations.
 * x2 is regarded as the latest estimate at the beginning. */
int zNLE1RegulaFalsi(double (* f)(double,void*), double *answer, double x1, double x2, void *util, double tol_step, double tol_err, int iter)
{
  double f1, f2, fc, xc, s0;
  int count, ret, side = 1; /* the end point of the latest estimate (-1 for x1, 1 for x2) */

  ZITERINIT( iter );
  s0 = _zNLE1Scale( x1, x2 );
  if( ( ret = _zNLE1Bracket( f, answer, &x1, &x2, util, &f1, &f2, s0, tol_step, tol_err, iter, &count ) ) != 0 )
    return ret > 0 ? count : -1;
  while( count < iter ){
    xc = ( f2*x1 - f1*x2 ) / ( f2 - f1 );
    if( !_zNLE1Eval( f, util, fabs(f1) < fabs(f2) ? x1 : x2, &xc, &fc, s0, tol_step, iter, &count ) ) break;
    *answer = xc;
    if( zIsTol( fc, tol_err ) ) return count; /* succeed */
    if( f1*fc < 0 ){
      x2 = xc; f2 = fc;
      if( side == 1 ) f1 *= 0.5; /* the latest estimate x2 is replaced */
      side = 1;
    } else{
      x1 = xc; f1 = fc;
      if( side == -1 ) f2 *= 0.5; /* the latest estimate x1 is replaced */
      side = -1;
    }
    if( _zNLE1StepIsSmall( x2 - x1, xc, s0, tol_step ) ) return -1; /* stagnate */
  }
  return _zNLE1Fail( count, iter );
}

/* Van Wijngaarden-Dekker-Brent method (R. P. Brent, 1973).
 * b is the best estimate of the root, c is the counterpart of the bracketing range [b,c], and a is the
 * previous value of b. the step is chosen by the inverse quadratic interpolation of a, b and c, or by
 * the linear interpolation of b and c if a = c. the interpolation is rejected and the bisection is taken
 * instead, if the interpolated point is not between b and the 3/4 point from b to c, or the step is not
 * less than a half of the step before last. the step is not less than the minimum step d_min. */
int zNLE1VDB(double (* f)(double,void*), double *answer, double x1, double x2, void *util, double tol_step, double tol_err, int iter)
{
  double a, b, c, fa, fb, fc, s0;
  double d, e, dmin, m, p, q, r, s;
  int count, ret;

  ZITERINIT( iter );
  s0 = _zNLE1Scale( x1, x2 );
  if( ( ret = _zNLE1Bracket( f, answer, &x1, &x2, util, &fa, &fb, s0, tol_step, tol_err, iter, &count ) ) != 0 )
    return ret > 0 ? count : -1;
  a = x1; b = x2;
  c = a; fc = fa;
  d = e = b - a;
  while( 1 ){
    if( fb * fc > 0 ){ /* keep [b,c] bracketing the root */
      c = a; fc = fa;
      d = e = b - a;
    }
    if( fabs( fc ) < fabs( fb ) ){ /* b is the best estimate */
      a = b; b = c; c = a;
      fa = fb; fb = fc; fc = fa;
    }
    *answer = b;
    if( zIsTol( fb, tol_err ) ) return count; /* succeed */
    dmin = 0.5 * tol_step * ( fabs( b ) + s0 );
    m = 0.5 * ( c - b );
    if( fabs( m ) <= dmin ) return -1; /* stagnate */
    if( count >= iter ) break;
    if( fabs( e ) >= dmin && fabs( fa ) > fabs( fb ) ){
      s = fb / fa;
      if( a == c ){ /* linear interpolation */
        p = 2 * m * s;
        q = 1 - s;
      } else{ /* inverse quadratic interpolation */
        q = fa / fc;
        r = fb / fc;
        p = s * ( 2 * m * q * ( q - r ) - ( b - a ) * ( r - 1 ) );
        q = ( q - 1 ) * ( r - 1 ) * ( s - 1 );
      }
      if( p > 0 ) q = -q; else p = -p; /* the step is p/q */
      s = e; e = d;
      if( 2 * p < 3 * m * q - fabs( dmin * q ) && p < fabs( 0.5 * s * q ) ){
        d = p / q; /* accept the interpolation */
      } else{
        d = e = m; /* bisection */
      }
    } else{
      d = e = m; /* bisection */
    }
    a = b; fa = fb;
    b += fabs( d ) > dmin ? d : ( m > 0 ? dmin : -dmin );
    if( !_zNLE1Eval( f, util, a, &b, &fb, s0, tol_step, iter, &count ) ) break;
  }
  return _zNLE1Fail( count, iter );
}
