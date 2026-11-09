#include <zm/zm_nle.h>

/* test equations of a single variable, where util points to the counter of evaluations */
#define count_up(util) do{ if( util ) (*(int *)(util))++; } while(0)

/* exp(-x) - x = 0, root: 0.567143290409783873 (omega constant) */
double omega_f(double x, void *util)
{
  count_up( util );
  return exp( -x ) - x;
}

/* x^3 - 2 x - 5 = 0 (the example by Wallis and Newton), root: 2.09455148154232659 */
double wallis_f(double x, void *util)
{
  count_up( util );
  return ( x*x - 2 )*x - 5;
}

/* tan(x) - x = 0 on [4,4.6], root: 4.49340945790906418 */
double tan_f(double x, void *util)
{
  count_up( util );
  return tan( x ) - x;
}

/* x^2 + 1 = 0, no root */
double noroot_f(double x, void *util)
{
  count_up( util );
  return x*x + 1;
}

typedef struct{
  const char *name;
  double (* f)(double,void*);
  double x1, x2;
  double root;
} problem_t;

problem_t problem[] = {
  { "omega",  omega_f, -10.0, 10.0, 0.567143290409783873 },
  { "Wallis", wallis_f,  2.0,  3.0, 2.09455148154232659 },
  { "tan",    tan_f,     4.0,  4.6, 4.49340945790906418 },
  { NULL, NULL, 0, 0, 0 },
};

typedef struct{
  const char *label;
  int (* solve)(double (*)(double,void*),double*,double,double,void*,double,double,int);
  bool bracketing; /* true if a range that brackets a root is required */
  bool robust; /* true if the method is expected to solve all problems */
} method_t;

method_t method[] = {
  { "bisection",    zNLE1Bisec,       true,  true },
  { "secant",       zNLE1Secant,      false, true },
  { "Regula-Falsi", zNLE1RegulaFalsi, true,  true },
  { "VDB",          zNLE1VDB,         true,  true },
  { NULL, NULL, false, false },
};

#define TOL_STEP zTOL
#define TOL_ERR  1.0e-10

/* every method has to succeed only if the residual becomes small, robust methods have to solve all
 * problems, and the returned value has to be the number of evaluations except the end points. */
void assert_nle1_solve(void)
{
  problem_t *prob;
  method_t *meth;
  double x;
  int ret, count;
  bool result1 = true, result2 = true, result3 = true;

  zEchoOff();
  for( prob=problem; prob->name; prob++ ){
    for( meth=method; meth->label; meth++ ){
      count = 0;
      ret = meth->solve( prob->f, &x, prob->x1, prob->x2, &count, TOL_STEP, TOL_ERR, 0 );
      if( ret >= 0 && !zIsTol( prob->f( x, NULL ), TOL_ERR ) ){
        eprintf( "%s, %s: succeeded with f = %g\n", prob->name, meth->label, prob->f( x, NULL ) );
        result1 = false;
      }
      if( meth->robust && ( ret < 0 || !zIsTol( x - prob->root, 1.0e-8 ) ) ){
        eprintf( "%s, %s: returned %d, x = %.16g\n", prob->name, meth->label, ret, x );
        result2 = false;
      }
      if( ret >= 0 && ret != count - 2 ){
        eprintf( "%s, %s: returned %d for %d evaluations\n", prob->name, meth->label, ret, count );
        result3 = false;
      }
    }
  }
  zEchoOn();
  zAssert( zNLE1 (success only with small residual), result1 );
  zAssert( zNLE1 (solution of test problems), result2 );
  zAssert( zNLE1 (returned number of iterations), result3 );
}

/* the returned value has to be 0 if an end point is a root, and -1 if the range does not bracket a
 * root or the number of iterations exceeds the limit. */
void assert_nle1_return(void)
{
  method_t *meth;
  double x;
  int ret;
  bool result1 = true, result2 = true, result3 = true;

  zEchoOff();
  for( meth=method; meth->label; meth++ ){
    /* end points */
    x = HUGE_VAL;
    if( ( ret = meth->solve( wallis_f, &x, problem[1].root, 3.0, NULL, TOL_STEP, TOL_ERR, 0 ) ) != 0 || x != problem[1].root ){
      eprintf( "%s: returned %d, x = %g for the first end point at the root\n", meth->label, ret, x );
      result1 = false;
    }
    x = HUGE_VAL;
    if( ( ret = meth->solve( wallis_f, &x, 2.0, problem[1].root, NULL, TOL_STEP, TOL_ERR, 0 ) ) != 0 || x != problem[1].root ){
      eprintf( "%s: returned %d, x = %g for the second end point at the root\n", meth->label, ret, x );
      result1 = false;
    }
    /* no root */
    if( ( ret = meth->solve( noroot_f, &x, -1.0, 2.0, NULL, TOL_STEP, TOL_ERR, 0 ) ) != -1 ){
      eprintf( "%s: returned %d for an equation without root\n", meth->label, ret );
      result2 = false;
    }
    /* range without a sign change */
    if( meth->bracketing &&
        ( ret = meth->solve( wallis_f, &x, 3.0, 4.0, NULL, TOL_STEP, TOL_ERR, 0 ) ) != -1 ){
      eprintf( "%s: returned %d for a range without a sign change\n", meth->label, ret );
      result2 = false;
    }
    /* excess of iterations */
    if( ( ret = meth->solve( omega_f, &x, -10.0, 10.0, NULL, TOL_STEP, TOL_ERR, 2 ) ) != -1 ){
      eprintf( "%s: returned %d within 2 iterations\n", meth->label, ret );
      result3 = false;
    }
  }
  zEchoOn();
  zAssert( zNLE1 (returned value for an end point at the root), result1 );
  zAssert( zNLE1 (returned value for an equation without root), result2 );
  zAssert( zNLE1 (returned value for excess of iterations), result3 );
}

/* sin(x) = 0 on [1,9], where the function has the same sign at the end points, roots: pi, 2pi */
double sin_f(double x, void *util)
{
  count_up( util );
  return sin( x );
}

/* (x-1)(x-2.5) = 0 on [0,3], where the function has the same sign at the end points, roots: 1, 2.5 */
double quad_f(double x, void *util)
{
  count_up( util );
  return ( x - 1 )*( x - 2.5 );
}

/* x^(1/3) - 3^(1/3) = 0, which is not a number for x < 0, root: 3 */
double cbrt_f(double x, void *util)
{
  count_up( util );
  return pow( x, 1.0/3 ) - pow( 3, 1.0/3 );
}

/* x - 0.75 = 0, which is not a number on (0.4,0.6), root: 0.75 */
double hole_f(double x, void *util)
{
  count_up( util );
  return x > 0.4 && x < 0.6 ? NAN : x - 0.75;
}

/* a root has to be found by the bracketing methods even if the function has the same sign at the end points of the given range. */
void assert_nle1_search(void)
{
  struct{
    double (* f)(double,void*);
    double x1, x2, root1, root2;
  } prob[] = {
    { sin_f,  1.0, 9.0, zPI, 2*zPI },
    { quad_f, 0.0, 3.0, 1.0, 2.5 },
    { NULL, 0, 0, 0, 0 },
  }, *pp;
  method_t *meth;
  double x;
  int ret;
  bool result = true;

  for( pp=prob; pp->f; pp++ )
    for( meth=method; meth->label; meth++ ){
      if( !meth->bracketing ) continue;
      ret = meth->solve( pp->f, &x, pp->x1, pp->x2, NULL, TOL_STEP, TOL_ERR, 0 );
      if( ret < 0 || !( zIsTol( x - pp->root1, 1.0e-8 ) || zIsTol( x - pp->root2, 1.0e-8 ) ) ){
        eprintf( "%s: returned %d, x = %.16g\n", meth->label, ret, x );
        result = false;
      }
    }
  zAssert( zNLE1 (search of a bracketing range), result );
}

/* (x-1)^2 = 0, double root: 1 */
double double_f(double x, void *util)
{
  count_up( util );
  return zSqr( x - 1 );
}

/* the secant method has to find a root that is not bracketed by the initial points, including a double root. */
void assert_nle1_secant(void)
{
  double x;
  int ret;
  bool result = true;

  /* a root outside of the initial points */
  if( ( ret = zNLE1Secant( wallis_f, &x, 2.5, 2.6, NULL, TOL_STEP, TOL_ERR, 0 ) ) < 0 || !zIsTol( x - problem[1].root, 1.0e-8 ) ){
    eprintf( "Wallis: returned %d, x = %.16g\n", ret, x );
    result = false;
  }
  if( ( ret = zNLE1Secant( omega_f, &x, 0.0, 0.1, NULL, TOL_STEP, TOL_ERR, 0 ) ) < 0 || !zIsTol( x - problem[0].root, 1.0e-8 ) ){
    eprintf( "omega: returned %d, x = %.16g\n", ret, x );
    result = false;
  }
  /* double root, which can be found only with the accuracy of the square root of TOL_ERR */
  if( ( ret = zNLE1Secant( double_f, &x, 0.0, 0.5, NULL, TOL_STEP, TOL_ERR, 0 ) ) < 0 || !zIsTol( x - 1, 2*sqrt(TOL_ERR) ) ){
    eprintf( "double root: returned %d, x = %.16g\n", ret, x );
    result = false;
  }
  zAssert( zNLE1Secant (root not bracketed), result );
}

/* the step has to be shortened to avoid points where the function is not a number. */
void assert_nle1_nan(void)
{
  method_t *meth;
  double x;
  int ret;
  bool result = true;

  /* the secant method goes to x < 0 without shortening the step */
  if( ( ret = zNLE1Secant( cbrt_f, &x, 1.0, 100.0, NULL, TOL_STEP, TOL_ERR, 0 ) ) < 0 || !zIsTol( x - 3, 1.0e-8 ) ){
    eprintf( "secant: returned %d, x = %.16g\n", ret, x );
    result = false;
  }
  /* the bisection method evaluates the function at 0.5 */
  for( meth=method; meth->label; meth++ ){
    if( !meth->robust ) continue;
    if( ( ret = meth->solve( hole_f, &x, 0.0, 1.0, NULL, TOL_STEP, TOL_ERR, 0 ) ) < 0 || !zIsTol( x - 0.75, 1.0e-8 ) ){
      eprintf( "%s: returned %d, x = %.16g\n", meth->label, ret, x );
      result = false;
    }
  }
  zAssert( zNLE1 (avoidance of not-a-number), result );
}

/* the solvers must not reinitialize the random number generator. */
void assert_nle1_rand(void)
{
  method_t *meth;
  double x, r1, r2;

  zRandInitSeed( 1 );
  zRandF( 0, 1 );
  r1 = zRandF( 0, 1 );
  zRandInitSeed( 1 );
  zRandF( 0, 1 );
  zEchoOff();
  for( meth=method; meth->label; meth++ ){
    meth->solve( omega_f, &x, -10.0, 10.0, NULL, TOL_STEP, TOL_ERR, 0 );
    meth->solve( noroot_f, &x, -1.0, 1.0, NULL, TOL_STEP, TOL_ERR, 0 );
  }
  zEchoOn();
  r2 = zRandF( 0, 1 );
  zAssert( zNLE1 (random number generator not reinitialized), r1 == r2 );
}

int main(void)
{
  assert_nle1_solve();
  assert_nle1_return();
  assert_nle1_search();
  assert_nle1_secant();
  assert_nle1_nan();
  assert_nle1_rand();
  return 0;
}
