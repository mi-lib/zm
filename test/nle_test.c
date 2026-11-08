#include <zm/zm_nle.h>

/* test problems of simultaneous nonlinear equations from
 * J. J. More, B. S. Garbow and K. E. Hillstrom, Testing Unconstrained Optimization Software,
 * ACM Transactions on Mathematical Software, 7(1), pp. 17-41, 1981.
 * the number of variables n and that of equations m are given by the sizes of vectors and matrices.
 */

/* extended Rosenbrock function (m = n, n even), solution: x_i = 1 */
zVec rosenbrock_f(const zVec var, zVec f, void *dummy)
{
  int i;
  double *x = zVecBufNC(var);

  for( i=0; i<zVecSizeNC(var); i+=2 ){
    zVecSetElemNC( f, i,   10*( x[i+1] - x[i]*x[i] ) );
    zVecSetElemNC( f, i+1, 1 - x[i] );
  }
  return f;
}

zMat rosenbrock_jac(const zVec var, zMat jac, void *dummy)
{
  int i;
  double *x = zVecBufNC(var);

  zMatZero( jac );
  for( i=0; i<zVecSizeNC(var); i+=2 ){
    zMatSetElemNC( jac, i,   i,  -20*x[i] );
    zMatSetElemNC( jac, i,   i+1, 10 );
    zMatSetElemNC( jac, i+1, i,  -1 );
  }
  return jac;
}

void rosenbrock_init(zVec var)
{
  int i;

  for( i=0; i<zVecSizeNC(var); i+=2 ){
    zVecSetElemNC( var, i,  -1.2 );
    zVecSetElemNC( var, i+1, 1.0 );
  }
}

/* extended Powell singular function (m = n, n multiple of 4), solution: x_i = 0, where the Jacobian matrix is singular */
zVec powell_f(const zVec var, zVec f, void *dummy)
{
  int i;
  double *x = zVecBufNC(var);

  for( i=0; i<zVecSizeNC(var); i+=4 ){
    zVecSetElemNC( f, i,   x[i] + 10*x[i+1] );
    zVecSetElemNC( f, i+1, sqrt(5)*( x[i+2] - x[i+3] ) );
    zVecSetElemNC( f, i+2, zSqr( x[i+1] - 2*x[i+2] ) );
    zVecSetElemNC( f, i+3, sqrt(10)*zSqr( x[i] - x[i+3] ) );
  }
  return f;
}

zMat powell_jac(const zVec var, zMat jac, void *dummy)
{
  int i;
  double *x = zVecBufNC(var);

  zMatZero( jac );
  for( i=0; i<zVecSizeNC(var); i+=4 ){
    zMatSetElemNC( jac, i,   i,    1 );
    zMatSetElemNC( jac, i,   i+1, 10 );
    zMatSetElemNC( jac, i+1, i+2, sqrt(5) );
    zMatSetElemNC( jac, i+1, i+3,-sqrt(5) );
    zMatSetElemNC( jac, i+2, i+1, 2*( x[i+1] - 2*x[i+2] ) );
    zMatSetElemNC( jac, i+2, i+2,-4*( x[i+1] - 2*x[i+2] ) );
    zMatSetElemNC( jac, i+3, i,   2*sqrt(10)*( x[i] - x[i+3] ) );
    zMatSetElemNC( jac, i+3, i+3,-2*sqrt(10)*( x[i] - x[i+3] ) );
  }
  return jac;
}

void powell_init(zVec var)
{
  int i;

  for( i=0; i<zVecSizeNC(var); i+=4 ){
    zVecSetElemNC( var, i,    3 );
    zVecSetElemNC( var, i+1, -1 );
    zVecSetElemNC( var, i+2,  0 );
    zVecSetElemNC( var, i+3,  1 );
  }
}

/* trigonometric function (m = n) */
zVec trigonometric_f(const zVec var, zVec f, void *dummy)
{
  int i, n = zVecSizeNC(var);
  double *x = zVecBufNC(var), s = 0;

  for( i=0; i<n; i++ ) s += cos( x[i] );
  for( i=0; i<n; i++ )
    zVecSetElemNC( f, i, n - s + (i+1)*( 1 - cos( x[i] ) ) - sin( x[i] ) );
  return f;
}

zMat trigonometric_jac(const zVec var, zMat jac, void *dummy)
{
  int i, j, n = zVecSizeNC(var);
  double *x = zVecBufNC(var);

  for( i=0; i<n; i++ )
    for( j=0; j<n; j++ )
      zMatSetElemNC( jac, i, j, sin( x[j] ) + ( i == j ? (i+1)*sin( x[i] ) - cos( x[i] ) : 0 ) );
  return jac;
}

void trigonometric_init(zVec var)
{
  zVecSetAll( var, 1.0/zVecSizeNC(var) );
}

/* variably dimensioned function (m = n+2), solution: x_i = 1 */
zVec variably_dimensioned_f(const zVec var, zVec f, void *dummy)
{
  int i, n = zVecSizeNC(var);
  double *x = zVecBufNC(var), s = 0;

  for( i=0; i<n; i++ ){
    zVecSetElemNC( f, i, x[i] - 1 );
    s += (i+1)*( x[i] - 1 );
  }
  zVecSetElemNC( f, n,   s );
  zVecSetElemNC( f, n+1, s*s );
  return f;
}

zMat variably_dimensioned_jac(const zVec var, zMat jac, void *dummy)
{
  int i, n = zVecSizeNC(var);
  double *x = zVecBufNC(var), s = 0;

  zMatZero( jac );
  for( i=0; i<n; i++ ) s += (i+1)*( x[i] - 1 );
  for( i=0; i<n; i++ ){
    zMatSetElemNC( jac, i,   i, 1 );
    zMatSetElemNC( jac, n,   i, i+1 );
    zMatSetElemNC( jac, n+1, i, 2*s*(i+1) );
  }
  return jac;
}

void variably_dimensioned_init(zVec var)
{
  int i, n = zVecSizeNC(var);

  for( i=0; i<n; i++ ) zVecSetElemNC( var, i, 1 - (double)(i+1)/n );
}

/* Brown almost-linear function (m = n) */
zVec brown_almost_linear_f(const zVec var, zVec f, void *dummy)
{
  int i, n = zVecSizeNC(var);
  double *x = zVecBufNC(var), s = 0, p = 1;

  for( i=0; i<n; i++ ){
    s += x[i];
    p *= x[i];
  }
  for( i=0; i<n-1; i++ ) zVecSetElemNC( f, i, x[i] + s - (n+1) );
  zVecSetElemNC( f, n-1, p - 1 );
  return f;
}

zMat brown_almost_linear_jac(const zVec var, zMat jac, void *dummy)
{
  int i, j, n = zVecSizeNC(var);
  double *x = zVecBufNC(var), p;

  for( i=0; i<n-1; i++ )
    for( j=0; j<n; j++ )
      zMatSetElemNC( jac, i, j, i == j ? 2 : 1 );
  for( j=0; j<n; j++ ){
    for( p=1, i=0; i<n; i++ )
      if( i != j ) p *= x[i];
    zMatSetElemNC( jac, n-1, j, p );
  }
  return jac;
}

void brown_almost_linear_init(zVec var)
{
  zVecSetAll( var, 0.5 );
}

/* Broyden tridiagonal function (m = n) */
zVec broyden_tridiagonal_f(const zVec var, zVec f, void *dummy)
{
  int i, n = zVecSizeNC(var);
  double *x = zVecBufNC(var);

  for( i=0; i<n; i++ )
    zVecSetElemNC( f, i, ( 3 - 2*x[i] )*x[i] - ( i > 0 ? x[i-1] : 0 ) - 2*( i < n-1 ? x[i+1] : 0 ) + 1 );
  return f;
}

zMat broyden_tridiagonal_jac(const zVec var, zMat jac, void *dummy)
{
  int i, n = zVecSizeNC(var);
  double *x = zVecBufNC(var);

  zMatZero( jac );
  for( i=0; i<n; i++ ){
    zMatSetElemNC( jac, i, i, 3 - 4*x[i] );
    if( i > 0 )   zMatSetElemNC( jac, i, i-1, -1 );
    if( i < n-1 ) zMatSetElemNC( jac, i, i+1, -2 );
  }
  return jac;
}

void broyden_tridiagonal_init(zVec var)
{
  zVecSetAll( var, -1.0 );
}

/* discrete boundary value function (m = n) */
zVec discrete_boundary_value_f(const zVec var, zVec f, void *dummy)
{
  int i, n = zVecSizeNC(var);
  double *x = zVecBufNC(var), h = 1.0/(n+1);

  for( i=0; i<n; i++ )
    zVecSetElemNC( f, i, 2*x[i] - ( i > 0 ? x[i-1] : 0 ) - ( i < n-1 ? x[i+1] : 0 ) + 0.5*h*h*pow( x[i] + (i+1)*h + 1, 3 ) );
  return f;
}

zMat discrete_boundary_value_jac(const zVec var, zMat jac, void *dummy)
{
  int i, n = zVecSizeNC(var);
  double *x = zVecBufNC(var), h = 1.0/(n+1);

  zMatZero( jac );
  for( i=0; i<n; i++ ){
    zMatSetElemNC( jac, i, i, 2 + 1.5*h*h*zSqr( x[i] + (i+1)*h + 1 ) );
    if( i > 0 )   zMatSetElemNC( jac, i, i-1, -1 );
    if( i < n-1 ) zMatSetElemNC( jac, i, i+1, -1 );
  }
  return jac;
}

void discrete_boundary_value_init(zVec var)
{
  int i, n = zVecSizeNC(var);
  double t, h = 1.0/(n+1);

  for( i=0; i<n; i++ ){
    t = (i+1)*h;
    zVecSetElemNC( var, i, t*( t - 1 ) );
  }
}

/* penalty function I (m = n+1), which has no exact solution */
#define PENALTY_A 1.0e-5
zVec penalty1_f(const zVec var, zVec f, void *dummy)
{
  int i, n = zVecSizeNC(var);
  double *x = zVecBufNC(var), s = 0;

  for( i=0; i<n; i++ ){
    zVecSetElemNC( f, i, sqrt(PENALTY_A)*( x[i] - 1 ) );
    s += x[i]*x[i];
  }
  zVecSetElemNC( f, n, s - 0.25 );
  return f;
}

zMat penalty1_jac(const zVec var, zMat jac, void *dummy)
{
  int i, n = zVecSizeNC(var);
  double *x = zVecBufNC(var);

  zMatZero( jac );
  for( i=0; i<n; i++ ){
    zMatSetElemNC( jac, i, i, sqrt(PENALTY_A) );
    zMatSetElemNC( jac, n, i, 2*x[i] );
  }
  return jac;
}

void penalty1_init(zVec var)
{
  int i;

  for( i=0; i<zVecSizeNC(var); i++ ) zVecSetElemNC( var, i, i+1 );
}

/* underdetermined systems (m < n) */
/* a point on the unit sphere (m = 1, n = 3) */
zVec sphere_f(const zVec var, zVec f, void *dummy)
{
  zVecSetElemNC( f, 0, zVecSqrNorm(var) - 1 );
  return f;
}

zMat sphere_jac(const zVec var, zMat jac, void *dummy)
{
  int i;

  for( i=0; i<zVecSizeNC(var); i++ )
    zMatSetElemNC( jac, 0, i, 2*zVecElemNC(var,i) );
  return jac;
}

void sphere_init(zVec var)
{
  zVecSetElemList( var, 2.0, 1.0, 0.5 );
}

/* a point on the intersection of two unit spheres centered at the origin and (1,0,0) (m = 2, n = 3) */
zVec two_spheres_f(const zVec var, zVec f, void *dummy)
{
  double *x = zVecBufNC(var);

  zVecSetElemNC( f, 0, x[0]*x[0] + x[1]*x[1] + x[2]*x[2] - 1 );
  zVecSetElemNC( f, 1, zSqr(x[0]-1) + x[1]*x[1] + x[2]*x[2] - 1 );
  return f;
}

zMat two_spheres_jac(const zVec var, zMat jac, void *dummy)
{
  double *x = zVecBufNC(var);

  zMatSetElemNC( jac, 0, 0, 2*x[0] );
  zMatSetElemNC( jac, 0, 1, 2*x[1] );
  zMatSetElemNC( jac, 0, 2, 2*x[2] );
  zMatSetElemNC( jac, 1, 0, 2*(x[0]-1) );
  zMatSetElemNC( jac, 1, 1, 2*x[1] );
  zMatSetElemNC( jac, 1, 2, 2*x[2] );
  return jac;
}

void two_spheres_init(zVec var)
{
  zVecSetElemList( var, 2.0, 1.0, 0.5 );
}

typedef struct{
  const char *name;
  int nv, ne;
  zVec (* f)(const zVec,zVec,void*);
  zMat (* jac)(const zVec,zMat,void*);
  void (* init)(zVec);
  bool consistent; /* true if the system has an exact solution */
} problem_t;

problem_t problem[] = {
  { "extended Rosenbrock",     10, 10, rosenbrock_f,              rosenbrock_jac,              rosenbrock_init,              true },
  { "extended Powell singular",12, 12, powell_f,                  powell_jac,                  powell_init,                  true },
  { "trigonometric",           10, 10, trigonometric_f,           trigonometric_jac,           trigonometric_init,           true },
  { "variably dimensioned",    10, 12, variably_dimensioned_f,    variably_dimensioned_jac,    variably_dimensioned_init,    true },
  { "Brown almost-linear",     10, 10, brown_almost_linear_f,     brown_almost_linear_jac,     brown_almost_linear_init,     true },
  { "Broyden tridiagonal",     10, 10, broyden_tridiagonal_f,     broyden_tridiagonal_jac,     broyden_tridiagonal_init,     true },
  { "discrete boundary value", 10, 10, discrete_boundary_value_f, discrete_boundary_value_jac, discrete_boundary_value_init, true },
  { "penalty I",               10, 11, penalty1_f,                penalty1_jac,                penalty1_init,                false },
  { "sphere",                   3,  1, sphere_f,                  sphere_jac,                  sphere_init,                  true },
  { "two spheres",              3,  2, two_spheres_f,             two_spheres_jac,             two_spheres_init,             true },
  { NULL, 0, 0, NULL, NULL, NULL, false },
};

typedef struct{
  const char *label;
  zNLE *(* assign)(zNLE*);
  bool robust; /* true if the method is expected to solve all consistent problems */
} method_t;

zNLE *assign_nr(zNLE *nle){ return zNLEAssignNR( nle ); }
zNLE *assign_broyden(zNLE *nle){ return zNLEAssignBroyden( nle ); }
zNLE *assign_lm_grad(zNLE *nle){ return zNLEAssignLM( nle, "grad", NULL, 0 ); }
zNLE *assign_lm_dec(zNLE *nle){ return zNLEAssignLM( nle, "dec", NULL, 0 ); }
zNLE *assign_vm(zNLE *nle){ return zNLEAssignVM( nle, NULL, NULL, 0 ); }
zNLE *assign_sd(zNLE *nle){ return zNLEAssignSD( nle, NULL, 0 ); }
zNLE *assign_cg(zNLE *nle){ return zNLEAssignCG( nle, NULL, 0 ); }

method_t method[] = {
  { "NR",       assign_nr,      false },
  { "Broyden",  assign_broyden, false },
  { "LM (grad)",assign_lm_grad, true },
  { "LM (dec)", assign_lm_dec,  true },
  { "VM",       assign_vm,      false },
  { "SD",       assign_sd,      false },
  { "CG",       assign_cg,      false },
  { NULL, NULL, false },
};

#define TOL_STEP zTOL
#define TOL_ERR  zTOL
#define ITER     1000

/* solve a problem by a method, and return the returned value of zNLESolve. the residual vector at the result is stored in res. */
int solve(problem_t *prob, method_t *meth, bool use_jac, zVec var, zVec res, int iter, double *err)
{
  zNLE nle;
  int ret;

  zNLECreate( &nle, prob->nv, prob->ne, prob->f, use_jac ? prob->jac : NULL );
  meth->assign( &nle );
  zEchoOff(); /* to suppress warnings of unconverged iterations */
  ret = zNLESolve( &nle, var, NULL, TOL_STEP, TOL_ERR, iter, err );
  zEchoOn();
  zNLEDestroy( &nle );
  prob->f( var, res, NULL );
  return ret;
}

/* the success has to be reported only if the residual is small, and the residual through the last argument has to be that at the result. */
bool check_result(problem_t *prob, method_t *meth, const char *jaclabel, int ret, zVec res, double err)
{
  bool result = true;

  if( ret >= 0 && !zVecIsTol( res, TOL_ERR ) ){
    eprintf( "%s, %s (%s): succeeded with |f| = %g\n", prob->name, meth->label, jaclabel, zVecElemAbsMax(res,NULL) );
    result = false;
  }
  if( !zIsTol( err - 0.5*zVecSqrNorm(res), zTOL*( 1 + err ) ) ){
    eprintf( "%s, %s (%s): residual returned %g, actual %g\n", prob->name, meth->label, jaclabel, err, 0.5*zVecSqrNorm(res) );
    result = false;
  }
  return result;
}

/* every method has to succeed only if the residual becomes small, robust methods have to solve all consistent
 * problems, and every method has to fail for an inconsistent problem. */
void assert_nle_solve(bool use_jac)
{
  problem_t *prob;
  method_t *meth;
  zVec var, res;
  double err;
  int ret;
  bool result1 = true, result2 = true, result3 = true;
  const char *jaclabel = use_jac ? "analytic Jacobian" : "numerical Jacobian";

  for( prob=problem; prob->name; prob++ ){
    var = zVecAlloc( prob->nv );
    res = zVecAlloc( prob->ne );
    for( meth=method; meth->label; meth++ ){
      if( meth->assign == assign_broyden && prob->nv != prob->ne ) continue; /* see assert_nle_broyden_nonsquare */
      prob->init( var );
      ret = solve( prob, meth, use_jac, var, res, ITER, &err );
      if( !check_result( prob, meth, jaclabel, ret, res, err ) ) result1 = false;
      if( prob->consistent && meth->robust && ret < 0 ){
        eprintf( "%s, %s (%s): failed with |f| = %g\n", prob->name, meth->label, jaclabel, zVecElemAbsMax(res,NULL) );
        result2 = false;
      }
      if( !prob->consistent && ret >= 0 ){
        eprintf( "%s, %s (%s): succeeded for an inconsistent problem\n", prob->name, meth->label, jaclabel );
        result3 = false;
      }
    }
    zVecFreeAtOnce( 2, var, res );
  }
  if( use_jac ){
    zAssert( zNLESolve (success only with small residual, analytic Jacobian), result1 );
    zAssert( zNLESolve (Levenberg-Marquardt for consistent problems, analytic Jacobian), result2 );
    zAssert( zNLESolve (failure for inconsistent problem, analytic Jacobian), result3 );
  } else{
    zAssert( zNLESolve (success only with small residual, numerical Jacobian), result1 );
    zAssert( zNLESolve (Levenberg-Marquardt for consistent problems, numerical Jacobian), result2 );
    zAssert( zNLESolve (failure for inconsistent problem, numerical Jacobian), result3 );
  }
}

/* underdetermined systems have to be solved by Newton-Raphson's method and Levenberg-Marquardt's method. */
void assert_nle_underdetermined(void)
{
  problem_t *prob;
  method_t *meth;
  zVec var, res;
  double err;
  int ret;
  bool result = true;

  for( prob=problem; prob->name; prob++ ){
    if( prob->ne >= prob->nv ) continue;
    var = zVecAlloc( prob->nv );
    res = zVecAlloc( prob->ne );
    for( meth=method; meth->label; meth++ ){
      if( meth->assign != assign_nr && !meth->robust ) continue;
      prob->init( var );
      ret = solve( prob, meth, false, var, res, ITER, &err );
      if( ret < 0 || !zVecIsTol( res, TOL_ERR ) ){
        eprintf( "%s, %s: returned %d with |f| = %g\n", prob->name, meth->label, ret, zVecElemAbsMax(res,NULL) );
        result = false;
      }
    }
    zVecFreeAtOnce( 2, var, res );
  }
  zAssert( zNLESolve (underdetermined system), result );
}

/* Broyden's method does not accept non-square systems in the current implementation. */
void assert_nle_broyden_nonsquare(void)
{
  problem_t *prob;
  method_t meth = { "Broyden", assign_broyden, false };
  zVec var, res;
  double err;
  bool result = true;

  for( prob=problem; prob->name; prob++ ){
    if( prob->nv == prob->ne ) continue;
    var = zVecAlloc( prob->nv );
    res = zVecAlloc( prob->ne );
    prob->init( var );
    if( solve( prob, &meth, false, var, res, ITER, &err ) != -1 ){
      eprintf( "%s: Broyden's method did not fail\n", prob->name );
      result = false;
    }
    zVecFreeAtOnce( 2, var, res );
  }
  zAssert( zNLESolve (Broyden method for non-square system), result );
}

/* the analytic Jacobian matrix has to coincide with the numerical one. */
void assert_nle_jacobian(void)
{
  problem_t *prob;
  zNLE nle;
  zVec var;
  zMat jn, ja;
  bool result = true;
  const double tol = 1.0e-6;

  for( prob=problem; prob->name; prob++ ){
    var = zVecAlloc( prob->nv );
    jn = zMatAlloc( prob->ne, prob->nv );
    ja = zMatAlloc( prob->ne, prob->nv );
    zNLECreate( &nle, prob->nv, prob->ne, prob->f, NULL );
    prob->init( var );
    nle._jac( &nle, var, jn, NULL );
    prob->jac( var, ja, NULL );
    if( !zMatIsTol( zMatSubDRC( jn, ja ), tol * ( 1 + zMatElemAbsMax( ja, NULL, NULL ) ) ) ){
      eprintf( "%s: |Jn-Ja| = %g, |Ja| = %g\n", prob->name, zMatElemAbsMax( jn, NULL, NULL ), zMatElemAbsMax( ja, NULL, NULL ) );
      result = false;
    }
    zNLEDestroy( &nle );
    zMatFreeAtOnce( 2, jn, ja );
    zVecFree( var );
  }
  zAssert( zNLE (numerical Jacobian), result );
}

/* the returned value has to be 0 if the initial guess is the solution, and -1 if the number of iterations exceeds the limit. */
void assert_nle_return(void)
{
  problem_t *prob = &problem[0]; /* extended Rosenbrock function */
  method_t *meth;
  zVec var, res;
  double err;
  int ret;
  bool result1 = true, result2 = true, result3 = true;

  var = zVecAlloc( prob->nv );
  res = zVecAlloc( prob->ne );
  for( meth=method; meth->label; meth++ ){
    zVecSetAll( var, 1.0 );
    if( ( ret = solve( prob, meth, false, var, res, ITER, &err ) ) != 0 ){
      eprintf( "%s: returned %d for the initial guess at the solution\n", meth->label, ret );
      result1 = false;
    }
    prob->init( var );
    if( ( ret = solve( prob, meth, false, var, res, 2, &err ) ) != -1 ){
      eprintf( "%s: returned %d within 2 iterations\n", meth->label, ret );
      result2 = false;
    }
    prob->init( var );
    if( ( ret = solve( prob, meth, false, var, res, ITER, &err ) ) == 0 || ret > ITER ){
      eprintf( "%s: returned %d\n", meth->label, ret );
      result3 = false;
    }
  }
  zVecFreeAtOnce( 2, var, res );
  zAssert( zNLESolve (returned value for the initial guess at the solution), result1 );
  zAssert( zNLESolve (returned value for excess of iterations), result2 );
  zAssert( zNLESolve (returned number of iterations), result3 );
}

int main(void)
{
  assert_nle_jacobian();
  assert_nle_solve( false );
  assert_nle_solve( true );
  assert_nle_underdetermined();
  assert_nle_broyden_nonsquare();
  assert_nle_return();
  return 0;
}
