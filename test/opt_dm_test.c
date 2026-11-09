#include <zm/zm_opt.h>

#define DIM 4

/* quadratic function: f(x) = k/2 sum_i (i+1)(x_i-i)^2, minimum at x_i = i, where k is given via util (1 if null) */
#define quad_scale(u) ( (u) ? *(double *)(u) : 1.0 )
double quad_eval(const zVec var, void *util)
{
  int i;
  double val = 0;

  for( i=0; i<zVecSizeNC(var); i++ )
    val += (i+1)*zSqr( zVecElemNC(var,i) - i );
  return 0.5*quad_scale(util)*val;
}

zVec quad_grad(const zVec var, zVec g, void *util)
{
  int i;

  for( i=0; i<zVecSizeNC(var); i++ )
    zVecSetElemNC( g, i, quad_scale(util)*(i+1)*( zVecElemNC(var,i) - i ) );
  return g;
}

/* Rosenbrock function: f(x) = sum_i 100(x_2i+1-x_2i^2)^2+(1-x_2i)^2, minimum at x_i = 1 */
double rosenbrock_eval(const zVec var, void *dummy)
{
  int i;
  double x, y, val = 0;

  for( i=0; i<zVecSizeNC(var); i+=2 ){
    x = zVecElemNC(var,i);
    y = zVecElemNC(var,i+1);
    val += 100*zSqr( y - x*x ) + zSqr( 1 - x );
  }
  return val;
}

zVec rosenbrock_grad(const zVec var, zVec g, void *dummy)
{
  int i;
  double x, y;

  for( i=0; i<zVecSizeNC(var); i+=2 ){
    x = zVecElemNC(var,i);
    y = zVecElemNC(var,i+1);
    zVecSetElemNC( g, i,  -400*x*(y-x*x) + 2*(x-1) );
    zVecSetElemNC( g, i+1, 200*(y-x*x) );
  }
  return g;
}

typedef struct{
  const char *label;
  zOptDM *(* assign)(zOptDM*,const char*,const char*);
  const char *stepmethod;
  const char *updatemethod;
} method_t;

zOptDM *assign_sd(zOptDM *opt, const char *stepmethod, const char *updatemethod){ return zOptDMAssignSD( opt, stepmethod, 0 ); }
zOptDM *assign_lm(zOptDM *opt, const char *stepmethod, const char *updatemethod){ return zOptDMAssignLM( opt, NULL, stepmethod, 0 ); }
zOptDM *assign_vm(zOptDM *opt, const char *stepmethod, const char *updatemethod){ return zOptDMAssignVM( opt, updatemethod, stepmethod, 0 ); }
zOptDM *assign_cg(zOptDM *opt, const char *stepmethod, const char *updatemethod){ return zOptDMAssignCG( opt, stepmethod, 0 ); }

method_t method[] = {
  { "SD (MT)",         assign_sd, "MT",    NULL },
  { "SD (DS)",         assign_sd, "DS",    NULL },
  { "SD (Brent)",      assign_sd, "Brent", NULL },
  { "LM (MT)",         assign_lm, "MT",    NULL },
  { "LM (DS)",         assign_lm, "DS",    NULL },
  { "VM (MT, BFGS)",   assign_vm, "MT",    "BFGS" },
  { "VM (DS, BFGS)",   assign_vm, "DS",    "BFGS" },
  { "VM (Brent, BFGS)",assign_vm, "Brent", "BFGS" },
  { "CG (Brent)",      assign_cg, NULL,    NULL },
  { "CG (MT)",         assign_cg, "MT",    NULL },
  { NULL, NULL, NULL, NULL },
};

/* the value returned through the last argument of zOptDMSolve has to be the evaluation at the returned solution. */
void assert_opt_dm_eval(void)
{
  zOptDM opt;
  zVec var;
  method_t *m;
  double eval;
  bool result = true;

  var = zVecAlloc( DIM );
  zOptDMCreate( &opt, DIM, quad_eval, quad_grad, NULL );
  for( m=method; m->label; m++ ){
    m->assign( &opt, m->stepmethod, m->updatemethod );
    zVecSetAll( var, 1.0 );
    if( zOptDMSolve( &opt, var, NULL, zTOL, 0, &eval ) < 0 ){
      eprintf( "%s: not converged\n", m->label );
      result = false;
    }
    if( eval != quad_eval( var, NULL ) ){
      eprintf( "%s: returned evaluation = %g, actual evaluation = %g\n", m->label, eval, quad_eval( var, NULL ) );
      result = false;
    }
  }
  zOptDMDestroy( &opt );
  zVecFree( var );
  zAssert( zOptDMSolve (returned evaluation), result );
}

/* Dennis-Schnabel's backtracking has to terminate and converge regardless of the scale of the function. */
void assert_opt_dm_ds(void)
{
  zOptDM opt;
  zVec var, ans;
  double k[] = { 1.0, 1.0e2, 1.0e4, 1.0e6, 0 }, *kp;
  bool result = true;
  int i;
  const double tol = 1.0e-6;

  var = zVecAlloc( DIM );
  ans = zVecAlloc( DIM );
  /* scaled quadratic functions */
  for( i=0; i<DIM; i++ ) zVecSetElemNC( ans, i, i );
  zOptDMCreate( &opt, DIM, quad_eval, quad_grad, NULL );
  for( kp=k; *kp>0; kp++ ){
    for( i=0; i<2; i++ ){
      if( i == 0 )
        zOptDMAssignSD( &opt, "DS", 0 );
      else
        zOptDMAssignVM( &opt, "BFGS", "DS", 0 );
      zVecSetAll( var, 1.0 );
      if( zOptDMSolve( &opt, var, kp, zTOL, 0, NULL ) < 0 || zVecDist( var, ans ) > tol ){
        eprintf( "%s (DS), k=%g: |x-x*| = %g\n", i == 0 ? "SD" : "VM", *kp, zVecDist( var, ans ) );
        result = false;
      }
    }
  }
  zOptDMDestroy( &opt );
  /* Rosenbrock function, where the original implementation fell into an infinite loop */
  zVecSetAll( ans, 1.0 );
  zOptDMCreate( &opt, DIM, rosenbrock_eval, rosenbrock_grad, NULL );
  zOptDMAssignVM( &opt, "BFGS", "DS", 0 );
  zVecSetAll( var, -1.2 );
  if( zOptDMSolve( &opt, var, NULL, zTOL, 0, NULL ) < 0 || zVecDist( var, ans ) > tol ){
    eprintf( "VM (DS), Rosenbrock: |x-x*| = %g\n", zVecDist( var, ans ) );
    result = false;
  }
  zOptDMDestroy( &opt );
  zVecFreeAtOnce( 2, var, ans );
  zAssert( zOptDMSolve (Dennis-Schnabel), result );
}

/* the inverse Hessian updated by DFP and BFGS formulae from the identity has to coincide with the closed forms, and the update has to be skipped if p^T q <= 0. */
void assert_opt_dm_update(void)
{
  zOptDM opt;
  zMat h, hp;
  double pq, qq;
  int i, j;
  bool result1 = true, result2 = true, result3 = true;

  h = zMatAllocSqr( DIM );
  hp = zMatAllocSqr( DIM );
  zOptDMCreate( &opt, DIM, quad_eval, quad_grad, NULL );
  for( i=0; i<DIM; i++ ){
    zVecSetElemNC( opt._p, i, zRandF(-1,1) );
    zVecSetElemNC( opt._q, i, zRandF(-1,1) );
  }
  if( ( pq = zVecInnerProd( opt._p, opt._q ) ) < 0 ){
    zVecRevDRC( opt._q );
    pq = -pq;
  }
  qq = zVecSqrNorm( opt._q );
  /* DFP: H = I + p p^T / p^T q - q q^T / q^T q */
  zOptDMAssignVM( &opt, "DFP", NULL, 0 );
  zMatIdent( opt._h );
  opt._update_h( &opt, 0 );
  for( i=0; i<DIM; i++ )
    for( j=0; j<DIM; j++ )
      zMatElemNC(h,i,j) = ( i == j ? 1 : 0 ) + zVecElemNC(opt._p,i)*zVecElemNC(opt._p,j)/pq - zVecElemNC(opt._q,i)*zVecElemNC(opt._q,j)/qq;
  if( !zMatEqual( opt._h, h, zTOL ) ){
    eprintf( "DFP update differs from the closed form\n" );
    result1 = false;
  }
  /* BFGS: H = I + ( ( 1 + q^T q / p^T q ) p p^T - p q^T - q p^T ) / p^T q */
  zOptDMAssignVM( &opt, "BFGS", NULL, 0 );
  zMatIdent( opt._h );
  opt._update_h( &opt, 0 );
  for( i=0; i<DIM; i++ )
    for( j=0; j<DIM; j++ )
      zMatElemNC(h,i,j) = ( i == j ? 1 : 0 ) + ( ( 1 + qq/pq )*zVecElemNC(opt._p,i)*zVecElemNC(opt._p,j) - zVecElemNC(opt._p,i)*zVecElemNC(opt._q,j) - zVecElemNC(opt._q,i)*zVecElemNC(opt._p,j) ) / pq;
  if( !zMatEqual( opt._h, h, zTOL ) ){
    eprintf( "BFGS update differs from the closed form\n" );
    result2 = false;
  }
  /* skip if p^T q < 0 */
  zVecRevDRC( opt._q );
  zMatCopy( opt._h, hp );
  opt._update_h( &opt, 0 );
  if( !zMatEqual( opt._h, hp, 0 ) ){
    eprintf( "BFGS update not skipped for p^T q < 0\n" );
    result3 = false;
  }
  zOptDMAssignVM( &opt, "DFP", NULL, 0 );
  opt._update_h( &opt, 0 );
  if( !zMatEqual( opt._h, hp, 0 ) ){
    eprintf( "DFP update not skipped for p^T q < 0\n" );
    result3 = false;
  }
  zOptDMDestroy( &opt );
  zMatFreeAtOnce( 2, h, hp );
  zAssert( zOptDM (DFP update), result1 );
  zAssert( zOptDM (BFGS update), result2 );
  zAssert( zOptDM (update skipped for non-positive curvature), result3 );
}

/* saddle function: f(x) = sum_i x_2i^2 - x_2i+1^2 + x_2i+1^4, saddle at the origin, minima at x_2i = 0, x_2i+1 = +-1/sqrt(2) */
double saddle_eval(const zVec var, void *dummy)
{
  int i;
  double x, y, val = 0;

  for( i=0; i<zVecSizeNC(var); i+=2 ){
    x = zVecElemNC(var,i);
    y = zVecElemNC(var,i+1);
    val += x*x - y*y + y*y*y*y;
  }
  return val;
}

zVec saddle_grad(const zVec var, zVec g, void *dummy)
{
  int i;
  double x, y;

  for( i=0; i<zVecSizeNC(var); i+=2 ){
    x = zVecElemNC(var,i);
    y = zVecElemNC(var,i+1);
    zVecSetElemNC( g, i,   2*x );
    zVecSetElemNC( g, i+1,-2*y + 4*y*y*y );
  }
  return g;
}

zMat saddle_hess(const zVec var, zMat h, void *dummy)
{
  int i;
  double y;

  zMatZero( h );
  for( i=0; i<zVecSizeNC(var); i+=2 ){
    y = zVecElemNC(var,i+1);
    zMatSetElemNC( h, i,   i,   2 );
    zMatSetElemNC( h, i+1, i+1,-2 + 12*y*y );
  }
  return h;
}

/* the search has to restart from the steepest descent direction if the given direction is not descent.
   near the saddle point, the direction of Levenberg-Marquardt's method is not descent since the Hessian matrix is indefinite. */
void assert_opt_dm_descent(void)
{
  zOptDM opt;
  zVec var;
  const char *stepmethod[] = { "MT", "DS", "Brent", NULL }, **sp;
  bool result = true;
  const double tol = 1.0e-6;

  var = zVecAlloc( DIM );
  zOptDMCreate( &opt, DIM, saddle_eval, saddle_grad, saddle_hess );
  for( sp=stepmethod; *sp; sp++ ){
    zOptDMAssignLM( &opt, NULL, *sp, 0 );
    zVecSetElemList( var, 0.001, 0.1, 0.001, -0.1 );
    if( zOptDMSolve( &opt, var, NULL, zTOL, 0, NULL ) < 0 || !zEqual( saddle_eval( var, NULL ), -0.5, tol ) ){
      eprintf( "LM (%s), saddle: f = %g (to be -0.5)\n", *sp, saddle_eval( var, NULL ) );
      result = false;
    }
  }
  zOptDMDestroy( &opt );
  zVecFree( var );
  zAssert( zOptDMSolve (restart for non-descent direction), result );
}

/* More-Thuente's line search has to find a step length that satisfies the strong Wolfe's conditions. */
void assert_opt_dm_mt_wolfe(void)
{
  zOptDM opt;
  zVec x0, x1, g0, g1, dx;
  double a, s0;
  int i;
  bool result = true;
  const int n = 100;

  x0 = zVecAlloc( DIM );
  x1 = zVecAlloc( DIM );
  g0 = zVecAlloc( DIM );
  g1 = zVecAlloc( DIM );
  dx = zVecAlloc( DIM );
  zOptDMCreate( &opt, DIM, rosenbrock_eval, rosenbrock_grad, NULL );
  zOptDMAssignSD( &opt, "MT", 0 );
  zEchoOff(); /* to suppress warnings of unconverged iterations */
  for( i=0; i<n; i++ ){
    zVecRandUniform( x0, -2, 2 );
    zVecCopy( x0, x1 );
    zOptDMSolve( &opt, x1, NULL, zTOL, 1, NULL ); /* one step along the steepest descent direction */
    rosenbrock_grad( x0, g0, NULL );
    rosenbrock_grad( x1, g1, NULL );
    zVecSub( x1, x0, dx );
    s0 = zVecSqrNorm( g0 );
    a = -zVecInnerProd( dx, g0 ) / s0;
    if( rosenbrock_eval( x1, NULL ) > rosenbrock_eval( x0, NULL ) - opt._df * a * s0 ||
        fabs( zVecInnerProd( g1, g0 ) ) > opt._cf * s0 ){
      result = false;
      break;
    }
  }
  zEchoOn();
  zOptDMDestroy( &opt );
  zVecFreeAtOnce( 5, x0, x1, g0, g1, dx );
  zAssert( zOptDM (More-Thuente's line search: strong Wolfe's conditions), result );
}

/* More-Thuente's line search has to converge regardless of the scale of the function. */
void assert_opt_dm_mt(void)
{
  zOptDM opt;
  zVec var, ans;
  double k[] = { 1.0e-2, 1.0, 1.0e2, 1.0e4, 1.0e6, 0 }, *kp;
  bool result = true;
  int i;
  const double tol = 1.0e-6;

  var = zVecAlloc( DIM );
  ans = zVecAlloc( DIM );
  /* scaled quadratic functions */
  for( i=0; i<DIM; i++ ) zVecSetElemNC( ans, i, i );
  zOptDMCreate( &opt, DIM, quad_eval, quad_grad, NULL );
  for( kp=k; *kp>0; kp++ ){
    for( i=0; i<3; i++ ){
      if( i == 0 )
        zOptDMAssignSD( &opt, "MT", 0 );
      else if( i == 1 )
        zOptDMAssignLM( &opt, NULL, "MT", 0 );
      else
        zOptDMAssignVM( &opt, "BFGS", "MT", 0 );
      zVecSetAll( var, 1.0 );
      if( zOptDMSolve( &opt, var, kp, zTOL, 0, NULL ) < 0 || zVecDist( var, ans ) > tol ){
        eprintf( "%s (MT), k=%g: |x-x*| = %g\n", i == 0 ? "SD" : i == 1 ? "LM" : "VM", *kp, zVecDist( var, ans ) );
        result = false;
      }
    }
  }
  zOptDMDestroy( &opt );
  /* Rosenbrock function */
  zVecSetAll( ans, 1.0 );
  zOptDMCreate( &opt, DIM, rosenbrock_eval, rosenbrock_grad, NULL );
  for( i=0; i<2; i++ ){
    if( i == 0 )
      zOptDMAssignLM( &opt, NULL, "MT", 0 );
    else
      zOptDMAssignVM( &opt, "BFGS", "MT", 0 );
    zVecSetAll( var, -1.2 );
    if( zOptDMSolve( &opt, var, NULL, zTOL, 0, NULL ) < 0 || zVecDist( var, ans ) > tol ){
      eprintf( "%s (MT), Rosenbrock: |x-x*| = %g\n", i == 0 ? "LM" : "VM", zVecDist( var, ans ) );
      result = false;
    }
  }
  zOptDMDestroy( &opt );
  zVecFreeAtOnce( 2, var, ans );
  zAssert( zOptDMSolve (More-Thuente), result );
}

/* the fixed step without line search has to be the descent vector itself in the opposite direction. */
void assert_opt_dm_fixed(void)
{
  zOptDM opt;
  zVec var, ans, x0, g0;
  double k;
  int i;
  bool result1 = true, result2 = true;
  const double tol = 1.0e-6;

  var = zVecAlloc( DIM );
  ans = zVecAlloc( DIM );
  x0 = zVecAlloc( DIM );
  g0 = zVecAlloc( DIM );
  for( i=0; i<DIM; i++ ) zVecSetElemNC( ans, i, i );
  zOptDMCreate( &opt, DIM, quad_eval, quad_grad, NULL );
  /* one step of steepest descent method */
  zOptDMAssignSD( &opt, "fixed", 0 );
  zVecSetAll( x0, 1.0 );
  zVecCopy( x0, var );
  zEchoOff(); /* to suppress a warning of unconverged iteration */
  zOptDMSolve( &opt, var, NULL, zTOL, 1, NULL );
  zEchoOn();
  quad_grad( x0, g0, NULL );
  zVecSubDRC( x0, g0 );
  if( !zVecEqual( var, x0, zTOL ) ){
    eprintf( "SD (fixed): wrong step\n" );
    result1 = false;
  }
  /* one step with a given step length */
  zOptDMAssignSD( &opt, "fixed", 0.5 );
  zVecSetAll( x0, 1.0 );
  zVecCopy( x0, var );
  zEchoOff();
  zOptDMSolve( &opt, var, NULL, zTOL, 1, NULL );
  zEchoOn();
  quad_grad( x0, g0, NULL );
  zVecCatNCDRC( x0, -0.5, g0 );
  if( !zVecEqual( var, x0, zTOL ) ){
    eprintf( "SD (fixed, 0.5): wrong step\n" );
    result1 = false;
  }
  /* convergence. the steepest descent method with the unit step converges only if the curvature is less than 2. */
  for( i=0; i<4; i++ ){
    k = 1.0;
    if( i == 0 ){
      zOptDMAssignSD( &opt, "fixed", 0 );
      k = 0.25;
    } else if( i == 1 )
      zOptDMAssignLM( &opt, NULL, "fixed", 0 );
    else
      zOptDMAssignVM( &opt, i == 2 ? "BFGS" : "DFP", "fixed", 0 );
    zVecSetAll( var, 1.0 );
    if( zOptDMSolve( &opt, var, &k, zTOL, 0, NULL ) < 0 || zVecDist( var, ans ) > tol ){
      eprintf( "%s (fixed): |x-x*| = %g\n", i == 0 ? "SD" : i == 1 ? "LM" : i == 2 ? "VM (BFGS)" : "VM (DFP)", zVecDist( var, ans ) );
      result2 = false;
    }
  }
  zOptDMDestroy( &opt );
  zVecFreeAtOnce( 4, var, ans, x0, g0 );
  zAssert( zOptDM (fixed step), result1 );
  zAssert( zOptDMSolve (fixed step), result2 );
}

/* scaled quadratic function: f(y) = k q(y/s), where q is the quadratic function above and (k,s) is given via util */
double quad_scaled_eval(const zVec var, void *util)
{
  double *ks = util;
  int i;
  double val = 0;

  for( i=0; i<zVecSizeNC(var); i++ )
    val += (i+1)*zSqr( zVecElemNC(var,i)/ks[1] - i );
  return 0.5*ks[0]*val;
}

zVec quad_scaled_grad(const zVec var, zVec g, void *util)
{
  double *ks = util;
  int i;

  for( i=0; i<zVecSizeNC(var); i++ )
    zVecSetElemNC( g, i, ks[0]/ks[1]*(i+1)*( zVecElemNC(var,i)/ks[1] - i ) );
  return g;
}

/* a function that cannot be evaluated except at the origin */
double nan_eval(const zVec var, void *dummy)
{
  int i;

  for( i=0; i<zVecSizeNC(var); i++ )
    if( zVecElemNC(var,i) != 0 ) return NAN;
  return 0;
}

zVec nan_grad(const zVec var, zVec g, void *dummy)
{
  return zVecSetAll( g, 1.0 );
}

/* the step has to be judged small relatively to the scale of the initial guess even for a small variable,
 * and the step has to be invariant with respect to the scale of the variable. */
void assert_opt_dm_step_scale(void)
{
  zOptDM opt;
  zVec var, ans;
  double ks[] = { 1.0, 0 };
  const char *label[] = { "SD (MT)", "LM (MT)", "VM (BFGS, MT)", "CG (MT)", "CG (Brent)" };
  int i, k;
  bool result = true;
  const double tol = 1.0e-6;

  var = zVecAlloc( DIM );
  ans = zVecAlloc( DIM );
  zOptDMCreate( &opt, DIM, quad_scaled_eval, quad_scaled_grad, NULL );
  zOptDMSetTolGrad( &opt, 0 ); /* terminate only by the step */
  for( ks[1]=1.0e-8; ks[1]<=1.0e8; ks[1]*=1.0e4 ){
    for( i=0; i<DIM; i++ ) zVecSetElemNC( ans, i, i*ks[1] );
    for( k=0; k<5; k++ ){
      switch( k ){
      case 0:  zOptDMAssignSD( &opt, "MT", 0 ); break;
      case 1:  zOptDMAssignLM( &opt, NULL, "MT", 0 ); break;
      case 2:  zOptDMAssignVM( &opt, "BFGS", "MT", 0 ); break;
      case 3:  zOptDMAssignCG( &opt, "MT", 0 ); break;
      default: zOptDMAssignCG( &opt, "Brent", 0 );
      }
      zVecSetAll( var, DIM*ks[1] );
      if( zOptDMSolve( &opt, var, ks, zTOL, 0, NULL ) < 0 || zVecDist( var, ans ) > tol * ks[1] ){
        eprintf( "%s, scale = %g: |x-x*|/s = %g\n", label[k], ks[1], zVecDist( var, ans ) / ks[1] );
        result = false;
      }
    }
  }
  zOptDMDestroy( &opt );
  zVecFreeAtOnce( 2, var, ans );
  zAssert( zOptDMSolve (termination by the step relative to the initial guess), result );
}

/* the threshold of the gradient vector relative to the initial one has to be adjustable, and the criterion has to be disabled with zero. */
void assert_opt_dm_tol_grad(void)
{
  zOptDM opt;
  zVec var;
  double ks[] = { 1.0, 1.0 };
  const double tol_grad[] = { 1.0e-3, ZM_OPT_DM_GRAD_MIN, 0 };
  double g0, g[3];
  int i, n[3];
  bool result1, result2 = true;

  var = zVecAlloc( DIM );
  zOptDMCreate( &opt, DIM, quad_scaled_eval, quad_scaled_grad, NULL );
  result1 = opt.tol_grad == ZM_OPT_DM_GRAD_MIN;
  zOptDMAssignSD( &opt, "MT", 0 ); /* linear convergence */
  zVecSetAll( var, 1.0 );
  g0 = zVecElemAbsMax( quad_scaled_grad( var, opt._g, ks ), NULL );
  for( i=0; i<3; i++ ){
    zOptDMSetTolGrad( &opt, tol_grad[i] );
    zVecSetAll( var, 1.0 );
    n[i] = zOptDMSolve( &opt, var, ks, zTOL, 0, NULL );
    g[i] = zVecElemAbsMax( quad_scaled_grad( var, opt._g, ks ), NULL );
    if( n[i] <= 0 || ( tol_grad[i] > 0 && g[i] > tol_grad[i] * g0 ) ){
      eprintf( "tol_grad = %g: returned %d, |g|/|g0| = %g\n", tol_grad[i], n[i], g[i] / g0 );
      result2 = false;
    }
  }
  if( !( n[0] < n[1] && n[1] < n[2] && g[2] < g[1] ) ){
    eprintf( "number of iterations: %d, %d, %d, |g|/|g0| = %g, %g, %g\n", n[0], n[1], n[2], g[0]/g0, g[1]/g0, g[2]/g0 );
    result2 = false;
  }
  zOptDMDestroy( &opt );
  zVecFree( var );
  zAssert( zOptDM (default threshold of the gradient vector), result1 );
  zAssert( zOptDMSolve (threshold of the gradient vector), result2 );
}

/* the iteration has to terminate with the gradient relatively small to the initial one, and has to fail if the line search fails. */
void assert_opt_dm_termination(void)
{
  zOptDM opt;
  zVec var, ans;
  double ks[] = { 1.0e-6, 1.0e3 }; /* a small and flat function, which the original absolute criterion misjudged as converged */
  const char *stepmethod[] = { "fixed", "MT", "DS" };
  int i;
  bool result1 = true, result2 = true;
  const double tol = 1.0e-6;

  var = zVecAlloc( DIM );
  ans = zVecAlloc( DIM );
  for( i=0; i<DIM; i++ ) zVecSetElemNC( ans, i, i*ks[1] );
  zOptDMCreate( &opt, DIM, quad_scaled_eval, quad_scaled_grad, NULL );
  for( i=0; i<3; i++ ){
    if( i == 0 )
      zOptDMAssignSD( &opt, "MT", 0 );
    else if( i == 1 )
      zOptDMAssignLM( &opt, NULL, "MT", 0 );
    else
      zOptDMAssignVM( &opt, "BFGS", "MT", 0 );
    zVecSetAll( var, ks[1] );
    if( zOptDMSolve( &opt, var, ks, zTOL, 0, NULL ) < 0 || zVecDist( var, ans ) > tol * ks[1] ){
      eprintf( "%s (MT), small and flat function: |x-x*|/s = %g\n", i == 0 ? "SD" : i == 1 ? "LM" : "VM", zVecDist( var, ans ) / ks[1] );
      result1 = false;
    }
  }
  zOptDMDestroy( &opt );
  /* failure of line search */
  zOptDMCreate( &opt, DIM, nan_eval, nan_grad, NULL );
  zEchoOff(); /* to suppress error messages of line search */
  for( i=0; i<3; i++ ){
    zOptDMAssignSD( &opt, stepmethod[i], 0 );
    zVecZero( var );
    if( zOptDMSolve( &opt, var, NULL, zTOL, 0, NULL ) != -1 ){
      eprintf( "SD (%s): line search failure not detected\n", stepmethod[i] );
      result2 = false;
    }
  }
  zEchoOn();
  zOptDMDestroy( &opt );
  zVecFreeAtOnce( 2, var, ans );
  zAssert( zOptDMSolve (termination by relative gradient), result1 );
  zAssert( zOptDMSolve (failure of line search), result2 );
}

/* the returned value has to be the number of iterations. */
void assert_opt_dm_return(void)
{
  zOptDM opt;
  zVec var;
  double k = 1.0;
  int i, n;
  bool result = true;

  var = zVecAlloc( DIM );
  zOptDMCreate( &opt, DIM, quad_eval, quad_grad, NULL );
  zOptDMAssignVM( &opt, "BFGS", "MT", 0 );
  /* no iteration at the minimum */
  for( i=0; i<DIM; i++ ) zVecSetElemNC( var, i, i );
  if( zOptDMSolve( &opt, var, &k, zTOL, 0, NULL ) != 0 ){
    eprintf( "nonzero iterations at the minimum\n" );
    result = false;
  }
  /* the iteration converged within n iterations has to be converged also with the upper limit n */
  zVecSetAll( var, 1.0 );
  if( ( n = zOptDMSolve( &opt, var, &k, zTOL, 0, NULL ) ) <= 0 ){
    eprintf( "not converged\n" );
    result = false;
  } else{
    zVecSetAll( var, 1.0 );
    zEchoOff();
    if( zOptDMSolve( &opt, var, &k, zTOL, n, NULL ) != n ||
        ( zVecSetAll( var, 1.0 ), zOptDMSolve( &opt, var, &k, zTOL, n-1, NULL ) ) != -1 ){
      eprintf( "the returned value differs from the number of iterations\n" );
      result = false;
    }
    zEchoOn();
  }
  zOptDMDestroy( &opt );
  zVecFree( var );
  zAssert( zOptDMSolve (returned number of iterations), result );
}

/* quadratic function that records whether it is evaluated outside the range limited by the scale from x0 */
typedef struct{
  zVec x0;
  double scale;
  bool outside;
} scale_check_t;

double quad_scale_check_eval(const zVec var, void *util)
{
  scale_check_t *sc = util;
  int i;

  for( i=0; i<zVecSizeNC(var); i++ )
    if( fabs( zVecElemNC(var,i) - zVecElemNC(sc->x0,i) ) > sc->scale * ( 1 + zTOL ) ) sc->outside = true;
  return quad_eval( var, NULL );
}

zVec quad_scale_check_grad(const zVec var, zVec g, void *util)
{
  return quad_grad( var, g, NULL );
}

/* the line search must not evaluate the function outside the range limited by the scale. */
void assert_opt_dm_scale(void)
{
  zOptDM opt;
  zVec var;
  scale_check_t sc;
  const char *stepmethod[] = { "MT", "DS", "Brent" };
  int i;
  bool result = true;

  var = zVecAlloc( DIM );
  sc.x0 = zVecAlloc( DIM );
  sc.scale = 0.1; /* smaller than the steepest descent step */
  zOptDMCreate( &opt, DIM, quad_scale_check_eval, quad_scale_check_grad, NULL );
  zEchoOff(); /* to suppress a warning of unconverged iteration */
  for( i=0; i<3; i++ ){
    zOptDMAssignSD( &opt, stepmethod[i], sc.scale );
    zVecSetAll( sc.x0, 1.0 );
    zVecCopy( sc.x0, var );
    sc.outside = false;
    zOptDMSolve( &opt, var, &sc, zTOL, 1, NULL ); /* one step */
    if( sc.outside ){
      eprintf( "SD (%s): evaluated outside the range of the scale\n", stepmethod[i] );
      result = false;
    }
  }
  zEchoOn();
  zOptDMDestroy( &opt );
  zVecFreeAtOnce( 2, var, sc.x0 );
  zAssert( zOptDM (evaluation within the range of the scale), result );
}

/* a function from the tutorial of Mathematica, which has a local minimum at (1.37638,1.67868) */
double mma_eval(const zVec var, void *dummy)
{
  int i;
  double x, y, val = 0;

  for( i=0; i<zVecSizeNC(var); i+=2 ){
    x = zVecElemNC(var,i);
    y = zVecElemNC(var,i+1);
    val += cos(x*x-3*y) + sin(x*x+y*y);
  }
  return val;
}

zVec mma_grad(const zVec var, zVec g, void *dummy)
{
  int i;
  double x, y;

  for( i=0; i<zVecSizeNC(var); i+=2 ){
    x = zVecElemNC(var,i);
    y = zVecElemNC(var,i+1);
    zVecSetElemNC( g, i,   2*x*( cos(x*x+y*y) - sin(x*x-3*y) ) );
    zVecSetElemNC( g, i+1, 3*sin(x*x-3*y) + 2*y*cos(x*x+y*y) );
  }
  return g;
}

/* the iteration has to be regarded as converged if no further decrease of the function is detected near the minimum. */
void assert_opt_dm_stay(void)
{
  zOptDM opt;
  zVec var, g;
  bool result = true;
  const double tol = 1.0e-6;

  var = zVecAlloc( DIM );
  g = zVecAlloc( DIM );
  zOptDMCreate( &opt, DIM, mma_eval, mma_grad, NULL );
  /* the original implementation of Dennis-Schnabel's method failed after reaching the minimum. */
  zOptDMAssignLM( &opt, NULL, "DS", 0 );
  zVecSetAll( var, 1.0 );
  if( zOptDMSolve( &opt, var, NULL, zTOL, 0, NULL ) < 0 || zVecElemAbsMax( mma_grad( var, g, NULL ), NULL ) > tol ){
    eprintf( "LM (DS): |g| = %g\n", zVecElemAbsMax( mma_grad( var, g, NULL ), NULL ) );
    result = false;
  }
  zOptDMDestroy( &opt );
  zVecFreeAtOnce( 2, var, g );
  zAssert( zOptDMSolve (convergence without further decrease), result );
}

/* conjugate gradient method has to compute Fletcher-Reeves' coefficient with any line search, and converge. */
void assert_opt_dm_cg(void)
{
  zOptDM opt;
  zVec var, ans, g0, g1;
  const char *stepmethod[] = { "Brent", "MT", "DS" };
  int i, n;
  bool result1 = true, result2 = true;
  const double tol = 1.0e-6;

  var = zVecAlloc( DIM );
  ans = zVecAlloc( DIM );
  g0 = zVecAlloc( DIM );
  g1 = zVecAlloc( DIM );
  for( i=0; i<DIM; i++ ) zVecSetElemNC( ans, i, i );
  zOptDMCreate( &opt, DIM, quad_eval, quad_grad, NULL );
  for( i=0; i<3; i++ ){
    zOptDMAssignCG( &opt, stepmethod[i], 0 );
    /* Fletcher-Reeves' coefficient after one iteration */
    zVecSetAll( var, 1.0 );
    quad_grad( var, g0, NULL );
    zEchoOff(); /* to suppress a warning of unconverged iteration */
    zOptDMSolve( &opt, var, NULL, zTOL, 1, NULL );
    zEchoOn();
    quad_grad( var, g1, NULL );
    if( !zEqual( opt._b, zVecSqrNorm( g1 ) / zVecSqrNorm( g0 ), zTOL ) ){
      eprintf( "CG (%s): coefficient = %g (to be %g)\n", stepmethod[i], opt._b, zVecSqrNorm( g1 ) / zVecSqrNorm( g0 ) );
      result1 = false;
    }
    /* convergence. with the strict line search, it converges within the dimension for a quadratic function. */
    zVecSetAll( var, 1.0 );
    if( ( n = zOptDMSolve( &opt, var, NULL, zTOL, 0, NULL ) ) < 0 || zVecDist( var, ans ) > tol || ( i == 0 && n > DIM ) ){
      eprintf( "CG (%s): %d iterations, |x-x*| = %g\n", stepmethod[i], n, zVecDist( var, ans ) );
      result2 = false;
    }
  }
  zOptDMDestroy( &opt );
  /* Rosenbrock function */
  zVecSetAll( ans, 1.0 );
  zOptDMCreate( &opt, DIM, rosenbrock_eval, rosenbrock_grad, NULL );
  zOptDMSetTolGrad( &opt, 1.0e-10 ); /* see assert_opt_dm_lm */
  zOptDMAssignCG( &opt, "MT", 0 );
  zVecSetElemList( var, -1.2, 1.0, -1.2, 1.0 ); /* the standard initial guess */
  if( zOptDMSolve( &opt, var, NULL, zTOL, 0, NULL ) < 0 || zVecDist( var, ans ) > tol ){
    eprintf( "CG (MT), Rosenbrock: |x-x*| = %g\n", zVecDist( var, ans ) );
    result2 = false;
  }
  zOptDMDestroy( &opt );
  zVecFreeAtOnce( 4, var, ans, g0, g1 );
  zAssert( zOptDM (coefficient of Fletcher-Reeves formula), result1 );
  zAssert( zOptDMSolve (conjugate gradient method), result2 );
}

/* Brent's line search has to bracket a minimum regardless of the scale of the function, and must not increase the function. */
void assert_opt_dm_brent(void)
{
  zOptDM opt;
  zVec var, ans;
  double k[] = { 1.0e-3, 1.0, 1.0e3, 0 }, *kp;
  int i, n;
  bool result = true;
  const double tol = 1.0e-6;

  var = zVecAlloc( DIM );
  ans = zVecAlloc( DIM );
  /* scaled quadratic functions, for which the original search within [0,1] did not reach the minimizer */
  for( i=0; i<DIM; i++ ) zVecSetElemNC( ans, i, i );
  zOptDMCreate( &opt, DIM, quad_eval, quad_grad, NULL );
  zOptDMAssignCG( &opt, "Brent", 0 );
  for( kp=k; *kp>0; kp++ ){
    zVecSetAll( var, 1.0 );
    if( ( n = zOptDMSolve( &opt, var, kp, zTOL, 0, NULL ) ) < 0 || n > DIM || zVecDist( var, ans ) > tol ){
      eprintf( "CG (Brent), k=%g: %d iterations, |x-x*| = %g\n", *kp, n, zVecDist( var, ans ) );
      result = false;
    }
  }
  zOptDMDestroy( &opt );
  /* Rosenbrock function, where the original search converged to a local minimum along the line higher than the initial point */
  zVecSetAll( ans, 1.0 );
  zOptDMCreate( &opt, DIM, rosenbrock_eval, rosenbrock_grad, NULL );
  zOptDMAssignCG( &opt, "Brent", 0 );
  zVecSetElemList( var, -1.2, 1.0, -1.2, 1.0 );
  if( zOptDMSolve( &opt, var, NULL, zTOL, 0, NULL ) < 0 || zVecDist( var, ans ) > tol ){
    eprintf( "CG (Brent), Rosenbrock: |x-x*| = %g\n", zVecDist( var, ans ) );
    result = false;
  }
  zOptDMDestroy( &opt );
  zVecFreeAtOnce( 2, var, ans );
  zAssert( zOptDMSolve (Brent), result );
}

/* Levenberg-Marquardt's method has to converge regardless of the scale of the function with both damping coefficients, even with the fixed step. */
void assert_opt_dm_lm(void)
{
  zOptDM opt;
  zVec var, ans;
  double k[] = { 1.0e-3, 1.0, 1.0e3, 0 }, *kp;
  const char *dampmethod[] = { "grad", "dec" };
  const char *stepmethod[] = { "fixed", "DS", "MT" };
  int i, j;
  bool result = true;
  const double tol = 1.0e-6;

  var = zVecAlloc( DIM );
  ans = zVecAlloc( DIM );
  /* scaled quadratic functions, for which the original damping coefficient |g|^2 stagnated */
  for( i=0; i<DIM; i++ ) zVecSetElemNC( ans, i, i );
  zOptDMCreate( &opt, DIM, quad_eval, quad_grad, NULL );
  for( kp=k; *kp>0; kp++ )
    for( i=0; i<2; i++ )
      for( j=0; j<3; j++ ){
        zOptDMAssignLM( &opt, dampmethod[i], stepmethod[j], 0 );
        zVecSetAll( var, 1.0 );
        if( zOptDMSolve( &opt, var, kp, zTOL, 0, NULL ) < 0 || zVecDist( var, ans ) > tol ){
          eprintf( "LM (%s, %s), k=%g: |x-x*| = %g\n", dampmethod[i], stepmethod[j], *kp, zVecDist( var, ans ) );
          result = false;
        }
      }
  zOptDMDestroy( &opt );
  /* Rosenbrock function */
  zVecSetAll( ans, 1.0 );
  zOptDMCreate( &opt, DIM, rosenbrock_eval, rosenbrock_grad, NULL );
  /* the default threshold of the gradient vector ZM_OPT_DM_GRAD_MIN |g0| = 2.2e-6 could leave the error
   * up to 5e-6 with the minimum eigenvalue 0.4 of the Hessian matrix at the solution. */
  zOptDMSetTolGrad( &opt, 1.0e-10 );
  for( i=0; i<2; i++ )
    for( j=1; j<3; j++ ){
      zOptDMAssignLM( &opt, dampmethod[i], stepmethod[j], 0 );
      zVecSetElemList( var, -1.2, 1.0, -1.2, 1.0 );
      if( zOptDMSolve( &opt, var, NULL, zTOL, 0, NULL ) < 0 || zVecDist( var, ans ) > tol ){
        eprintf( "LM (%s, %s), Rosenbrock: |x-x*| = %g\n", dampmethod[i], stepmethod[j], zVecDist( var, ans ) );
        result = false;
      }
    }
  zOptDMDestroy( &opt );
  zVecFreeAtOnce( 2, var, ans );
  zAssert( zOptDMSolve (Levenberg-Marquardt), result );
}

/* a smooth function f(x) = sum_i (i+1) sin(x_i/s) exp(0.3 x_{i+1}/s) with the length scale s given via util */
double smooth_eval(const zVec var, void *util)
{
  double s = *(double *)util, val = 0;
  int i;

  for( i=0; i<DIM; i++ )
    val += (i+1) * sin( zVecElemNC(var,i)/s ) * exp( 0.3*zVecElemNC(var,(i+1)%DIM)/s );
  return val;
}

zVec smooth_grad(const zVec var, zVec g, void *util)
{
  double s = *(double *)util, e;
  int i, j;

  zVecZero( g );
  for( i=0; i<DIM; i++ ){
    j = ( i + 1 ) % DIM;
    e = exp( 0.3*zVecElemNC(var,j)/s );
    zVecElemNC(g,i) += (i+1) * cos( zVecElemNC(var,i)/s ) / s * e;
    zVecElemNC(g,j) += (i+1) * sin( zVecElemNC(var,i)/s ) * 0.3 / s * e;
  }
  return g;
}

zMat smooth_hess(const zVec var, zMat h, void *util)
{
  double s = *(double *)util, c, e, si, ci;
  int i, j;

  zMatZero( h );
  for( i=0; i<DIM; i++ ){
    j = ( i + 1 ) % DIM;
    c = (i+1) / ( s*s );
    e = exp( 0.3*zVecElemNC(var,j)/s );
    si = sin( zVecElemNC(var,i)/s );
    ci = cos( zVecElemNC(var,i)/s );
    zMatElemNC(h,i,i) -= c * si * e;
    zMatElemNC(h,i,j) += c * ci * 0.3 * e;
    zMatElemNC(h,j,i) += c * ci * 0.3 * e;
    zMatElemNC(h,j,j) += c * si * 0.09 * e;
  }
  return h;
}

/* the numerical gradient vector and Hessian matrix have to be accurate regardless of the scale of the variable. */
void assert_opt_dm_numdiff(void)
{
  zOptDM opt1, opt2;
  zVec var, ga, gn;
  zMat ha, h1, h2;
  double s[] = { 1.0e-6, 1.0e-3, 1.0, 1.0e3, 1.0e6, 0 }, *sp;
  int i, j;
  bool result1 = true, result2 = true, result3 = true, result4 = true;
  const double tol1 = 1.0e-8, tol2 = 1.0e-5;

  var = zVecAlloc( DIM );
  ga = zVecAlloc( DIM );
  gn = zVecAlloc( DIM );
  ha = zMatAllocSqr( DIM );
  h1 = zMatAllocSqr( DIM );
  h2 = zMatAllocSqr( DIM );
  zOptDMCreate( &opt1, DIM, smooth_eval, smooth_grad, NULL ); /* Hessian matrix from the gradient vector */
  zOptDMCreate( &opt2, DIM, smooth_eval, NULL, NULL );        /* gradient vector and Hessian matrix from the function value */
  for( sp=s; *sp>0; sp++ ){
    for( j=0; j<2; j++ ){
      for( i=0; i<DIM; i++ ) zVecSetElemNC( var, i, *sp*(0.3+0.2*i) );
      if( j == 1 ) zVecSetElemNC( var, 2, 1.0e-12 ); /* a component near zero */
      smooth_grad( var, ga, sp );
      smooth_hess( var, ha, sp );
      opt2._grad( &opt2, var, gn, sp );
      opt1._hess( &opt1, var, h1, sp );
      opt2._hess( &opt2, var, h2, sp );
      if( zVecDist( gn, ga ) > tol1 * zVecNorm( ga ) ){
        eprintf( "numerical gradient, s=%g: relative error = %g\n", *sp, zVecDist( gn, ga ) / zVecNorm( ga ) );
        result1 = false;
      }
      zMatSubDRC( h1, ha );
      if( zMatNorm( h1 ) > tol1 * zMatNorm( ha ) ){
        eprintf( "numerical Hessian from gradient, s=%g: relative error = %g\n", *sp, zMatNorm( h1 ) / zMatNorm( ha ) );
        result2 = false;
      }
      if( !zMatIsSymmetric( h2 ) ) result4 = false;
      zMatSubDRC( h2, ha );
      if( zMatNorm( h2 ) > tol2 * zMatNorm( ha ) ){
        eprintf( "numerical Hessian from function, s=%g: relative error = %g\n", *sp, zMatNorm( h2 ) / zMatNorm( ha ) );
        result3 = false;
      }
    }
  }
  opt1._hess( &opt1, var, h1, sp-1 );
  if( !zMatIsSymmetric( h1 ) ) result4 = false;
  /* the ratio of the step can be set */
  zOptDMSetEPSMin( &opt2, 1.0e-3 );
  if( opt2.eps_min != 1.0e-3 ) result4 = false;
  zOptDMDestroy( &opt1 );
  zOptDMDestroy( &opt2 );
  zVecFreeAtOnce( 3, var, ga, gn );
  zMatFreeAtOnce( 3, ha, h1, h2 );
  zAssert( zOptDM (numerical gradient), result1 );
  zAssert( zOptDM (numerical Hessian from gradient), result2 );
  zAssert( zOptDM (numerical Hessian from function), result3 );
  zAssert( zOptDM (symmetry of numerical Hessian and the ratio of the step), result4 );
}

int main(void)
{
  zRandInit();
  assert_opt_dm_eval();
  assert_opt_dm_ds();
  assert_opt_dm_update();
  assert_opt_dm_descent();
  assert_opt_dm_mt_wolfe();
  assert_opt_dm_mt();
  assert_opt_dm_fixed();
  assert_opt_dm_termination();
  assert_opt_dm_tol_grad();
  assert_opt_dm_step_scale();
  assert_opt_dm_return();
  assert_opt_dm_scale();
  assert_opt_dm_stay();
  assert_opt_dm_cg();
  assert_opt_dm_brent();
  assert_opt_dm_lm();
  assert_opt_dm_numdiff();
  return 0;
}
