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

zOptDM *assign_sd(zOptDM *opt, const char *stepmethod, const char *updatemethod){ return zOptDMAssignSD( opt, stepmethod ); }
zOptDM *assign_lm(zOptDM *opt, const char *stepmethod, const char *updatemethod){ return zOptDMAssignLM( opt, stepmethod ); }
zOptDM *assign_vm(zOptDM *opt, const char *stepmethod, const char *updatemethod){ return zOptDMAssignVM( opt, stepmethod, updatemethod ); }
zOptDM *assign_cg(zOptDM *opt, const char *stepmethod, const char *updatemethod){ return zOptDMAssignCG( opt ); }

method_t method[] = {
  { "SD (MT)",         assign_sd, "MT",    NULL },
  { "SD (DS)",         assign_sd, "DS",    NULL },
  { "SD (Noc)",        assign_sd, "Noc",   NULL },
  { "SD (Brent)",      assign_sd, "Brent", NULL },
  { "LM (MT)",         assign_lm, "MT",    NULL },
  { "LM (DS)",         assign_lm, "DS",    NULL },
  { "LM (Noc)",        assign_lm, "Noc",   NULL },
  { "VM (MT, BFGS)",   assign_vm, "MT",    "BFGS" },
  { "VM (DS, BFGS)",   assign_vm, "DS",    "BFGS" },
  { "VM (Noc, BFGS)",  assign_vm, "Noc",   "BFGS" },
  { "VM (Brent, BFGS)",assign_vm, "Brent", "BFGS" },
  { "CG",              assign_cg, NULL,    NULL },
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
  zOptDMCreate( &opt, DIM, 0, quad_eval, quad_grad, NULL );
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
  const double tol = 1.0e-8;

  var = zVecAlloc( DIM );
  ans = zVecAlloc( DIM );
  /* scaled quadratic functions */
  for( i=0; i<DIM; i++ ) zVecSetElemNC( ans, i, i );
  zOptDMCreate( &opt, DIM, 0, quad_eval, quad_grad, NULL );
  for( kp=k; *kp>0; kp++ ){
    for( i=0; i<2; i++ ){
      if( i == 0 )
        zOptDMAssignSD( &opt, "DS" );
      else
        zOptDMAssignVM( &opt, "DS", "BFGS" );
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
  zOptDMCreate( &opt, DIM, 0, rosenbrock_eval, rosenbrock_grad, NULL );
  zOptDMAssignVM( &opt, "DS", "BFGS" );
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
  zOptDMCreate( &opt, DIM, 0, quad_eval, quad_grad, NULL );
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
  zOptDMAssignVM( &opt, NULL, "DFP" );
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
  zOptDMAssignVM( &opt, NULL, "BFGS" );
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
  zOptDMAssignVM( &opt, NULL, "DFP" );
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
  const char *stepmethod[] = { "MT", "DS", "Noc", "Brent", NULL }, **sp;
  bool result = true;
  const double tol = 1.0e-8;

  var = zVecAlloc( DIM );
  zOptDMCreate( &opt, DIM, 0, saddle_eval, saddle_grad, saddle_hess );
  for( sp=stepmethod; *sp; sp++ ){
    zOptDMAssignLM( &opt, *sp );
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

int main(void)
{
  zRandInit();
  assert_opt_dm_eval();
  assert_opt_dm_ds();
  assert_opt_dm_update();
  assert_opt_dm_descent();
  return 0;
}
