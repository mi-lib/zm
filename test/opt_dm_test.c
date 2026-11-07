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

int main(void)
{
  assert_opt_dm_eval();
  assert_opt_dm_ds();
  return 0;
}
