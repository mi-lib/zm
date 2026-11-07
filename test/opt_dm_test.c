#include <zm/zm_opt.h>

#define DIM 4

/* quadratic function: f(x) = 1/2 sum_i (i+1)(x_i-i)^2, minimum at x_i = i */
double quad_eval(const zVec var, void *dummy)
{
  int i;
  double val = 0;

  for( i=0; i<zVecSizeNC(var); i++ )
    val += (i+1)*zSqr( zVecElemNC(var,i) - i );
  return 0.5*val;
}

zVec quad_grad(const zVec var, zVec g, void *dummy)
{
  int i;

  for( i=0; i<zVecSizeNC(var); i++ )
    zVecSetElemNC( g, i, (i+1)*( zVecElemNC(var,i) - i ) );
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

int main(void)
{
  assert_opt_dm_eval();
  return 0;
}
