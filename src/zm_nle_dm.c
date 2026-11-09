/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_nle_dm - nonlinear equation:
 * simultaneous nonlinear equation solver based on descent method.
 */

#include <zm/zm_nle.h>

static zMat _zNLEJacobi(zNLE *nle, const zVec var, zMat j, void *util){
  return nle->jac( var, j, util );
}
/* numerical Jacobian matrix by the central difference with the step relative to the magnitude of the variable. */
static zMat _zNLEJacobiNG(zNLE *nle, const zVec var, zMat j, void *util)
{
  int i;
  double org, h, xmax;

  xmax = zVecElemAbsMax( var, NULL );
  for( i=0; i<zVecSizeNC(var); i++ ){
    h = zOptDMDiffStep( var, i, xmax, nle->_opt.eps_min, ZM_OPT_DM_DIFF_SCALE_RATIO );
    org = zVecElemNC(var,i);
    zVecSetElemNC( var, i, org + h );
    nle->f( var, nle->_adg, util );
    zVecSetElemNC( var, i, org - h );
    nle->f( var, nle->_prg, util );
    zVecSetElemNC( var, i, org );
    zVecSubDRC( nle->_adg, nle->_prg );
    zMatPutColNC( j, i, zVecMulDRC( nle->_adg, 0.5 / h ) ); /* the i-th column is the derivative with respect to the i-th variable */
  }
  return j;
}

#define _zm_nle(u) ( (zNLE *)u )

static double _zNLEEval(const zVec var, void *util)
{
  _zm_nle(util)->f( var, _zm_nle(util)->_f, _zm_nle(util)->util );
  zVecAmpNC( _zm_nle(util)->_f, _zm_nle(util)->we, _zm_nle(util)->_fw );
  return 0.5 * zVecInnerProd( _zm_nle(util)->_f, _zm_nle(util)->_fw );
}

static zVec _zNLEGrad(const zVec var, zVec grad, void *util)
{
  _zm_nle(util)->_jac( _zm_nle(util), var, _zm_nle(util)->_j, _zm_nle(util)->util );
  return zMulMatTVecNC( _zm_nle(util)->_j, _zm_nle(util)->_fw, grad );
}

static zMat _zNLEHess(const zVec var, zMat h, void *util)
{
  return zMatTQuadNC( _zm_nle(util)->_j, _zm_nle(util)->we, h );
}

#undef _zm_nle

/* check if the step a d is small relatively to the variable. */
static bool _zNLEStepIsSmall(zNLE *nle, const zVec var, double a, double tol_step)
{
  return a * zVecElemAbsMax( nle->_opt._d, NULL ) <= tol_step * ( 1 + zVecElemAbsMax( var, NULL ) );
}

static void _zNLEEvalRes(zNLE *nle, const zVec var, void *util, double *err)
{
  if( !err ) return;
  nle->util = util;
  *err = _zNLEEval( var, nle );
}

/* the iteration succeeds if the residual becomes small, namely, |f|inf <= tol_err, and fails if
 * - the step becomes small relatively to the variable, namely, |p|inf <= tol_step ( 1 + |x|inf ),
 * - the number of iterations exceeds the given upper limit iter, or
 * - the function cannot be evaluated along the step.
 */

/* descent method. the criterion of zOptDMSolve() with the gradient vector is disabled in zNLECreate(). */
static int _zNLESolveDM(zNLE *nle, zVec var, void *util, double tol_step, double tol_err, int iter, double *err)
{
  int ret;

  nle->f( var, nle->_f, util );
  if( zVecIsTol( nle->_f, tol_err ) ){
    _zNLEEvalRes( nle, var, util, err );
    return 0;
  }
  nle->util = util;
  ret = zOptDMSolve( &nle->_opt, var, nle, tol_step, iter, err );
  nle->f( var, nle->_f, util );
  return ret >= 0 && zVecIsTol( nle->_f, tol_err ) ? ret : -1;
}

/* Newton-Raphson's method. */

#define ZM_NLE_NR_ARMIJO_COEFF    ( 1.0e-4 )
#define ZM_NLE_NR_BACKTRACK_RATIO ( 0.5 )

/* a half of the squared norm of the residual vector f weighted by we. */
static double _zNLEResNorm(zNLE *nle, const zVec f)
{
  zVecAmpNC( f, nle->we, nle->_fw );
  return 0.5 * zVecInnerProd( f, nle->_fw );
}

/* damping vector of Newton-Raphson's method, which is wn scaled by nr_damp_coeff times the maximum
 * diagonal component of J^T diag(we) J, so that it is invariant with respect to the scales of the
 * variable and the residual. */
static zVec _zNLENRDamp(zNLE *nle, zVec damp)
{
  int i, j;
  double d, dmax = 0;

  for( j=0; j<zMatColSizeNC(nle->_j); j++ ){
    for( d=0, i=0; i<zMatRowSizeNC(nle->_j); i++ )
      d += zVecElemNC(nle->we,i) * zSqr( zMatElemNC(nle->_j,i,j) );
    if( d > dmax ) dmax = d;
  }
  return zVecMulNC( nle->wn, nle->nr_damp_coeff * dmax, damp );
}

/* update of the Jacobian matrix at the next point nle->_opt._x from the current point var, which returns
 * the true value if the updated Jacobian matrix is exact. */

/* Newton-Raphson's method: the exact Jacobian matrix. */
static bool _zNLEUpdateJacobiNR(zNLE *nle, const zVec var, void *util)
{
  nle->_jac( nle, nle->_opt._x, nle->_j, util );
  return true;
}

/* Broyden's method: J += ( y - J s ) s^T / s^T s, where s is the step and y = f(x+s) - f(x). */
static bool _zNLEUpdateJacobiBroyden(zNLE *nle, const zVec var, void *util)
{
  zVecSubNC( nle->_opt._x, var, nle->_opt._p ); /* s */
  zMulMatVecNC( nle->_j, nle->_opt._p, nle->_fw );
  zVecSubNC( nle->_fp, nle->_fw, nle->_fw );
  zVecSubNCDRC( nle->_fw, nle->_f ); /* y - J s */
  zMatCatDyadNC( nle->_j, 1.0/zVecSqrNorm(nle->_opt._p), nle->_fw, nle->_opt._p );
  return false;
}

/* the descent direction d is computed by the SR-inverse of the Jacobian matrix J, namely,
 * ( J^T W J + diag(damp) ) d = J^T W f, and the step -a d is chosen by the backtracking line search
 * to satisfy Armijo's condition on the squared norm of the residual. J is updated by update_jac.
 * if J is not exact, i.e., approximated by Broyden's rule, and it fails to give a step that decreases
 * the residual or the step stagnates, J is recomputed exactly. */
static int _zNLESolveNewton(zNLE *nle, zVec var, void *util, double tol_step, double tol_err, int iter, double *err, bool (* update_jac)(zNLE*,const zVec,void*))
{
  int i, ret = -1;
  double a, e, ep, slope;
  bool is_exact; /* true if J is the exact Jacobian matrix at var */
  bool terminate = false;
  zLEWorkspace workspace;

  workspace.m_regular = nle->_opt._h;
  workspace.v_ispace = nle->_opt._p;
  workspace.index_le = nle->_opt._idx;
  workspace.v_scale = nle->_opt._q;
  nle->f( var, nle->_f, util );
  if( zVecIsTol( nle->_f, tol_err ) ){
    ret = 0;
    goto TERMINATE;
  }
  e = _zNLEResNorm( nle, nle->_f );
  nle->_jac( nle, var, nle->_j, util );
  is_exact = true;
  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
   RETRY:
    _zNLENRDamp( nle, nle->_opt._r );
    zVecCopyNC( nle->_f, nle->_fw ); /* since zLESolveSRDST destroys the right-hand-side vector */
    if( !zLESolveSRDST( nle->_j, nle->_fw, nle->_opt._r, nle->we, nle->_opt._d, &workspace ) )
      goto RESTART;
    /* directional derivative of the squared norm of the residual along -d */
    zVecAmpNC( nle->_f, nle->we, nle->_fw );
    zMulMatTVecNC( nle->_j, nle->_fw, nle->_opt._g );
    slope = zVecInnerProd( nle->_opt._g, nle->_opt._d );
    for( a=1.0; ; a*=ZM_NLE_NR_BACKTRACK_RATIO ){
      zVecCatNC( var, -a, nle->_opt._d, nle->_opt._x );
      nle->f( nle->_opt._x, nle->_fp, util );
      if( !zVecIsNan( nle->_fp ) &&
          ( ep = _zNLEResNorm( nle, nle->_fp ) ) <= e - ZM_NLE_NR_ARMIJO_COEFF * a * slope ) break;
      if( _zNLEStepIsSmall( nle, var, a, tol_step ) ) goto RESTART; /* line search fails */
    }
    /* the Jacobian matrix is updated before var and f are updated, but not if the iteration terminates */
    if( zVecIsTol( nle->_fp, tol_err ) ){
      ret = i + 1; /* succeed */
      terminate = true;
    } else
    if( !_zNLEStepIsSmall( nle, var, a, tol_step ) ){
      is_exact = update_jac( nle, var, util );
    } else
    if( is_exact ){
      terminate = true; /* stagnate */
    } else{ /* the stagnation might be due to the approximate Jacobian matrix */
      nle->_jac( nle, nle->_opt._x, nle->_j, util );
      is_exact = true;
    }
    zVecCopyNC( nle->_opt._x, var );
    zVecCopyNC( nle->_fp, nle->_f );
    e = ep;
    if( terminate ) goto TERMINATE;
    continue;
   RESTART:
    if( is_exact ) goto TERMINATE; /* fails even with the exact Jacobian matrix */
    /* retry with the exact Jacobian matrix if the approximate one does not give a descent step */
    nle->_jac( nle, var, nle->_j, util );
    is_exact = true;
    goto RETRY;
  }
  ZITERWARN( iter );
 TERMINATE:
  _zNLEEvalRes( nle, var, util, err );
  return ret;
}

/* Newton-Raphson's method. */
static int _zNLESolveNR(zNLE *nle, zVec var, void *util, double tol_step, double tol_err, int iter, double *err)
{
  return _zNLESolveNewton( nle, var, util, tol_step, tol_err, iter, err, _zNLEUpdateJacobiNR );
}

/* Broyden's method. */
static int _zNLESolveBroyden(zNLE *nle, zVec var, void *util, double tol_step, double tol_err, int iter, double *err)
{
  return _zNLESolveNewton( nle, var, util, tol_step, tol_err, iter, err, _zNLEUpdateJacobiBroyden );
}

/* create a simultaneous nonlinear equations solver. */
zNLE *zNLECreate(zNLE *nle, int nv, int ne, zVec (* f)(const zVec,zVec,void*), zMat (* jac)(const zVec,zMat,void*))
{
  nle->util = NULL;
  nle->wn = zVecAlloc( nv );
  nle->we = zVecAlloc( ne );
  nle->_f = zVecAlloc( ne );
  nle->_fw = zVecAlloc( ne );
  nle->_fp = zVecAlloc( ne );
  nle->_j = zMatAlloc( ne, nv );
  if( !nle->wn || !nle->we || !nle->_f ||
      !nle->_fw || !nle->_fp || !nle->_j ||
      !zOptDMCreate( &nle->_opt, nv, _zNLEEval, _zNLEGrad, _zNLEHess ) ){
    zNLEDestroy( nle );
    return NULL;
  }
  /* the gradient vector can be small even far from the solution if the Jacobian matrix is nearly singular. */
  zOptDMSetTolGrad( &nle->_opt, 0 );
  nle->f = f; /* residual function */
  /* Jacobian matrix function */
  if( ( nle->jac = jac ) ){
    nle->_jac = _zNLEJacobi;
    nle->_adg = nle->_prg = NULL;
  } else{
    nle->_jac = _zNLEJacobiNG;
    nle->_adg = zVecAlloc( ne );
    nle->_prg = zVecAlloc( ne );
    if( !nle->_adg || !nle->_prg ){
      zNLEDestroy( nle );
      return NULL;
    }
  }
  zVecSetAll( nle->wn, 1.0 );
  zVecSetAll( nle->we, 1.0 );
  zNLESetNRDampCoeff( nle, ZM_NLE_NR_DEFAULT_DAMP_COEFF );
  return zNLEAssignNR( nle );
}

/* destroy a simultaneous nonlinear equations solver. */
void zNLEDestroy(zNLE *nle)
{
  nle->f = NULL;
  nle->jac = NULL;
  nle->util = NULL;
  zVecFreeAtOnce( 7, nle->wn, nle->we, nle->_f, nle->_fw, nle->_fp, nle->_adg, nle->_prg );
  zMatFree( nle->_j );
  zOptDMDestroy( &nle->_opt );
}

zNLE *zNLEAssignSD(zNLE *nle, const char *stepmethod, double scale){
  if( !zOptDMAssignSD( &nle->_opt, stepmethod, scale ) )
    return NULL;
  nle->_solve = _zNLESolveDM;
  return nle;
}
zNLE *zNLEAssignLM(zNLE *nle, const char *dampmethod, const char *stepmethod, double scale){
  if( !zOptDMAssignLM( &nle->_opt, dampmethod, stepmethod, scale ) )
    return NULL;
  nle->_solve = _zNLESolveDM;
  return nle;
}
zNLE *zNLEAssignVM(zNLE *nle, const char *updatemethod, const char *stepmethod, double scale){
  if( !zOptDMAssignVM( &nle->_opt, updatemethod, stepmethod, scale ) )
    return NULL;
  nle->_solve = _zNLESolveDM;
  return nle;
}
zNLE *zNLEAssignCG(zNLE *nle, const char *stepmethod, double scale){
  if( !zOptDMAssignCG( &nle->_opt, stepmethod, scale ) )
    return NULL;
  nle->_solve = _zNLESolveDM;
  return nle;
}
zNLE *zNLEAssignNR(zNLE *nle){
  nle->_solve = _zNLESolveNR;
  return nle;
}
zNLE *zNLEAssignBroyden(zNLE *nle){
  nle->_solve = _zNLESolveBroyden;
  return nle;
}
