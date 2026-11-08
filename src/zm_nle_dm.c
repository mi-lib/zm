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

/* update variable vector toward descent direction, and return the step length, or zero if it fails. */
static double _zNLEStep(zNLE *nle, zVec var, void *util)
{
  double a;

  for( a=1.0; a>zTOL; a*=0.5 ){
    zVecCatNC( var, -a, nle->_opt._d, nle->_opt._x );
    nle->f( nle->_opt._x, nle->_fp, util );
    if( !zVecIsNan( nle->_fp ) ){
      /* now, f( var ) is stored in _fp. if need it, copy it to _f. */
      zVecCopyNC( nle->_opt._x, var );
      return a;
    }
  }
  ZRUNWARN( ZM_WARN_OPT_BADINI ); /* dummy */
  return 0;
}

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
static int _zNLESolveNR(zNLE *nle, zVec var, void *util, double tol_step, double tol_err, int iter, double *err)
{
  int i, ret = -1;
  double a;
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
  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    nle->_jac( nle, var, nle->_j, util );
    /* uses SR-inverse instead of MP-inverse, similar to LM method. */
    zLESolveSRDST( nle->_j, nle->_f, nle->wn, nle->we, nle->_opt._d, &workspace );
    if( ( a = _zNLEStep( nle, var, util ) ) == 0 ) goto TERMINATE;
    zVecCopyNC( nle->_fp, nle->_f );
    if( zVecIsTol( nle->_f, tol_err ) ){
      ret = i + 1; /* succeed */
      goto TERMINATE;
    }
    if( _zNLEStepIsSmall( nle, var, a, tol_step ) ) goto TERMINATE; /* stagnate */
  }
  ZITERWARN( iter );
 TERMINATE:
  _zNLEEvalRes( nle, var, util, err );
  return ret;
}

/* solve simultaneous nonlinear equations by Broyden's method. */
static int _zNLESolveBroyden(zNLE *nle, zVec var, void *util, double tol_step, double tol_err, int iter, double *err)
{
  int i, ret = -1;
  double a;

  nle->f( var, nle->_f, util );
  if( zVecIsTol( nle->_f, tol_err ) ){
    ret = 0;
    goto TERMINATE;
  }
  nle->_jac( nle, var, nle->_j, util );
  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    if( !zLESolveGauss( nle->_j, nle->_f, nle->_opt._d ) ){
      nle->_jac( nle, var, nle->_j, util );
      if( !zLESolveGauss( nle->_j, nle->_f, nle->_opt._d ) ){
        ZRUNWARN( ZM_WARN_OPT_BADINI ); /* dummy */
        goto TERMINATE;
      }
    }
    if( ( a = _zNLEStep( nle, var, util ) ) == 0 ) goto TERMINATE;
    zVecCopyNC( nle->_fp, nle->_f );
    if( zVecIsTol( nle->_f, tol_err ) ){
      ret = i + 1; /* succeed */
      goto TERMINATE;
    }
    if( _zNLEStepIsSmall( nle, var, a, tol_step ) ) goto TERMINATE; /* stagnate */
    /* Broyden's update rule */
    zMatCatDyadNC( nle->_j, -1.0/zVecSqrNorm(nle->_opt._d), nle->_f, nle->_opt._d );
  }
  ZITERWARN( iter );
 TERMINATE:
  _zNLEEvalRes( nle, var, util, err );
  return ret;
}

#define ZM_NLE_WN_DEFAULT ( 1.0e-4 )
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
  zVecSetAll( nle->wn, ZM_NLE_WN_DEFAULT );
  zVecSetAll( nle->we, 1.0 );
  return zNLEAssignNR( nle );
}

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
