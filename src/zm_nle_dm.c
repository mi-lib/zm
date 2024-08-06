/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_nle_dm - nonlinear equation:
 * simultaneous nonlinear equation solver based on descent method.
 */

#include <zm/zm_nle.h>

static zMat _zNLEJacobi(zNLE *nle, zVec var, zMat j, void *util){
  return nle->jac( var, j, util );
}
static zMat _zNLEJacobiNG(zNLE *nle, zVec var, zMat j, void *util)
{
  int i;
  double org;

  for( i=0; i<zVecSizeNC(var); i++ ){
    org = zVecElem(var,i);
    zVecSetElemNC( var, i, org+Z_OPT_EPS );
    nle->f( var, nle->_adg, util );
    zVecSetElemNC( var, i, org-Z_OPT_EPS );
    nle->f( var, nle->_prg, util );
    zRawVecSub( zVecBuf(nle->_adg), zVecBuf(nle->_prg),
      &zMatBuf(j)[zMatColSizeNC(j)*i], zMatRowSizeNC(j) );
    zVecSetElemNC( var, i, org );
  }
  zMatMulDRC( j, 0.5/Z_OPT_EPS );
  return j;
}

#define _zm_nle(u) ( (zNLE *)u )

static double _zNLEEval(zVec var, void *util)
{
  _zm_nle(util)->f( var, _zm_nle(util)->_f, _zm_nle(util)->util );
  zVecAmpNC( _zm_nle(util)->_f, _zm_nle(util)->we, _zm_nle(util)->_fw );
  return 0.5 * zVecInnerProd( _zm_nle(util)->_f, _zm_nle(util)->_fw );
}

static zVec _zNLEGrad(zVec var, zVec grad, void *util)
{
  _zm_nle(util)->_jac( _zm_nle(util), var, _zm_nle(util)->_j, _zm_nle(util)->util );
  return zMulMatTVecNC( _zm_nle(util)->_j, _zm_nle(util)->_fw, grad );
}

static zMat _zNLEHess(zVec var, zMat h, void *util)
{
  return zMatTQuadNC( _zm_nle(util)->_j, _zm_nle(util)->we, h );
}

#undef _zm_nle

/* update variable vector toward descent direction. */
static zVec _zNLEStep(zNLE *nle, zVec var, void *util)
{
  double a;

  for( a=1.0; a>zTOL; a*=0.5 ){
    zVecCatNC( var, -a, nle->_opt._d, nle->_opt._x );
    nle->f( nle->_opt._x, nle->_fp, util );
    if( !zVecIsNan( nle->_fp ) )
      /* now, f( var ) is stored in _q. if need it, copy it to _f. */
      return zVecCopyNC( nle->_opt._x, var );
  }
  ZRUNWARN( ZM_WARN_OPT_BADINI ); /* dummy */
  return NULL;
}

static void _zNLEEvalRes(zNLE *nle, zVec var, void *util, double *err)
{
  if( err ){
    nle->util = util;
    *err = _zNLEEval( var, nle );
  }
}

static int _zNLESolveDM(zNLE *nle, zVec var, void *util, double tol, int iter, double *err)
{
  nle->util = util;
  return zOptDMSolve( &nle->_opt, var, nle, tol, iter, err );
}

/* Newton-Raphson's method. */
static int _zNLESolveNR(zNLE *nle, zVec var, void *util, double tol, int iter, double *err)
{
  int i, ret = -1;
  zLEWorkspace workspace;

  workspace.m = nle->_opt._h;
  workspace.v1 = nle->_opt._p;
  workspace.idx1 = nle->_opt._idx;
  workspace.s = nle->_opt._q;
  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    nle->f( var, nle->_f, util );
    if( zVecIsTol( nle->_f, tol ) ){
      _zNLEEvalRes( nle, var, util, err );
      ret = i;
      goto TERMINATE; /* succeed */
    }
    nle->_jac( nle, var, nle->_j, util );
    /* uses SR-inverse instead of MP-inverse, similar to LM method. */
    zLESolveSRDST( nle->_j, nle->_f, nle->wn, nle->we, nle->_opt._d, &workspace );
    if( !_zNLEStep( nle, var, util ) ) break;
  }
  ZITERWARN( iter );
 TERMINATE:
eprintf("eheraehera.\n");
  zLEWorkspaceFree( &workspace );
eprintf("ungyaasu.\n");
  return ret;
}

/* solve simultaneous nonlinear equations by Broyden's method. */
static int _zNLESolveBroyden(zNLE *nle, zVec var, void *util, double tol, int iter, double *err)
{
  int i;

  nle->f( var, nle->_f, util );
  if( zVecIsTol( nle->_f, tol ) ){
    _zNLEEvalRes( nle, var, util, err );
    return 0;
  }
  nle->_jac( nle, var, nle->_j, util );
  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    if( !zLESolveGauss( nle->_j, nle->_f, nle->_opt._d ) ){
      nle->_jac( nle, var, nle->_j, util );
      if( !zLESolveGauss( nle->_j, nle->_f, nle->_opt._d ) ){
        ZRUNWARN( ZM_WARN_OPT_BADINI ); /* dummy */
        break;
      }
    }
    if( !_zNLEStep( nle, var, util ) ) break;
    zVecCopyNC( nle->_fp, nle->_f );
    if( zVecIsTol( nle->_f, tol ) ){
      _zNLEEvalRes( nle, var, util, err );
      return i; /* succeeded. */
    }
    /* Broyden's update rule */
    zMatCatDyadNC( nle->_j, -1.0/zVecSqrNorm(nle->_opt._d), nle->_f, nle->_opt._d );
  }
  ZITERWARN( iter );
  return -1;
}

#define Z_NLE_WN_DEFAULT ( 1.0e-4 )
zNLE *zNLECreate(zNLE *nle, int nv, int ne, double scale, zVec (*f)(zVec,zVec,void*), zMat (*jac)(zVec,zMat,void*))
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
      !zOptDMCreate( &nle->_opt, nv, scale, _zNLEEval, _zNLEGrad, _zNLEHess ) ){
    zNLEDestroy( nle );
    return NULL;
  }
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
  zVecSetAll( nle->wn, Z_NLE_WN_DEFAULT );
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

zNLE *zNLEAssignSD(zNLE *nle, const char *stepmethod){
  if( !zOptDMAssignSD( &nle->_opt, stepmethod ) )
    return NULL;
  nle->_solve = _zNLESolveDM;
  return nle;
}
zNLE *zNLEAssignLM(zNLE *nle, const char *stepmethod){
  if( !zOptDMAssignLM( &nle->_opt, stepmethod ) )
    return NULL;
  nle->_solve = _zNLESolveDM;
  return nle;
}
zNLE *zNLEAssignVM(zNLE *nle, const char *stepmethod, const char *updatemethod){
  if( !zOptDMAssignVM( &nle->_opt, stepmethod, updatemethod ) )
    return NULL;
  nle->_solve = _zNLESolveDM;
  return nle;
}
zNLE *zNLEAssignCG(zNLE *nle){
  if( !zOptDMAssignCG( &nle->_opt ) )
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
