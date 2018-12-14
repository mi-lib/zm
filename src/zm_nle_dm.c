/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_nle_dm - nonlinear equation:
 * simultaneous nonlinear equation solver based on descent method.
 */

#include <zm/zm_nle.h>

static zMat _zNLEJacobi(zNLE *nle, zVec var, zMat j, void *util);
static zMat _zNLEJacobiNG(zNLE *nle, zVec var, zMat j, void *util);

static double _zNLEEval(zVec var, void *util);
static zVec _zNLEGrad(zVec var, zVec grad, void *util);
static zMat _zNLEHess(zVec var, zMat h, void *util);

static zVec _zNLEStep(zNLE *nle, zVec var, void *util);
static void _zNLEEvalRes(zNLE *nle, zVec var, void *util, double *err);
static int _zNLESolveDM(zNLE *nle, zVec var, void *util, double tol, int iter, double *err);
static int _zNLESolveNR(zNLE *nle, zVec var, void *util, double tol, int iter, double *err);
static int _zNLESolveBroyden(zNLE *nle, zVec var, void *util, double tol, int iter, double *err);

zMat _zNLEJacobi(zNLE *nle, zVec var, zMat j, void *util){
  return nle->jac( var, j, util );
}
zMat _zNLEJacobiNG(zNLE *nle, zVec var, zMat j, void *util)
{
  register int i;
  double org;

  for( i=0; i<zVecSizeNC(var); i++ ){
    org = zVecElem(var,i);
    zVecSetElem( var, i, org+Z_OPT_EPS );
    nle->f( var, nle->_adg, util );
    zVecSetElem( var, i, org-Z_OPT_EPS );
    nle->f( var, nle->_prg, util );
    zRawVecSub( zVecBuf(nle->_adg), zVecBuf(nle->_prg),
      &zMatBuf(j)[zMatColSizeNC(j)*i], zMatRowSizeNC(j) );
    zVecSetElem( var, i, org );
  }
  zMatMulDRC( j, 0.5/Z_OPT_EPS );
  return j;
}

double _zNLEEval(zVec var, void *util)
{
  zNLE *nle;

  nle = util;
  nle->f( var, nle->_f, nle->util );
  zVecAmpNC( nle->_f, nle->we, nle->_fw );
  return 0.5 * zVecInnerProd( nle->_f, nle->_fw );
}

zVec _zNLEGrad(zVec var, zVec grad, void *util)
{
  zNLE *nle;

  nle = util;
  nle->_jac( nle, var, nle->_j, nle->util );
  return zMulMatTVecNC( nle->_j, nle->_fw, grad );
}

zMat _zNLEHess(zVec var, zMat h, void *util)
{
  zNLE *nle;

  nle = util;
  return zMatTQuadNC( nle->_j, nle->we, h );
}

/* (static)
 * _zNLEStep
 * - update variable vector toward descent direction.
 */
zVec _zNLEStep(zNLE *nle, zVec var, void *util)
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

void _zNLEEvalRes(zNLE *nle, zVec var, void *util, double *err)
{
  if( err ){
    nle->util = util;
    *err = _zNLEEval( var, nle );
  }
}

int _zNLESolveDM(zNLE *nle, zVec var, void *util, double tol, int iter, double *err)
{
  nle->util = util;
  return zOptDMSolve( &nle->_opt, var, nle, tol, iter, err );
}

/* zNLESolveNR
 * - Newton-Raphson's method.
 */
int _zNLESolveNR(zNLE *nle, zVec var, void *util, double tol, int iter, double *err)
{
  register int i;

  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    nle->f( var, nle->_f, util );
    if( zVecIsTol( nle->_f, tol ) ){
      _zNLEEvalRes( nle, var, util, err );
      return i; /* succeed */
    }
    nle->_jac( nle, var, nle->_j, util );
    /* uses SR-inverse instead of MP-inverse,
       rather similar to LM method. */
    zLESolveSRDST( nle->_j, nle->_f, nle->wn, nle->we, nle->_opt._d,
      nle->_opt._h, nle->_opt._p, nle->_opt._idx, nle->_opt._q );
    if( !_zNLEStep( nle, var, util ) ) break;
  }
  ZITERWARN( iter );
  return -1;
}

/* zNLESolveBroyden
 * - solve simultaneous nonlinear equations by Broyden's method.
 */
int _zNLESolveBroyden(zNLE *nle, zVec var, void *util, double tol, int iter, double *err)
{
  register int i;

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
  zVecFreeAO( 7,
    nle->wn, nle->we, nle->_f, nle->_fw, nle->_fp, nle->_adg, nle->_prg );
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
