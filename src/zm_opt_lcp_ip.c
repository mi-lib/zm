/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_lcp_ip - optimization tools:
 * Potra's predictor-corrector algorithm on
 * infeasible-interior-point method.
 */

#include <zm/zm_opt.h>

/* working memory for Potra's predictor-corrector algorithm */
typedef struct{
  zMat g;       /* Newton's gradient */
  zVec r;       /* residual vector */
  zVec dz, dw;  /* gradient vectors */
  zVec b, c, s; /* working space */
  zIndex idx;
  double e;     /* complementarity error */
  double t;     /* approximate analytical center point */
} _zLCPIP_PC;

#define _Z_LCPIP_PC_A 0.25
#define _Z_LCPIP_PC_B 0.50
#define _Z_LCPIP_PC_R 100.0

static bool _zLCPIP_PCAlloc(_zLCPIP_PC *wm, int n);
static void _zLCPIP_PCFree(_zLCPIP_PC *wm);
static bool _zLCPIP_PCIni(_zLCPIP_PC *wm, zMat m, zVec q, zVec w, zVec z);
static void _zLCPIP_PCErr(_zLCPIP_PC *wm, zMat m, zVec q, zVec w, zVec z);
static bool _zLCPIP_PCGrad(_zLCPIP_PC *wm, zMat m, zVec w, zVec z);
static bool _zLCPIP_PCPred(_zLCPIP_PC *wm, zMat m, zVec w, zVec z);
static bool _zLCPIP_PCCorr(_zLCPIP_PC *wm, zMat m, zVec w, zVec z);
static bool _zLCPIP_PCStep(_zLCPIP_PC *wm, zVec w, zVec z, double *step);

/* (static)
 * _zLCPIP_PCAlloc
 * - allocate internal working memory for LCP-IP-PC.
 */
bool _zLCPIP_PCAlloc(_zLCPIP_PC *wm, int n)
{
  wm->g = zMatAllocSqr( n );
  wm->r = zVecAlloc( n );
  wm->dz = zVecAlloc( n );
  wm->dw = zVecAlloc( n );
  wm->b = zVecAlloc( n );
  wm->c = zVecAlloc( n );
  wm->s = zVecAlloc( n );
  wm->idx = zIndexCreate( n );
  if( !wm->g || !wm->r || !wm->dz || !wm->dw ||
      !wm->b || !wm->c || !wm->s || !wm->idx ){
    ZALLOCERROR();
    _zLCPIP_PCFree( wm );
    return false;
  }
  return true;
}

/* (static)
 * _zLCPIP_PCFree
 * - free internal working memory for LCP-IP-PC.
 */
void _zLCPIP_PCFree(_zLCPIP_PC *wm)
{
  zMatFree( wm->g );
  zVecFreeAO( 6, wm->r, wm->dz, wm->dw, wm->b, wm->c, wm->s );
  zIndexFree( wm->idx );
}

/* (static)
 * _zLCPIP_PCIni
 * - heuristicly initialize vectors for LCP-IP-PC.
 */
bool _zLCPIP_PCIni(_zLCPIP_PC *wm, zMat m, zVec q, zVec w, zVec z)
{
  zVecSetAll( z, _Z_LCPIP_PC_R / zVecSizeNC(z) );
  zVecSetAll( w, _Z_LCPIP_PC_R / zVecSizeNC(w) );
  /* approximate analytical center point */
  wm->t = zVecInnerProd( w, z ) / ( zVecSizeNC(z) + _Z_LCPIP_PC_A * sqrt(zVecSizeNC(z)) );
  return true;
}

/* (static)
 * _zLCPIP_PCErr
 * - compute residual vector and complementarity gap for LCP-IP-PC.
 */
void _zLCPIP_PCErr(_zLCPIP_PC *wm, zMat m, zVec q, zVec w, zVec z)
{
  /* residual vector */
  zMulMatVecNC( m, z, wm->r );
  zVecAddNCDRC( wm->r, q );
  zVecSubNCDRC( wm->r, w );
  /* complementarity gap */
  wm->e = zVecInnerProd( w, z );
}

/* (static)
 * _zLCPIP_PCGrad
 * - gradient matrix and its inverse for LCP-IP-PC.
 */
bool _zLCPIP_PCGrad(_zLCPIP_PC *wm, zMat m, zVec w, zVec z)
{
  register int i;

  for( i=0; i<zVecSizeNC(w); i++ ){
    zRawVecMul( zMatRowBuf(m,i), zVecElem(z,i), zMatRowBuf(wm->g,i), zMatColSizeNC(m) );
    zMatElem(wm->g,i,i) += zVecElem(w,i);
  }
  return zLESolveGaussDST( wm->g, wm->b, wm->dz, wm->idx, wm->s ) ? true : false;
}

/* (static)
 * _zLCPIP_PCPred
 * - predictor direction vector for LCP-IP-PC.
 */
bool _zLCPIP_PCPred(_zLCPIP_PC *wm, zMat m, zVec w, zVec z)
{
  zVecAddNC( w, wm->r, wm->b );
  zVecAmpNCDRC( wm->b, z );
  zVecRevDRC( wm->b );
  if( !_zLCPIP_PCGrad( wm, m, w, z ) ) return false;
  zMulMatVecNC( m, wm->dz, wm->dw );
  zVecAddNCDRC( wm->dw, wm->r );
  return true;
}

/* (static)
 * _zLCPIP_PCCorr
 * - corrector (centering) direction vector for LCP-IP-PC.
 */
bool _zLCPIP_PCCorr(_zLCPIP_PC *wm, zMat m, zVec w, zVec z)
{
  register int i;

  zVecAmpNC( w, z, wm->b );
  for( i=0; i<zVecSizeNC(wm->b); i++ )
    zVecSetElem( wm->b, i, wm->t - zVecElem(wm->b,i) );
  if( !_zLCPIP_PCGrad( wm, m, w, z ) ) return false;
  zMulMatVecNC( m, wm->dz, wm->dw );
  return true;
}

/* (static)
 * _zLCPIP_PCStep
 * - updating step of LCP-IP-PC.
 */
bool _zLCPIP_PCStep(_zLCPIP_PC *wm, zVec w, zVec z, double *step)
{
  double k1, k2, a;
  register int i;

  zVecAmpNC( w, z, wm->b );
  for( i=0; i<zVecSizeNC(wm->b); i++ )
    zVecElem(wm->b,i) -= wm->t;
  zVecAmpNC( wm->dw, wm->dz, wm->c );
  k2 = zVecInnerProd( wm->b, wm->c );
  k1 = ( zSqr(_Z_LCPIP_PC_B*wm->t) - zVecSqrNorm(wm->b) ) / k2;
  a = k1 / ( 1 + sqrt( 1 + k1/k2*zVecSqrNorm(wm->c) ) );
  if( a < 0 && a > -4.0 ){
    ZRUNERROR( ZM_ERR_OPT_STEP );
    return false;
  }
  *step = 2.0 / ( 1 + sqrt( 1 + 4.0/a ) );
  wm->t *= 1 - *step;
  return true;
}

/* zLCPSolveIP
 * - solve linear complementarity problem by Potra's
 *   predictor-corrector algorithm on infeasible-interior-point
 *   method.
 */
bool zLCPSolveIP(zMat m, zVec q, zVec w, zVec z)
{
  _zLCPIP_PC wm;
  zVec worg;
  double a;
  int iter = 0;
  bool ret = false;
  register int i;

  if( !zMatIsSqr(m) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return false;
  }
  if( !zMatRowVecSizeIsEqual(m,q) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return false;
  }
  if( ( w && ( !zVecSizeIsEqual(w,q) || !zVecSizeIsEqual(w,z) ) ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );
    return false;
  }
  if( !_zLCPIP_PCAlloc( &wm, zVecSizeNC(q) ) )
    return false;
  if( !( worg = w ) && !( w = zVecAlloc( zVecSizeNC(q) ) ) ){
    ZALLOCERROR();
    goto TERMINATE;
  }
  /* initial interior point */
  _zLCPIP_PCIni( &wm, m, q, w, z );
  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    _zLCPIP_PCErr( &wm, m, q, w, z );
    if( zVecIsTiny( wm.r ) && zIsTiny( wm.e ) ) break;
    if( !_zLCPIP_PCPred( &wm, m, w, z ) ||
        !_zLCPIP_PCStep( &wm, w, z, &a ) ){
      ZRUNWARN( ZM_ERR_OPT_UNSOLVE );
      goto TERMINATE;
    }
    zVecCatNCDRC( w, a, wm.dw );
    zVecCatNCDRC( z, a, wm.dz );
    if( !_zLCPIP_PCCorr( &wm, m, w, z ) ){
      ZRUNWARN( ZM_ERR_OPT_UNSOLVE );
      goto TERMINATE;
    }
    zVecAddNCDRC( w, wm.dw );
    zVecAddNCDRC( z, wm.dz );
  }
  zVecTouchup( z );
  if( i == iter ) ZITERWARN( iter );
  ret = true;

 TERMINATE:
  if( !worg ) zVecFree( w );
  _zLCPIP_PCFree( &wm );
  return ret;
}
