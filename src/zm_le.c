/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le - linear equation.
 */

#include <zm/zm_le.h>

/* zBalancingColDST
 * - directly make a matrix column-balanced.
 */
void zBalancingColDST(zMat m, zVec s)
{
  register int i, j;
  double tmp;

  for( i=0; i<zMatColSizeNC(m); i++ ){
    zVecSetElem( s, i, fabs( zMatElem(m,0,i) ) );
    for( j=1; j<zMatRowSizeNC(m); j++ ){
      tmp = fabs( zMatElem(m,j,i) );
      if( tmp > zVecElem(s,i) ) zVecSetElem( s, i, tmp );
    }
    if( zVecElem(s,i) == 0 ) continue;
    /* inverse column-balancing factor */
    zVecSetElem( s, i, 1.0 / zVecElem(s,i) );
    for( j=0; j<zMatRowSizeNC(m); j++ )
      zMatElem(m,j,i) *= zVecElem(s,i);
  }
}

/* zBalancingDST
 * - directly make a pair of matrix and vector balanced.
 */
void zBalancingDST(zMat m, zVec v, zVec s)
{
  register int i;
  double *mp, max;

  if( s ) /* if necessary for column-balancing */
    zBalancingColDST( m, s );
  for( mp=zMatBuf(m), i=0; i<zMatRowSizeNC(m); mp+=zMatColSizeNC(m), i++ ){
    if( ( max = zDataAbsMax( mp, zMatColSizeNC(m), NULL ) ) == 0 )
      continue;
    zRawVecDivDRC( mp, max, zMatColSizeNC(m) );
    zVecElem(v,i) /= max;
  }
}

/* zBalancing
 * - make a pair of matrix and vector balanced.
 */
bool zBalancing(zMat morg, zVec vorg, zMat m, zVec v, zVec s)
{
  if( !zMatSizeIsEqual( morg, m ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return false;
  }
  if( !zVecSizeIsEqual( vorg, v ) ||
      ( s && !zVecSizeIsEqual( vorg, s ) ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );
    return false;
  }
  if( zMatRowSizeNC(m) != zVecSizeNC(v) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return false;
  }
  zMatCopyNC( morg, m );
  zVecCopyNC( vorg, v );
  zBalancingDST( m, v, s );
  return true;
}

/* zLEResidual
 * - residual b - a x.
 */
zVec zLEResidual(zMat a, zVec b, zVec x, zVec res)
{
  zMulMatVecNC( a, x, res );
  return zVecSubNC( b, res, res );
}

/* zLESolveGaussDST
 * - linear equation solver based on Gaussian elimination method
 *   (destructive).
 */
zVec zLESolveGaussDST(zMat a, zVec b, zVec ans, zIndex idx, zVec s)
{
  register int i, j, k;
  int n, p, q;
  double ahead;
  double x;

  n = zVecSizeNC( b );
  zBalancingDST( a, b, s );
  /* forward elimination */
  for( i=0; i<n; i++ ){
    p = zPivoting( a, idx, i, i );
    if( ( ahead = zMatElem(a,p,i) ) == 0 ){
      ZRUNERROR( ZM_ERR_LE_SINGULAR );
      return NULL;
    }
    ahead = 1.0 / ahead;
    zMatSetElem( a, p, i, 1 );
    for( j=i+1; j<n; j++ )
      zMatElem(a,p,j) *= ahead;
    zVecElem(b,p) *= ahead;
    for( j=i+1; j<n; j++ ){
      q = zIndexElem( idx, j );
      if( !zIsTiny( ahead = zMatElem(a,q,i) ) ){
        for( k=i+1; k<n; k++ )
          zMatElem(a,q,k) -= zMatElem(a,p,k) * ahead;
        zVecElem(b,q) -= zVecElem(b,p) * ahead;
      }
      zMatSetElem( a, q, i, 0 );
    }
  }
  /* backward elimination */
  for( i=n-1; i>=0; i-- ){
    p = zIndexElem( idx, i );
    x = zVecElem( b, p );
    for( j=n-1; j>i; j-- )
      x -= zMatElem(a,p,j)*zVecElem(ans,j);
    zVecSetElem( ans, i, x );
  }
  if( s ) zVecAmpDRC( ans, s );
  return ans;
}

/* zLESolveGauss
 * - linear equation solver based on Gaussian elimination method.
 */
zVec zLESolveGauss(zMat a, zVec b, zVec ans)
{
  int n;
  zMat acp;
  zVec bcp, s;
  zIndex idx;

  n = zVecSizeNC(b);
  if( !zMatIsSqr(a) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return NULL;
  }
  if( zMatColSize(a) != n || zVecSize(ans) != n ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  acp = zMatClone( a );
  bcp = zVecClone( b );
  s = zVecAlloc( n );
  idx = zIndexCreate( n );
  if( acp && bcp && idx && s )
    ans = zLESolveGaussDST( acp, bcp, ans, idx, s );
  else
    ZALLOCERROR();
  zMatFree( acp );
  zVecFree( bcp );
  zVecFree( s );
  zIndexFree( idx );
  return ans;
}

/* zLESolve_L
 * - inner solver of zLUSolve for L matrix.
 */
zVec zLESolve_L(zMat lmat, zVec b, zVec ans, zIndex idx)
{
  register int i, j;
  int p;
  double x;

  for( i=0; i<zArrayNum(idx); i++ ){
    x = zVecElem( b, (p=zIndexElem(idx,i)) );
    for( j=0; j<i; j++ )
      x -= zMatElem(lmat,p,j)*zVecElem(ans,j);
    zVecSetElem( ans, i, x/zMatElem(lmat,p,i) );
  }
  return ans;
}

/* zLESolve_U
 * - inner solver of zLUSolve for U matrix.
 */
zVec zLESolve_U(zMat umat, zVec b, zVec ans)
{
  register int i, j;
  double x;

  for( i=zVecSizeNC(b)-1; i>=0; i-- ){
    x = zVecElem( b, i );
    for( j=zVecSizeNC(b)-1; j>i; j-- )
      x -= zMatElem(umat,i,j)*zVecElem(ans,j);
    zVecSetElem( ans, i, x );
  }
  return ans;
}

/* zLESolve_L_U
 * - inner solver of zLESolveLU.
 */
zVec zLESolve_L_U(zMat lmat, zMat umat, zVec b, zVec ans, zIndex idx)
{
  zVec c;

  if( !zMatIsSqr(lmat) || !zMatIsSqr(umat) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return NULL;
  }
  if( zMatRowSize(lmat) != zArrayNum(idx) ||
      zMatRowSize(umat) != zArrayNum(idx) ||
      zVecSize(b) != zArrayNum(idx) ||
      zVecSize(ans) != zArrayNum(idx) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  if( !( c = zVecAlloc( zArrayNum(idx) ) ) ) return NULL;
  zLESolve_L( lmat, b, c, idx );
  zLESolve_U( umat, c, ans );
  zVecFree( c );
  return ans;
}

/* zLESolveLU
 * - linear equation solver based on LU decomposition method.
 */
zVec zLESolveLU(zMat a, zVec b, zVec ans)
{
  int n;
  zMat lmat, umat;
  zIndex idx;

  n = zVecSizeNC( b );
  lmat = zMatAllocSqr( n );
  umat = zMatAllocSqr( n );
  idx = zIndexCreate( n );
  if( !lmat || !umat || !idx ) goto TERMINATE;

  if( zLUDecomp( a, lmat, umat, idx ) < zMatRowSizeNC(a) ){
    ZRUNERROR( ZM_ERR_LE_SINGULAR );
    ans = NULL;
    goto TERMINATE;
  }
  zLESolve_L_U( lmat, umat, b, ans, idx );

 TERMINATE:
  zMatFree( lmat );
  zMatFree( umat );
  zIndexFree( idx );
  return ans;
}

/* zLESolveRI
 * - linear equation solver: Residual iteration on LU decomposition.
 */
zVec zLESolveRI(zMat a, zVec b, zVec ans)
{
  register int i;
  int n;
  zMat lmat, umat;
  zVec res, err;
  double err_norm, err_norm_old = HUGE_VAL;
  zIndex idx;

  n = zVecSizeNC( b );
  lmat = zMatAllocSqr( n );
  umat = zMatAllocSqr( n );
  res = zVecAlloc( n );
  err = zVecAlloc( n );
  idx = zIndexCreate( n );
  if( !lmat || !umat || !res || !err || !idx ) goto TERMINATE;

  if( zLUDecomp( a, lmat, umat, idx ) < zMatRowSizeNC(a) ){
    ZRUNERROR( ZM_ERR_LE_SINGULAR );
    ans = NULL;
    goto TERMINATE;
  }
  zLESolve_L_U( lmat, umat, b, ans, idx );
  for( i=0; i<Z_MAX_ITER_NUM; i++ ){
    zLEResidual( a, b, ans, res );
    err_norm = zVecNorm( res );
    if( err_norm >= err_norm_old ) goto TERMINATE;
    err_norm_old = err_norm;
    zLESolve_L_U( lmat, umat, res, err, idx );
    zVecAddNCDRC( ans, err );
  }
  ZITERWARN( Z_MAX_ITER_NUM );

 TERMINATE:
  zMatFree( lmat );
  zMatFree( umat );
  zVecFree( res );
  zVecFree( err );
  zIndexFree( idx );
  return ans;
}

/* zLESolveGS
 * - linear equation solver: Gauss-Seidel's method.
 */
zVec zLESolveGS(zMat a, zVec b, zVec ans)
{
  register int i, j, k;
  int p, count;
  double x;
  zIndex idx;

  if( !zMatIsSqr(a) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return NULL;
  }
  if( !zMatColVecSizeIsEqual(a,ans) ||
      !zMatRowVecSizeIsEqual(a,b) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  if( !( idx = zIndexCreate(zVecSizeNC(ans)) ) ) return NULL;
  for( i=0; i<zArrayNum(idx); i++ )
    zPivoting( a, idx, i, i );

  for( i=0; ; i++ ){
    for( count=0, j=0; j<zVecSizeNC(b); j++ ){
      p = zIndexElem(idx,j);
      x = zVecElem(b,p);
      for( k=0; k<zVecSizeNC(ans); k++ )
        if( k != j ) x -= zMatElem(a,p,k)*zVecElem(ans,k);
      x /= zMatElem(a,p,j);
      if( zIsTiny( x - zVecElem(ans,j) ) ) count++;
      zVecSetElem( ans, j, x );
    }
    if( count == zVecSizeNC(ans) ) break;
    if( i == Z_MAX_ITER_NUM ){
      ZITERWARN( Z_MAX_ITER_NUM );
      break;
    }
  }
  zIndexFree( idx );
  return ans;
}
