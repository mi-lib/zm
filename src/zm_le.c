/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le - linear equation.
 */

#include <zm/zm_le.h>

/* directly make a matrix column-balanced. */
void zMatColBalancingDST(zMat m, zVec s)
{
  int i, j;
  double tmp;

  for( i=0; i<zMatColSizeNC(m); i++ ){
    zVecSetElemNC( s, i, fabs( zMatElemNC(m,0,i) ) );
    for( j=1; j<zMatRowSizeNC(m); j++ ){
      tmp = fabs( zMatElemNC(m,j,i) );
      if( tmp > zVecElemNC(s,i) ) zVecSetElemNC( s, i, tmp );
    }
    if( zVecElemNC(s,i) == 0 ) continue;
    /* inverse column-balancing factor */
    zVecSetElemNC( s, i, 1.0 / zVecElemNC(s,i) );
    for( j=0; j<zMatRowSizeNC(m); j++ )
      zMatElemNC(m,j,i) *= zVecElemNC(s,i);
  }
}

/* directly make a pair of matrix and vector row-balanced. */
void zMatVecRowBalancingDST(zMat m, zVec v)
{
  int i;
  double *mp, tmp;

  for( mp=zMatBuf(m), i=0; i<zMatRowSizeNC(m); mp+=zMatColSizeNC(m), i++ ){
    if( ( tmp = zDataAbsMax( mp, zMatColSizeNC(m), NULL ) ) == 0 )
      continue;
    zRawVecDivDRC( mp, tmp, zMatColSizeNC(m) );
    zVecElemNC(v,i) /= tmp;
  }
}

/* directly make a pair of matrix and vector balanced. */
void zMatVecBalancingDST(zMat m, zVec v, zVec s)
{
  if( s )
    zMatColBalancingDST( m, s );
  zMatVecRowBalancingDST( m, v );
}

/* make a pair of matrix and vector balanced. */
bool zMatVecBalancing(const zMat morg, const zVec vorg, zMat m, zVec v, zVec s)
{
  if( !zMatSizeEqual( morg, m ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return false;
  }
  if( !zVecSizeEqual( vorg, v ) ||
      ( s && !zVecSizeEqual( vorg, s ) ) ){
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );
    return false;
  }
  if( zMatRowSizeNC(m) != zVecSizeNC(v) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return false;
  }
  zMatCopyNC( morg, m );
  zVecCopyNC( vorg, v );
  zMatVecBalancingDST( m, v, s );
  return true;
}

/* directly make a pair of matrices row-balanced. */
void zMatMatRowBalancingDST(zMat m1, zMat m2)
{
  int i;
  double *mp1, *mp2, tmp;

  for( mp1=zMatBuf(m1), mp2=zMatBuf(m2), i=0; i<zMatRowSizeNC(m1); mp1+=zMatColSizeNC(m1), mp2+=zMatColSizeNC(m2), i++ ){
    if( ( tmp = zDataAbsMax( mp1, zMatColSizeNC(m1), NULL ) ) == 0 )
      continue;
    zRawVecDivDRC( mp1, tmp, zMatColSizeNC(m1) );
    zRawVecDivDRC( mp2, tmp, zMatColSizeNC(m2) );
  }
}

/* directly make a pair of matrices balanced. */
void zMatMatBalancingDST(zMat m1, zMat m2, zVec s)
{
  if( s )
    zMatColBalancingDST( m1, s );
  zMatMatRowBalancingDST( m1, m2 );
}

/* directly make a pair of matrices and a vector row-balanced. */
void zMatMatVecRowBalancingDST(zMat m1, zMat m2, zVec v)
{
  int i;
  double *mp1, *mp2, tmp;

  for( mp1=zMatBuf(m1), mp2=zMatBuf(m2), i=0; i<zMatRowSizeNC(m1); mp1+=zMatColSizeNC(m1), mp2+=zMatColSizeNC(m2), i++ ){
    if( ( tmp = zDataAbsMax( mp1, zMatColSizeNC(m1), NULL ) ) == 0 )
      continue;
    zRawVecDivDRC( mp1, tmp, zMatColSizeNC(m1) );
    zRawVecDivDRC( mp2, tmp, zMatColSizeNC(m2) );
    zVecElemNC(v,i) /= tmp;
  }
}

/* directly make a pair of matrices and a vector balanced. */
void zMatMatVecBalancingDST(zMat m1, zMat m2, zVec v, zVec s)
{
  if( s )
    zMatColBalancingDST( m1, s );
  zMatMatVecRowBalancingDST( m1, m2, v );
}

/* residual b - a x. */
zVec zLEResidual(const zMat a, const zVec b, const zVec x, zVec res)
{
  zMulMatVecNC( a, x, res );
  return zVecSubNC( b, res, res );
}

/* linear equation solver based on Gauss's elimination method (destructive). */
zVec zLESolveGaussDST(zMat a, zVec b, zVec ans, zIndex idx, zVec s)
{
  int i, j, k, n, p, q;
  double ahead;
  double x;

  n = zVecSizeNC( b );
  zMatVecBalancingDST( a, b, s );
  /* forward elimination */
  for( i=0; i<n; i++ ){
    p = zMatPivoting( a, idx, i, i );
    if( ( ahead = zMatElemNC(a,p,i) ) == 0 ){
      ZRUNERROR( ZM_ERR_MAT_SINGULAR );
      return NULL;
    }
    ahead = 1.0 / ahead;
    zMatSetElemNC( a, p, i, 1 );
    for( j=i+1; j<n; j++ )
      zMatElemNC(a,p,j) *= ahead;
    zVecElemNC(b,p) *= ahead;
    for( j=i+1; j<n; j++ ){
      q = zIndexElemNC( idx, j );
      if( !zIsTiny( ahead = zMatElemNC(a,q,i) ) ){
        for( k=i+1; k<n; k++ )
          zMatElemNC(a,q,k) -= zMatElemNC(a,p,k) * ahead;
        zVecElemNC(b,q) -= zVecElemNC(b,p) * ahead;
      }
      zMatSetElemNC( a, q, i, 0 );
    }
  }
  /* backward elimination */
  for( i=n-1; i>=0; i-- ){
    x = zVecElemNC( b, ( p = zIndexElemNC(idx,i) ) );
    for( j=n-1; j>i; j-- )
      x -= zMatElemNC(a,p,j)*zVecElemNC(ans,j);
    zVecSetElemNC( ans, i, x );
  }
  if( s ) zVecAmpDRC( ans, s );
  return ans;
}

/* linear equation solver based on Gauss's elimination method. */
zVec zLESolveGauss(const zMat a, const zVec b, zVec ans)
{
  zMat acp;
  zVec bcp, s;
  zIndex idx;

  if( !zMatIsSqr( a ) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return NULL;
  }
  if( !zMatColVecSizeEqual( a, b ) || !zVecSizeEqual( ans, b ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  acp = zMatClone( a );
  bcp = zVecClone( b );
  s = zVecAlloc( zVecSizeNC(b) );
  idx = zIndexCreate( zVecSizeNC(b) );
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

/* a solver for Ly=b. */
zVec zLESolveL(const zMat l, const zVec b, zVec ans, const zIndex idx)
{
  int i, j, p;
  double x;

  for( i=0; i<zIndexSizeNC(idx); i++ ){
    x = zVecElemNC( b, (p=zIndexElemNC(idx,i)) );
    for( j=0; j<i; j++ )
      x -= zMatElemNC(l,p,j)*zVecElemNC(ans,j);
    zVecSetElemNC( ans, i, x/zMatElemNC(l,p,i) );
  }
  return ans;
}

/* a solver for Ux=y. */
zVec zLESolveU(const zMat u, const zVec b, zVec ans)
{
  int i, j;
  double x;

  for( i=zVecSizeNC(b)-1; i>=0; i-- ){
    x = zVecElemNC( b, i );
    for( j=zVecSizeNC(b)-1; j>i; j-- )
      x -= zMatElemNC(u,i,j)*zVecElemNC(ans,j);
    zVecSetElemNC( ans, i, x );
  }
  return ans;
}

/* a solver for LUx=b. */
zVec zLESolveLU(const zMat l, const zMat u, const zVec b, zVec ans, const zIndex idx)
{
  zVec c;

  if( !zMatIsSqr(l) || !zMatIsSqr(u) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return NULL;
  }
  if( zMatRowSize(l) != zIndexSizeNC(idx) ||
      zMatRowSize(u) != zIndexSizeNC(idx) ||
      zVecSize(b) != zIndexSizeNC(idx) ||
      zVecSize(ans) != zIndexSizeNC(idx) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  if( !( c = zVecAlloc( zIndexSizeNC(idx) ) ) ) return NULL;
  zLESolveL( l, b, c, idx );
  zLESolveU( u, c, ans );
  zVecFree( c );
  return ans;
}

/* linear equation solver: Residual iteration on LU decomposition. */
zVec zLESolveLURI(const zMat a, const zVec b, zVec ans)
{
  int i;
  zMat l, u;
  zVec res, err;
  double err_norm, err_norm_old = HUGE_VAL;
  zIndex idx;

  l = zMatAllocSqr( zVecSizeNC(b) );
  u = zMatAllocSqr( zVecSizeNC(b) );
  res = zVecAlloc( zVecSizeNC(b) );
  err = zVecAlloc( zVecSizeNC(b) );
  idx = zIndexCreate( zVecSizeNC(b) );
  if( !l || !u || !res || !err || !idx ) goto TERMINATE;

  if( zMatDecompLU( a, l, u, idx ) < zMatRowSizeNC(a) ){
    ZRUNERROR( ZM_ERR_MAT_SINGULAR );
    ans = NULL;
    goto TERMINATE;
  }
  zLESolveLU( l, u, b, ans, idx );
  for( i=0; i<Z_MAX_ITER_NUM; i++ ){
    zLEResidual( a, b, ans, res );
    err_norm = zVecNorm( res );
    if( err_norm >= err_norm_old ) goto TERMINATE;
    err_norm_old = err_norm;
    zLESolveLU( l, u, res, err, idx );
    zVecAddNCDRC( ans, err );
  }
  ZITERWARN( Z_MAX_ITER_NUM );

 TERMINATE:
  zMatFree( l );
  zMatFree( u );
  zVecFree( res );
  zVecFree( err );
  zIndexFree( idx );
  return ans;
}

/* linear equation solver: Gauss-Seidel's method. */
zVec zLESolveGaussSeidel(const zMat a, const zVec b, zVec ans)
{
  int i, j, k, p, count;
  double x;
  zIndex idx;

  if( !zMatIsSqr(a) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return NULL;
  }
  if( !zMatColVecSizeEqual(a,ans) ||
      !zMatRowVecSizeEqual(a,b) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  if( !( idx = zIndexCreate(zVecSizeNC(ans)) ) ) return NULL;
  for( i=0; i<zIndexSizeNC(idx); i++ )
    zMatPivoting( a, idx, i, i );

  for( i=0; ; i++ ){
    for( count=0, j=0; j<zVecSizeNC(b); j++ ){
      p = zIndexElemNC(idx,j);
      x = zVecElemNC(b,p);
      for( k=0; k<zVecSizeNC(ans); k++ )
        if( k != j ) x -= zMatElemNC(a,p,k)*zVecElemNC(ans,k);
      x /= zMatElemNC(a,p,j);
      if( zIsTiny( x - zVecElemNC(ans,j) ) ) count++;
      zVecSetElemNC( ans, j, x );
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
