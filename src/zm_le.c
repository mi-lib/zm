/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le - linear equation.
 */

#include <zm/zm_le.h>

/* directly make a matrix column-balanced, where each column is divided by the power of two not greater
 * than the absolute-maximum component in the column, so that the balancing is exact. */
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
    if( zVecElemNC(s,i) == 0 ){
      zVecSetElemNC( s, i, 1.0 );
      continue;
    }
    /* inverse column-balancing factor, which is also a power of two */
    zVecSetElemNC( s, i, 1.0 / zPow2Floor( zVecElemNC(s,i) ) );
    for( j=0; j<zMatRowSizeNC(m); j++ )
      zMatElemNC(m,j,i) *= zVecElemNC(s,i);
  }
}

/* the power of two not greater than the absolute-maximum component of a row of a matrix. */
#define _zMatRowBalancingFactor(mp,m) zPow2Floor( zDataAbsMax( mp, zMatColSizeNC(m), NULL ) )

/* directly make a pair of matrix and vector row-balanced. */
void zMatVecRowBalancingDST(zMat m, zVec v)
{
  int i;
  double *mp, tmp;

  for( mp=zMatBuf(m), i=0; i<zMatRowSizeNC(m); mp+=zMatColCapacity(m), i++ ){
    if( ( tmp = _zMatRowBalancingFactor( mp, m ) ) == 0 ) continue;
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

  for( mp1=zMatBuf(m1), mp2=zMatBuf(m2), i=0; i<zMatRowSizeNC(m1); mp1+=zMatColCapacity(m1), mp2+=zMatColCapacity(m2), i++ ){
    if( ( tmp = _zMatRowBalancingFactor( mp1, m1 ) ) == 0 ) continue;
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

  for( mp1=zMatBuf(m1), mp2=zMatBuf(m2), i=0; i<zMatRowSizeNC(m1); mp1+=zMatColCapacity(m1), mp2+=zMatColCapacity(m2), i++ ){
    if( ( tmp = _zMatRowBalancingFactor( mp1, m1 ) ) == 0 ) continue;
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
zVec zLESolveGaussDST(zMat a, zVec b, zVec ans, zIndex index, zVec s)
{
  int i, j, k, n, p, q;
  double ahead;
  double x;

  n = zVecSizeNC( b );
  zMatVecBalancingDST( a, b, s );
  /* forward elimination */
  for( i=0; i<n; i++ ){
    p = zMatPivoting( a, index, i, i );
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
      q = zIndexElemNC( index, j );
      ahead = zMatElemNC(a,q,i);
      for( k=i+1; k<n; k++ )
        zMatElemNC(a,q,k) -= zMatElemNC(a,p,k) * ahead;
      zVecElemNC(b,q) -= zVecElemNC(b,p) * ahead;
      zMatSetElemNC( a, q, i, 0 );
    }
  }
  /* backward elimination */
  for( i=n-1; i>=0; i-- ){
    x = zVecElemNC( b, ( p = zIndexElemNC(index,i) ) );
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
  zIndex index;

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
  index = zIndexCreate( zVecSizeNC(b) );
  if( acp && bcp && index && s )
    ans = zLESolveGaussDST( acp, bcp, ans, index, s );
  else
    ZALLOCERROR();
  zMatFree( acp );
  zVecFree( bcp );
  zVecFree( s );
  zIndexFree( index );
  return ans;
}

/* a solver for Ly=b. */
zVec zLESolveL(const zMat l, const zVec b, zVec ans, const zIndex index)
{
  int i, j, p;
  double x;

  for( i=0; i<zIndexSizeNC(index); i++ ){
    x = zVecElemNC( b, (p=zIndexElemNC(index,i)) );
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

/* a solver for Ux=y, where U is an upper triangular matrix with columns permuted by index. */
zVec zLESolveUCP(const zMat u, const zVec b, zVec ans, const zIndex index)
{
  int i, j, c;
  double x;

  for( i=zVecSizeNC(b)-1; i>=0; i-- ){
    x = zVecElemNC( b, i );
    for( j=zVecSizeNC(b)-1; j>i; j-- ){
      c = zIndexElemNC(index,j);
      x -= zMatElemNC(u,i,c)*zVecElemNC(ans,c);
    }
    zVecSetElemNC( ans, zIndexElemNC(index,i), x );
  }
  return ans;
}

/* a solver for L L^T x = b with the Cholesky factor L, where the rows index[0], index[1], ... of L are lower
   triangular in this order. the intermediate and final answers are put at the permuted positions of ans,
   so that b and ans can be the same vector. */
zVec zLESolveCholesky(const zMat l, const zVec b, zVec ans, const zIndex index)
{
  int i, j, n;
  double x;

  if( !zMatIsSqr(l) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return NULL;
  }
  if( zMatRowSize(l) != zIndexSizeNC(index) ||
      zVecSize(b) != zIndexSizeNC(index) ||
      zVecSize(ans) != zIndexSizeNC(index) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  n = zIndexSizeNC(index);
  for( i=0; i<n; i++ ){ /* forward substitution */
    x = zVecElemNC( b, zIndexElemNC(index,i) );
    for( j=0; j<i; j++ )
      x -= zMatElemNC(l,zIndexElemNC(index,i),j) * zVecElemNC(ans,zIndexElemNC(index,j));
    zVecSetElemNC( ans, zIndexElemNC(index,i), x / zMatElemNC(l,zIndexElemNC(index,i),i) );
  }
  for( i=n-1; i>=0; i-- ){ /* backward substitution */
    x = zVecElemNC( ans, zIndexElemNC(index,i) );
    for( j=i+1; j<n; j++ )
      x -= zMatElemNC(l,zIndexElemNC(index,j),i) * zVecElemNC(ans,zIndexElemNC(index,j));
    zVecSetElemNC( ans, zIndexElemNC(index,i), x / zMatElemNC(l,zIndexElemNC(index,i),i) );
  }
  return ans;
}

/* a solver for LUx=b. */
zVec zLESolveLU(const zMat l, const zMat u, const zVec b, zVec ans, const zIndex index)
{
  zVec c;

  if( !zMatIsSqr(l) || !zMatIsSqr(u) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return NULL;
  }
  if( zMatRowSize(l) != zIndexSizeNC(index) ||
      zMatRowSize(u) != zIndexSizeNC(index) ||
      zVecSize(b) != zIndexSizeNC(index) ||
      zVecSize(ans) != zIndexSizeNC(index) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  if( !( c = zVecAlloc( zIndexSizeNC(index) ) ) ) return NULL;
  zLESolveL( l, b, c, index );
  zLESolveU( u, c, ans );
  zVecFree( c );
  return ans;
}

/* residual b - a x computed as accurately as in twice the working precision, based on Dot2 algorithm
 * by T. Ogita, S. M. Rump and S. Oishi:
 *  T. Ogita, S. M. Rump and S. Oishi, Accurate Sum and Dot Product, SIAM J. Sci. Comput.,
 *  26(6):1955-1988, 2005.
 * The components of a and x are assumed to be less than about 1e+300 in the absolute values, so that
 * zMul() does not overflow. */
static zVec _zLEResidualDot2(const zMat a, const zVec b, const zVec x, zVec res)
{
  int i, j;
  double p, q, h, r, c;

  for( i=0; i<zMatRowSizeNC(a); i++ ){
    p = zVecElemNC(b,i);
    for( c=0, j=0; j<zMatColSizeNC(a); j++ ){
      h = zMul( -zMatElemNC(a,i,j), zVecElemNC(x,j), &r );
      p = zAdd( p, h, &q );
      c += q + r;
    }
    zVecSetElemNC( res, i, p + c );
  }
  return res;
}

/* linear equation solver: Residual iteration on LU decomposition. */
zVec zLESolveLURI(const zMat a, const zVec b, zVec ans)
{
  int i;
  zMat ab, l, u;
  zVec bb, s, res, err;
  double err_norm, err_norm_old = HUGE_VAL;
  zIndex index;
  zVec ret = NULL;

  ab = zMatClone( a );
  bb = zVecClone( b );
  s = zVecAlloc( zVecSizeNC(b) );
  l = zMatAllocSqr( zVecSizeNC(b) );
  u = zMatAllocSqr( zVecSizeNC(b) );
  res = zVecAlloc( zVecSizeNC(b) );
  err = zVecAlloc( zVecSizeNC(b) );
  index = zIndexCreate( zVecSizeNC(b) );
  if( !ab || !bb || !s || !l || !u || !res || !err || !index ) goto TERMINATE;
  /* the balancing by powers of two is exact, so that the residual of the balanced equation is that of
   * the original equation scaled exactly. */
  zMatVecBalancingDST( ab, bb, s );
  if( zMatDecompLU( ab, l, u, index ) < zMatRowSizeNC(ab) ){
    ZRUNERROR( ZM_ERR_MAT_SINGULAR );
    goto TERMINATE;
  }
  zLESolveLU( l, u, bb, ans, index );
  /* the iteration terminates when the correction becomes negligible or stops decreasing. */
  for( i=0; i<Z_MAX_ITER_NUM; i++ ){
    _zLEResidualDot2( ab, bb, ans, res );
    zLESolveLU( l, u, res, err, index );
    if( ( err_norm = zVecNorm( err ) ) >= err_norm_old ) break;
    zVecAddNCDRC( ans, err );
    if( err_norm <= DBL_EPSILON * zVecNorm( ans ) ) break;
    err_norm_old = err_norm;
  }
  if( i == Z_MAX_ITER_NUM ) ZITERWARN( Z_MAX_ITER_NUM );
  ret = zVecAmpNC( ans, s, ans );

 TERMINATE:
  zMatFreeAtOnce( 3, ab, l, u );
  zVecFreeAtOnce( 4, bb, s, res, err );
  zIndexFree( index );
  return ret;
}

/* linear equation solver: Gauss-Seidel's method. */
zVec zLESolveGaussSeidel(const zMat a, const zVec b, zVec ans)
{
  int i, j, k, p;
  double x, dx_max, x_max;
  zIndex index;

  if( !zMatIsSqr(a) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return NULL;
  }
  if( !zMatColVecSizeEqual(a,ans) ||
      !zMatRowVecSizeEqual(a,b) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  if( !( index = zIndexCreate(zVecSizeNC(ans)) ) ) return NULL;
  for( i=0; i<zIndexSizeNC(index); i++ )
    zMatPivoting( a, index, i, i );

  for( i=0; ; i++ ){
    for( dx_max=x_max=0, j=0; j<zVecSizeNC(b); j++ ){
      p = zIndexElemNC(index,j);
      x = zVecElemNC(b,p);
      for( k=0; k<zVecSizeNC(ans); k++ )
        if( k != j ) x -= zMatElemNC(a,p,k)*zVecElemNC(ans,k);
      x /= zMatElemNC(a,p,j);
      dx_max = zMax( dx_max, fabs( x - zVecElemNC(ans,j) ) );
      x_max = zMax( x_max, fabs( x ) );
      zVecSetElemNC( ans, j, x );
    }
    /* converged if the change is negligible relative to the solution */
    if( dx_max <= zTOL * x_max ) break;
    if( i == Z_MAX_ITER_NUM ){
      ZITERWARN( Z_MAX_ITER_NUM );
      break;
    }
  }
  zIndexFree( index );
  return ans;
}
