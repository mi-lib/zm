/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_lu - linear equation: LU decomposition.
 */

#include <zm/zm_le.h>

/* LU decomposition of a matrix (destructive). */
int zMatDecompLUDST(zMat m, zMat l, zMat u, zIndex idx)
{
  int r, c, i, j, p, q;
  double ahead;

  zMatZero( l );
  zMatZero( u );
  for( r=c=0; r<zMatRowSizeNC(m) && c<zMatColSizeNC(m); c++ ){
    p = zPivoting( m, idx, r, c );
    if( zIsTiny( ( ahead = zMatElemNC( m, p, c ) ) ) ) continue;

    zMatSetElemNC( l, p, r, ahead );
    zMatSetElemNC( u, r, c, 1 );
    for( i=r+1; i<zMatRowSizeNC(l); i++ ){ /* L column */
      q = zIndexElemNC( idx, i );
      zMatSetElemNC( l, q, c, zMatElemNC( m, q, c ) );
    }
    for( i=c+1; i<zMatColSizeNC(u); i++ ) /* U row */
      zMatSetElemNC( u, r, i, zMatElemNC( m, p, i ) / ahead );
    for( j=c+1; j<zMatColSizeNC(m); j++ ) /* other mod. */
      for( i=r+1; i<zMatRowSizeNC(m); i++ ){
        q = zIndexElemNC( idx, i );
        zMatElemNC(m,q,j) -= zMatElemNC(l,q,c) * zMatElemNC(u,r,j);
      }
    r++;
  }
  return r; /* rank */
}

/* LU decomposition of a matrix. */
int zMatDecompLU(zMat m, zMat l, zMat u, zIndex idx)
{
  int rank;
  zMat mc;

  if( !zMatIsSqr(l) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return -1;
  }
  if( !zMatRowSizeIsEqual(m,l) || !zMatSizeIsEqual(m,u) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return -1;
  }
  if( zMatRowSize(m) != zIndexSizeNC(idx) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return -1;
  }
  zIndexOrder( idx, 0 );
  if( !( mc = zMatClone(m) ) ) return -1;

  rank = zMatDecompLUDST( mc, l, u, idx );
  zMatFree( mc );
  return rank;
}

/* LU decomposition and regression of a matrix. */
int zMatDecompLUReg(zMat m, zMat l, zMat u, zIndex idx)
{
  int rank;

  if( ( rank = zMatDecompLU( m, l, u, idx ) ) < 0 ) return -1;
  if( rank < (int)zMatRowSizeNC(m) ){
    zMatColReg( l, rank ); /* lower triangle regression */
    zMatRowReg( u, rank ); /* upper triangle regression */
  }
  return rank;
}

/* LU decomposition with an automatic matrix allocation and resize. */
int zMatDecompLUAlloc(zMat m, zMat *l, zMat *u, zIndex *idx)
{
  *l = zMatAllocSqr( zMatRowSizeNC(m) );
  *u = zMatAlloc( zMatRowSizeNC(m), zMatColSizeNC(m) );
  *idx = zIndexCreate( zMatRowSizeNC(m) );
  if( !*l || !*u || !*idx ){
    zMatFree( *l );
    zMatFree( *u );
    zIndexFree( *idx );
    return -1;
  }
  return zMatDecompLUReg( m, *l, *u, *idx );
}

/* Cholesky decomposition of a matrix (destructive). */
int zMatDecompCholeskyDST(zMat m, zMat l, zIndex idx)
{
  int i, j, k, n, p, rank;
  double a;

  n = zIndexSizeNC(idx);
  zMatZero( l );
  for( rank=0, i=0; i<n; i++ ){
    p = zPivotingDiag( m, idx, i );
    if( zIsTiny( ( a = zMatElemNC(m,p,p) ) ) ){
      ZRUNWARN( ZM_ERR_LE_SINGULAR );
      a = 0;
    } else
    if( a < 0 ){
      ZRUNERROR( ZM_ERR_LE_CHOLESKY );
      return -1;
    } else{
      zMatSetElemNC( m, p, p, ( a = sqrt( a ) ) );
      a = 1.0 / a;
      rank++;
    }

    for( j=i+1; j<n; j++ ){
      zMatElemNC(m,p,zIndexElemNC(idx,j)) *= a;
      zMatElemNC(m,zIndexElemNC(idx,j),p) *= a;
    }
    for( j=i+1; j<n; j++ )
      for( k=i+1; k<n; k++ )
        zMatElemNC(m,zIndexElemNC(idx,j),zIndexElemNC(idx,k))
          -= zMatElemNC(m,zIndexElemNC(idx,j),p)
            * zMatElemNC(m,p,zIndexElemNC(idx,k));
  }
  for( i=0; i<n; i++ )
    for( j=0; j<=i; j++ )
      zMatSetElemNC( l, zIndexElemNC(idx,i), j,
        zMatElemNC( m, zIndexElemNC(idx,i), zIndexElemNC(idx,j) ) );
  return rank;
}

/* Cholesky decomposition of a matrix */
int zMatDecompCholesky(zMat m, zMat l, zIndex idx)
{
  int rank;
  zMat mc;

  if( !zMatIsSqr(m) || !zMatIsSqr(l) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return -1;
  }
  if( zMatRowSize(m) != zIndexSizeNC(idx) ||
      zMatRowSize(l) != zIndexSizeNC(idx) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return -1;
  }
  zIndexOrder( idx, 0 );
  if( !( mc = zMatClone( m ) ) ) return -1;
  rank = zMatDecompCholeskyDST( mc, l, idx );
  zMatFree( mc );
  return rank;
}

/* Cholesky decomposition and regression of a matrix. */
int zMatDecompCholeskyReg(zMat m, zMat l, zIndex idx)
{
  int rank;

  if( ( rank = zMatDecompCholesky( m, l, idx ) ) < 0 ) return -1;
  if( rank < (int)zMatRowSizeNC(m) )
    zMatColReg( l, rank ); /* lower triangle regression */
  return rank;
}

/* Cholesky decomposition with an automatic matrix allocation and resize. */
int zMatDecompCholeskyAlloc(zMat m, zMat *l, zIndex *idx)
{
  *l = zMatAllocSqr( zMatRowSizeNC(m) );
  *idx = zIndexCreate( zMatRowSizeNC(m) );
  if( !*l || !*idx ){
    zMatFree( *l );
    zIndexFree( *idx );
    return -1;
  }
  return zMatDecompCholeskyReg( m, *l, *idx );
}
