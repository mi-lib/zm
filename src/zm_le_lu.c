/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_lu - linear equation: LU decomposition.
 */

#include <zm/zm_le.h>

/* LU decomposition based on Crout method (destructive). u has to be initially set for the original matrix. */
static int _zMatDecompLUDST(zMat l, zMat u, zIndex idx, double tol)
{
  int r, c, i, p, size;
  double ahead, max, tmp, *ur, *ui;

  zMatZero( l );
  zIndexOrder( idx, 0 );
  for( r=c=0; r<zMatRowSizeNC(u) && c<zMatColSizeNC(u); c++ ){
    /* partial pivoting with a physical row swap */
    for( p=r, max=fabs( zMatElemNC(u,r,c) ), i=r+1; i<zMatRowSizeNC(u); i++ )
      if( ( tmp = fabs( zMatElemNC(u,i,c) ) ) > max ){
        max = tmp;
        p = i;
      }
    if( max <= tol ){ /* the column is dependent on the preceding ones */
      for( i=r; i<zMatRowSizeNC(u); i++ ) zMatSetElemNC( u, i, c, 0 );
      continue;
    }
    if( p != r ){
      zMatSwapRowNC( u, r, p );
      zIndexSwapNC( idx, r, p );
    }
    ur = zMatRowBufNC(u,r);
    ahead = ur[c];
    zMatSetElemNC( l, zIndexElemNC(idx,r), r, ahead );
    ur[c] = 1;
    size = zMatColSizeNC(u) - c - 1;
    for( i=c+1; i<zMatColSizeNC(u); i++ ) ur[i] /= ahead; /* U row */
    for( i=r+1; i<zMatRowSizeNC(u); i++ ){ /* L column and elimination of the other rows */
      ui = zMatRowBufNC(u,i);
      zMatSetElemNC( l, zIndexElemNC(idx,i), r, ui[c] );
      zRawVecCatDRC( ui+c+1, -ui[c], ur+c+1, size );
      ui[c] = 0;
    }
    r++;
  }
  if( r < zMatRowSizeNC(u) ){ /* automatically adjust sizes of factorized matrices */
    zMatColResize( l, r ); /* to be column-full-rank */
    zMatRowResize( u, r ); /* to be row-full-rank */
  }
  return r; /* rank */
}

/* LU decomposition of a matrix. */
int zMatDecompLU(const zMat m, zMat l, zMat u, zIndex idx)
{
  if( !zMatIsSqr( l ) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return -1;
  }
  if( !zMatRowSizeEqual(m,l) || !zMatSizeEqual(m,u) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return -1;
  }
  if( zMatRowSize(m) != zIndexSizeNC(idx) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return -1;
  }
  zMatCopyNC( m, u );
  return _zMatDecompLUDST( l, u, idx, ZM_LU_DECOMP_TOL * zMatElemAbsMax( m, NULL, NULL ) );
}

/* LU decomposition with an automatic matrix allocation. */
int zMatDecompLUAlloc(const zMat m, zMat *l, zMat *u, zIndex *idx)
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
  return zMatDecompLU( m, *l, *u, *idx );
}

/* Cholesky decomposition of a matrix (destructive). */
int zMatDecompCholeskyDST(zMat m, zMat l, zIndex idx)
{
  int i, j, k, n, p, rank;
  double a;

  n = zIndexSizeNC(idx);
  zMatZero( l );
  zIndexOrder( idx, 0 );
  for( rank=0, i=0; i<n; i++ ){
    p = zMatPivotingDiag( m, idx, i );
    if( zIsTiny( ( a = zMatElemNC(m,p,p) ) ) ){
      a = 0;
    } else
    if( a < 0 ){
      ZRUNERROR( ZM_ERR_MAT_NOTPOSITIVESEMIDEFINITE );
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
          -= zMatElemNC(m,zIndexElemNC(idx,j),p) * zMatElemNC(m,p,zIndexElemNC(idx,k));
  }
  for( i=0; i<n; i++ )
    for( j=0; j<=i; j++ )
      zMatSetElemNC( l, zIndexElemNC(idx,i), j,
        zMatElemNC( m, zIndexElemNC(idx,i), zIndexElemNC(idx,j) ) );
  if( rank < zMatRowSizeNC(m) ) /* automatically adjust size of the factorized matrix */
    zMatColResize( l, rank );
  return rank;
}

/* Cholesky decomposition of a matrix */
int zMatDecompCholesky(const zMat m, zMat l, zIndex idx)
{
  int rank;
  zMat mc;

  if( !zMatIsSqr( m ) || !zMatIsSqr( l ) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return -1;
  }
  if( zMatRowSize(m) != zIndexSizeNC(idx) || zMatRowSize(l) != zIndexSizeNC(idx) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return -1;
  }
  if( !( mc = zMatClone( m ) ) ) return -1;
  rank = zMatDecompCholeskyDST( mc, l, idx );
  zMatFree( mc );
  return rank;
}

/* Cholesky decomposition with an automatic matrix allocation. */
int zMatDecompCholeskyAlloc(const zMat m, zMat *l, zIndex *idx)
{
  *l = zMatAllocSqr( zMatRowSizeNC(m) );
  *idx = zIndexCreate( zMatRowSizeNC(m) );
  if( !*l || !*idx ){
    zMatFree( *l );
    zIndexFree( *idx );
    return -1;
  }
  return zMatDecompCholesky( m, *l, *idx );
}
