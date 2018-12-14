/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_lu - linear equation: LU decomposition.
 */

#include <zm/zm_le.h>

/* zLUDecompDST
 * - LU decomposition of matrix (destructive).
 */
int zLUDecompDST(zMat m, zMat l, zMat u, zIndex idx)
{
  register int r, c, i, j;
  int p, q;
  double ahead;

  zMatClear( l );
  zMatClear( u );
  for( r=c=0; r<zMatRowSizeNC(m) && c<zMatColSizeNC(m); c++ ){
    p = zPivoting( m, idx, r, c );
    if( zIsTiny( ( ahead = zMatElem( m, p, c ) ) ) ) continue;

    zMatSetElem( l, p, r, ahead );
    zMatSetElem( u, r, c, 1 );
    for( i=r+1; i<zMatRowSizeNC(l); i++ ){ /* L column */
      q = zIndexElem( idx, i );
      zMatSetElem( l, q, c, zMatElem( m, q, c ) );
    }
    for( i=c+1; i<zMatColSizeNC(u); i++ ) /* U row */
      zMatSetElem( u, r, i, zMatElem( m, p, i ) / ahead );
    for( j=c+1; j<zMatColSizeNC(m); j++ ) /* other mod. */
      for( i=r+1; i<zMatRowSizeNC(m); i++ ){
        q = zIndexElem( idx, i );
        zMatElem(m,q,j) -= zMatElem(l,q,c) * zMatElem(u,r,j);
      }
    r++;
  }
  return r; /* rank */
}

/* zLUDecomp
 * - LU decomposition of matrix.
 */
int zLUDecomp(zMat m, zMat l, zMat u, zIndex idx)
{
  int n, rank;
  zMat mc;

  if( !zMatIsSqr(l) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return -1;
  }
  if( !zMatRowSizeIsEqual(m,l) ||
      !zMatSizeIsEqual(m,u) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return -1;
  }
  if( zMatRowSize(m) != ( n = zArrayNum(idx) ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return -1;
  }
  zIndexOrder( idx, 0 );
  if( !( mc = zMatClone(m) ) ) return -1;

  rank = zLUDecompDST( mc, l, u, idx );
  zMatFree( mc );
  return rank;
}

/* zLUDecompReg
 * - LU decomposition and regression.
 */
int zLUDecompReg(zMat m, zMat l, zMat u, zIndex idx)
{
  int rank;

  if( ( rank = zLUDecomp( m, l, u, idx ) ) < zMatRowSizeNC(m) ){
    zMatColReg( l, rank ); /* lower triangle regeression */
    zMatRowReg( u, rank ); /* upper triangle regeression */
  }
  return rank;
}

/* zLUDecompAlloc
 * - LU decomposition with an automatic matrix allocation and resize.
 */
int zLUDecompAlloc(zMat m, zMat *l, zMat *u, zIndex *idx)
{
  *l = zMatAllocSqr( zMatRowSizeNC(m) );
  *u = zMatAlloc( zMatRowSizeNC(m), zMatColSizeNC(m) );
  *idx = zIndexCreate( zMatRowSizeNC(m) );
  if( !*l || !*u || !*idx ){
    zMatFree( *l );
    zMatFree( *u );
    zIndexFree( *idx );
    return 0;
  }
  return zLUDecompReg( m, *l, *u, *idx );
}

/* zCholeskyDecompDST
 * - Cholesky decomposition of matrix (destructive).
 */
int zCholeskyDecompDST(zMat m, zMat l, zIndex idx)
{
  register int i, j, k;
  int n, p, rank;
  double a;

  n = zArrayNum( idx );
  zMatClear( l );
  for( rank=0, i=0; i<n; i++ ){
    p = zPivotingDiag( m, idx, i );
    if( zIsTiny( ( a = zMatElem(m,p,p) ) ) ){
      ZRUNWARN( ZM_ERR_LE_SINGULAR );
      a = 0;
    } else
    if( a < 0 ){
      ZRUNERROR( ZM_ERR_LE_CHOLESKY );
      return -1;
    } else{
      zMatSetElem( m, p, p, ( a = sqrt( a ) ) );
      a = 1.0 / a;
      rank++;
    }

    for( j=i+1; j<n; j++ ){
      zMatElem(m,p,zIndexElem(idx,j)) *= a;
      zMatElem(m,zIndexElem(idx,j),p) *= a;
    }
    for( j=i+1; j<n; j++ )
      for( k=i+1; k<n; k++ )
        zMatElem(m,zIndexElem(idx,j),zIndexElem(idx,k))
          -= zMatElem(m,zIndexElem(idx,j),p)
            * zMatElem(m,p,zIndexElem(idx,k));
  }
  for( i=0; i<n; i++ )
    for( j=0; j<=i; j++ )
      zMatSetElem( l, zIndexElem(idx,i), j,
        zMatElem( m, zIndexElem(idx,i), zIndexElem(idx,j) ) );
  return rank;
}

/* zCholeskyDecomp
 * - Cholesky decomposition of matrix
 */
int zCholeskyDecomp(zMat m, zMat l, zIndex idx)
{
  int rank;
  zMat mc;

  if( !zMatIsSqr(m) || !zMatIsSqr(l) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return false;
  }
  if( zMatRowSize(m) != zArrayNum(idx) ||
      zMatRowSize(l) != zArrayNum(idx) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return false;
  }
  zIndexOrder( idx, 0 );
  if( !( mc = zMatClone( m ) ) ) return -1;
  rank = zCholeskyDecompDST( mc, l, idx );
  zMatFree( mc );
  return rank;
}

/* zCholeskyDecompReg
 * - Cholesky decomposition and regression.
 */
int zCholeskyDecompReg(zMat m, zMat l, zIndex idx)
{
  int rank;

  rank = zCholeskyDecomp( m, l, idx );
  if( rank > 0 && rank < zMatRowSizeNC(m) )
    zMatColReg( l, rank ); /* lower triangle regeression */
  return rank;
}

/* zCholeskyDecompAlloc
 * - Cholesky decomposition with an automatic matrix allocation and resize.
 */
int zCholeskyDecompAlloc(zMat m, zMat *l, zIndex *idx)
{
  *l = zMatAllocSqr( zMatRowSizeNC(m) );
  *idx = zIndexCreate( zMatRowSizeNC(m) );
  if( !*l || !*idx ){
    zMatFree( *l );
    zIndexFree( *idx );
    return 0;
  }
  return zCholeskyDecompReg( m, *l, *idx );
}
