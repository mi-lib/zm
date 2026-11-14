/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_mat_decomp_lu - matrix decomposition: LU decomposition.
 */

#include <zm/zm_mat_decomp.h>

/* find the pivot on the c-th column among the r-th and the following rows of m, where the i-th row
 * means the idx[i]-th row of m if idx is given, or the i-th row of m otherwise.
 * If the maximum absolute value of the candidates is not greater than tol, the column is regarded as
 * dependent on the preceding ones, the candidates are zeroed, and -1 is returned. */
static int _zMatDecompLUPivot(zMat m, zIndex idx, int r, int c, double tol)
{
  int i, p;
  double max, tmp;

#define _zmat_decomp_row_pivot(i) ( idx ? zIndexElemNC(idx,i) : (i) )
  for( p=r, max=fabs( zMatElemNC(m,_zmat_decomp_row_pivot(r),c) ), i=r+1; i<zMatRowSizeNC(m); i++ )
    if( ( tmp = fabs( zMatElemNC(m,_zmat_decomp_row_pivot(i),c) ) ) > max ){
      max = tmp;
      p = i;
    }
  if( max <= tol ){
    for( i=r; i<zMatRowSizeNC(m); i++ ) zMatSetElemNC( m, _zmat_decomp_row_pivot(i), c, 0 );
    p = -1;
  }
#undef _zmat_decomp_row_pivot
  return p;
}

/* LU decomposition of column-long rectangular matrices based on Crout method (destructive).
 * u has to be initially set for the original matrix. */
static int _zMatDecompLUDST_U(zMat l, zMat u, zIndex idx, double tol)
{
  int r, c, i, p, size;
  double ahead, *ur, *ui;

  zMatZero( l );
  zIndexOrder( idx, 0 );
  for( r=c=0; r<zMatRowSizeNC(u) && c<zMatColSizeNC(u); c++ ){
    /* partial pivoting with a physical row swap */
    if( ( p = _zMatDecompLUPivot( u, NULL, r, c, tol ) ) < 0 ) continue;
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
  return r; /* now, r means rank of the original matrix. */
}

/* LU decomposition of row-long rectangular matrices based on Crout method (destructive).
 * l has to be initially set for the original matrix. */
static int _zMatDecompLUDST_L(zMat l, zMat u, zIndex idx, double tol)
{
  int r, c, i, j, p, size;
  double ahead, *lp, *ur;

  zMatZero( u );
  zIndexOrder( idx, 0 );
  for( r=c=0; r<zMatRowSizeNC(l) && c<zMatColSizeNC(l); c++ ){
    /* partial pivoting (rows are not physically swapped to keep the correspondence of rows of l to m) */
    if( ( p = _zMatDecompLUPivot( l, idx, r, c, tol ) ) < 0 ) continue;
    zIndexSwapNC( idx, r, p );
    lp = zMatRowBufNC(l,zIndexElemNC(idx,r));
    ahead = lp[c];
    ur = zMatRowBufNC(u,r);
    ur[c] = 1;
    for( j=c+1; j<zMatColSizeNC(l); j++ ) ur[j] = lp[j] / ahead; /* U row */
    lp[r] = ahead; /* r <= c, where the column r has already been processed. */
    zRawVecZero( lp + r + 1, zMatColSizeNC(l) - r - 1 );
    size = zMatColSizeNC(l) - c - 1;
    for( i=r+1; i<zMatRowSizeNC(l); i++ ){ /* L column and elimination of the other rows */
      lp = zMatRowBufNC(l,zIndexElemNC(idx,i));
      lp[r] = lp[c];
      if( r < c ) lp[c] = 0;
      zRawVecCatDRC( lp+c+1, -lp[r], ur+c+1, size );
    }
    r++;
  }
  if( r < zMatColSizeNC(l) ){ /* automatically adjust sizes of factorized matrices */
    zMatColResize( l, r ); /* to be column-full-rank */
    zMatRowResize( u, r ); /* to be row-full-rank */
  }
  return r; /* now, r means rank of the original matrix. */
}

/* LU decomposition of a matrix. */
int zMatDecompLU(const zMat m, zMat l, zMat u, zIndex idx)
{
  double tol;

  if( zMatRowSize(m) != zIndexSizeNC(idx) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return -1;
  }
  tol = ZM_LU_DECOMP_TOL * zMatElemAbsMax( m, NULL, NULL );
  if( zMatRowSizeNC(m) > zMatColSizeNC(m) ){ /* row-long rectangular matrix */
    if( !zMatSizeEqual( m, l ) || !zMatRowColSizeEqual( u, l ) || !zMatIsSqr( u ) ){
      ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
      return -1;
    }
    zMatCopyNC( m, l );
    return _zMatDecompLUDST_L( l, u, idx, tol );
  }
  /* column-long rectangular matrix */
  if( !zMatSizeEqual( m, u ) || !zMatRowSizeEqual( m, l ) || !zMatIsSqr( l ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return -1;
  }
  zMatCopyNC( m, u );
  return _zMatDecompLUDST_U( l, u, idx, tol );
}

/* LU decomposition with an automatic matrix allocation. */
int zMatDecompLUAlloc(const zMat m, zMat *l, zMat *u, zIndex *idx)
{
  int minsize;

  minsize = zMin( zMatRowSizeNC(m), zMatColSizeNC(m) );
  *l = zMatAlloc( zMatRowSizeNC(m), minsize );
  *u = zMatAlloc( minsize, zMatColSizeNC(m) );
  *idx = zIndexCreate( zMatRowSizeNC(m) );
  if( !*l || !*u || !*idx ){
    zMatFree( *l );
    zMatFree( *u );
    zIndexFree( *idx );
    *l = *u = NULL;
    *idx = NULL;
    return -1;
  }
  return zMatDecompLU( m, *l, *u, *idx );
}

/* rearrange columns of a matrix as l'[i][j] = l[i][idx[j]], and zero the upper part of the pivoted rows. */
static void _zMatDecompCholeskyArrange(zMat l, zIndex idx)
{
  int i, j, j0, k, n;
  double tmp, *lp;

  n = zIndexSizeNC(idx);
  for( j0=0; j0<n; j0++ ){
    if( zIndexElemNC(idx,j0) < 0 ) continue; /* already rearranged */
    for( i=0; i<n; i++ ){ /* rotate a cycle of the permutation in each row */
      lp = zMatRowBufNC(l,i);
      tmp = lp[j0];
      for( j=j0; ( k = zIndexElemNC(idx,j) ) != j0; j=k ) lp[j] = lp[k];
      lp[j] = tmp;
    }
    for( j=j0; ( k = zIndexElemNC(idx,j) ) >= 0; j=k ) /* mark the cycle */
      zIndexSetElemNC( idx, j, ~k );
  }
  for( j=0; j<n; j++ ) /* unmark */
    zIndexSetElemNC( idx, j, ~zIndexElemNC(idx,j) );
  for( i=0; i<n; i++ )
    zRawVecZero( zMatRowBufNC(l,zIndexElemNC(idx,i)) + i + 1, n - i - 1 );
}

/* Cholesky decomposition of a matrix (destructive). l has to be initially set for the original matrix. */
static int _zMatDecompCholeskyDST(zMat l, zIndex idx, double tol)
{
  int i, j, k, n, p, rank;
  double a;

  n = zIndexSizeNC(idx);
  zIndexOrder( idx, 0 );
  for( rank=0, i=0; i<n; i++ ){
    p = zMatPivotingDiag( l, idx, i );
    if( fabs( ( a = zMatElemNC(l,p,p) ) ) <= tol ){
      a = 0;
    } else
    if( a < 0 ){
      ZRUNERROR( ZM_ERR_MAT_NOTPOSITIVESEMIDEFINITE );
      return -1;
    } else{
      zMatSetElemNC( l, p, p, ( a = sqrt( a ) ) );
      a = 1.0 / a;
      rank++;
    }
    for( j=i+1; j<n; j++ ){
      zMatElemNC(l,p,zIndexElemNC(idx,j)) *= a;
      zMatElemNC(l,zIndexElemNC(idx,j),p) *= a;
    }
    for( j=i+1; j<n; j++ )
      for( k=i+1; k<n; k++ )
        zMatElemNC(l,zIndexElemNC(idx,j),zIndexElemNC(idx,k))
          -= zMatElemNC(l,zIndexElemNC(idx,j),p) * zMatElemNC(l,p,zIndexElemNC(idx,k));
  }
  _zMatDecompCholeskyArrange( l, idx );
  if( rank < n ) /* automatically adjust size of the factorized matrix */
    zMatColResize( l, rank );
  return rank;
}

/* Cholesky decomposition of a matrix */
int zMatDecompCholesky(const zMat m, zMat l, zIndex idx)
{
  int i;
  double d, dmax = 0;

  if( !zMatIsSqr( m ) || !zMatIsSqr( l ) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return -1;
  }
  if( zMatRowSize(m) != zIndexSizeNC(idx) || zMatRowSize(l) != zIndexSizeNC(idx) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return -1;
  }
  for( i=0; i<zMatRowSizeNC(m); i++ ) /* maximum absolute diagonal component */
    if( ( d = fabs( zMatElemNC(m,i,i) ) ) > dmax ) dmax = d;
  zMatCopyNC( m, l );
  return _zMatDecompCholeskyDST( l, idx, ZM_CHOLESKY_DECOMP_TOL * dmax );
}

/* Cholesky decomposition with an automatic matrix allocation. */
int zMatDecompCholeskyAlloc(const zMat m, zMat *l, zIndex *idx)
{
  *l = zMatAllocSqr( zMatRowSizeNC(m) );
  *idx = zIndexCreate( zMatRowSizeNC(m) );
  if( !*l || !*idx ){
    ZALLOCERROR();
    zMatFree( *l );
    zIndexFree( *idx );
    *l = NULL;
    *idx = NULL;
    return -1;
  }
  return zMatDecompCholesky( m, *l, *idx );
}
