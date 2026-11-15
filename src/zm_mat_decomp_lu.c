/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_mat_decomp_lu - matrix decomposition: LU decomposition.
 */

#include <zm/zm_mat_decomp.h>

/* find the pivot on the c-th column among the r-th and the following rows of m, where the i-th row
 * means the index[i]-th row of m if index is given, or the i-th row of m otherwise.
 * If the maximum absolute value of the candidates is not greater than tol, the column is regarded as
 * dependent on the preceding ones, the candidates are zeroed, and -1 is returned. */
static int _zMatDecompLUPivot(zMat m, zIndex index, int r, int c, double tol)
{
  int i, p;
  double max, tmp;

#define _zmat_decomp_row_pivot(i) ( index ? zIndexElemNC(index,i) : (i) )
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
static int _zMatDecompLUDST_U(zMat l, zMat u, zIndex index, double tol)
{
  int r, c, i, p, size;
  double ahead, *ur, *ui;

  zMatZero( l );
  zIndexOrder( index, 0 );
  for( r=c=0; r<zMatRowSizeNC(u) && c<zMatColSizeNC(u); c++ ){
    /* partial pivoting with a physical row swap */
    if( ( p = _zMatDecompLUPivot( u, NULL, r, c, tol ) ) < 0 ) continue;
    if( p != r ){
      zMatSwapRowNC( u, r, p );
      zIndexSwapNC( index, r, p );
    }
    ur = zMatRowBufNC(u,r);
    ahead = ur[c];
    zMatSetElemNC( l, zIndexElemNC(index,r), r, ahead );
    ur[c] = 1;
    size = zMatColSizeNC(u) - c - 1;
    for( i=c+1; i<zMatColSizeNC(u); i++ ) ur[i] /= ahead; /* U row */
    for( i=r+1; i<zMatRowSizeNC(u); i++ ){ /* L column and elimination of the other rows */
      ui = zMatRowBufNC(u,i);
      zMatSetElemNC( l, zIndexElemNC(index,i), r, ui[c] );
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
static int _zMatDecompLUDST_L(zMat l, zMat u, zIndex index, double tol)
{
  int r, c, i, j, p, size;
  double ahead, *lp, *ur;

  zMatZero( u );
  zIndexOrder( index, 0 );
  for( r=c=0; r<zMatRowSizeNC(l) && c<zMatColSizeNC(l); c++ ){
    /* partial pivoting (rows are not physically swapped to keep the correspondence of rows of l to m) */
    if( ( p = _zMatDecompLUPivot( l, index, r, c, tol ) ) < 0 ) continue;
    zIndexSwapNC( index, r, p );
    lp = zMatRowBufNC(l,zIndexElemNC(index,r));
    ahead = lp[c];
    ur = zMatRowBufNC(u,r);
    ur[c] = 1;
    for( j=c+1; j<zMatColSizeNC(l); j++ ) ur[j] = lp[j] / ahead; /* U row */
    lp[r] = ahead; /* r <= c, where the column r has already been processed. */
    zRawVecZero( lp + r + 1, zMatColSizeNC(l) - r - 1 );
    size = zMatColSizeNC(l) - c - 1;
    for( i=r+1; i<zMatRowSizeNC(l); i++ ){ /* L column and elimination of the other rows */
      lp = zMatRowBufNC(l,zIndexElemNC(index,i));
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
int zMatDecompLU(const zMat m, zMat l, zMat u, zIndex index)
{
  double tol;

  if( zMatRowSize(m) != zIndexSizeNC(index) ){
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
    return _zMatDecompLUDST_L( l, u, index, tol );
  }
  /* column-long rectangular matrix */
  if( !zMatSizeEqual( m, u ) || !zMatRowSizeEqual( m, l ) || !zMatIsSqr( l ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return -1;
  }
  zMatCopyNC( m, u );
  return _zMatDecompLUDST_U( l, u, index, tol );
}

/* LU decomposition with an automatic matrix allocation. */
int zMatDecompLUAlloc(const zMat m, zMat *l, zMat *u, zIndex *index)
{
  int minsize;

  minsize = zMin( zMatRowSizeNC(m), zMatColSizeNC(m) );
  *l = zMatAlloc( zMatRowSizeNC(m), minsize );
  *u = zMatAlloc( minsize, zMatColSizeNC(m) );
  *index = zIndexCreate( zMatRowSizeNC(m) );
  if( !*l || !*u || !*index ){
    zMatFree( *l );
    zMatFree( *u );
    zIndexFree( *index );
    *l = *u = NULL;
    *index = NULL;
    return -1;
  }
  return zMatDecompLU( m, *l, *u, *index );
}

/* rearrange columns of a matrix as l'[i][j] = l[i][index[j]], and zero the upper part of the pivoted rows. */
static void _zMatDecompCholeskyArrange(zMat l, zIndex index)
{
  int i, j, j0, k, n;
  double tmp, *lp;

  n = zIndexSizeNC(index);
  for( j0=0; j0<n; j0++ ){
    if( zIndexElemNC(index,j0) < 0 ) continue; /* already rearranged */
    for( i=0; i<n; i++ ){ /* rotate a cycle of the permutation in each row */
      lp = zMatRowBufNC(l,i);
      tmp = lp[j0];
      for( j=j0; ( k = zIndexElemNC(index,j) ) != j0; j=k ) lp[j] = lp[k];
      lp[j] = tmp;
    }
    for( j=j0; ( k = zIndexElemNC(index,j) ) >= 0; j=k ) /* mark the cycle */
      zIndexSetElemNC( index, j, ~k );
  }
  for( j=0; j<n; j++ ) /* unmark */
    zIndexSetElemNC( index, j, ~zIndexElemNC(index,j) );
  for( i=0; i<n; i++ )
    zRawVecZero( zMatRowBufNC(l,zIndexElemNC(index,i)) + i + 1, n - i - 1 );
}

/* Cholesky decomposition of a matrix (destructive). l has to be initially set for the original matrix. */
static int _zMatDecompCholeskyDST(zMat l, zIndex index, double tol)
{
  int i, j, k, n, p, rank;
  double a;

  n = zIndexSizeNC(index);
  zIndexOrder( index, 0 );
  for( rank=0, i=0; i<n; i++ ){
    p = zMatPivotingDiag( l, index, i );
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
      zMatElemNC(l,p,zIndexElemNC(index,j)) *= a;
      zMatElemNC(l,zIndexElemNC(index,j),p) *= a;
    }
    for( j=i+1; j<n; j++ )
      for( k=i+1; k<n; k++ )
        zMatElemNC(l,zIndexElemNC(index,j),zIndexElemNC(index,k))
          -= zMatElemNC(l,zIndexElemNC(index,j),p) * zMatElemNC(l,p,zIndexElemNC(index,k));
  }
  _zMatDecompCholeskyArrange( l, index );
  if( rank < n ) /* automatically adjust size of the factorized matrix */
    zMatColResize( l, rank );
  return rank;
}

/* Cholesky decomposition of a matrix */
int zMatDecompCholesky(const zMat m, zMat l, zIndex index)
{
  int i;
  double d, dmax = 0;

  if( !zMatIsSqr( m ) || !zMatIsSqr( l ) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return -1;
  }
  if( zMatRowSize(m) != zIndexSizeNC(index) || zMatRowSize(l) != zIndexSizeNC(index) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return -1;
  }
  for( i=0; i<zMatRowSizeNC(m); i++ ) /* maximum absolute diagonal component */
    if( ( d = fabs( zMatElemNC(m,i,i) ) ) > dmax ) dmax = d;
  zMatCopyNC( m, l );
  return _zMatDecompCholeskyDST( l, index, ZM_CHOLESKY_DECOMP_TOL * dmax );
}

/* Cholesky decomposition with an automatic matrix allocation. */
int zMatDecompCholeskyAlloc(const zMat m, zMat *l, zIndex *index)
{
  *l = zMatAllocSqr( zMatRowSizeNC(m) );
  *index = zIndexCreate( zMatRowSizeNC(m) );
  if( !*l || !*index ){
    ZALLOCERROR();
    zMatFree( *l );
    zIndexFree( *index );
    *l = NULL;
    *index = NULL;
    return -1;
  }
  return zMatDecompCholesky( m, *l, *index );
}
