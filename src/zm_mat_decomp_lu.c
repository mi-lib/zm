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

/* find the pivot of the complete pivoting among the k-th and the following rows and columns of m,
 * where the i-th row and the j-th column mean the index_row[i]-th row and the index_col[j]-th column
 * of m if the index vectors are given, or the i-th row and the j-th column of m otherwise. */
static double _zMatDecompLUCPPivot(const zMat m, const zIndex index_row, const zIndex index_col, int k, int *p, int *q)
{
  int i, j;
  double max = 0, tmp;

#define _zmat_decomp_row(i) ( index_row ? zIndexElemNC(index_row,i) : (i) )
#define _zmat_decomp_col(j) ( index_col ? zIndexElemNC(index_col,j) : (j) )
  *p = *q = k;
  for( i=k; i<zMatRowSizeNC(m); i++ )
    for( j=k; j<zMatColSizeNC(m); j++ )
      if( ( tmp = fabs( zMatElemNC(m,_zmat_decomp_row(i),_zmat_decomp_col(j)) ) ) > max ){
        max = tmp;
        *p = i;
        *q = j;
      }
#undef _zmat_decomp_row
#undef _zmat_decomp_col
  return max;
}

/* LU decomposition of column-long rectangular matrices with complete pivoting (destructive).
 * u has to be initially set for the original matrix. Rows of u are physically swapped, while columns
 * are not, so that u keeps the correspondence of columns to the original matrix. */
static int _zMatDecompLUCPDST_U(zMat l, zMat u, zIndex index_row, zIndex index_col, double tol)
{
  int i, j, k, p, q, cr, cj;
  double pivot, f, *ur, *ui;

  zMatZero( l );
  zIndexOrder( index_row, 0 );
  zIndexOrder( index_col, 0 );
  for( k=0; k<zMatRowSizeNC(u); k++ ){
    if( !( _zMatDecompLUCPPivot( u, NULL, index_col, k, &p, &q ) > tol ) ) break; /* remaining rows are dependent */
    if( p != k ){
      zMatSwapRowNC( u, k, p );
      zIndexSwapNC( index_row, k, p );
    }
    zIndexSwapNC( index_col, k, q );
    cr = zIndexElemNC(index_col,k);
    ur = zMatRowBufNC(u,k);
    zMatSetElemNC( l, zIndexElemNC(index_row,k), k, ( pivot = ur[cr] ) );
    ur[cr] = 1;
    for( j=k+1; j<zMatColSizeNC(u); j++ ) ur[zIndexElemNC(index_col,j)] /= pivot; /* U row */
    for( i=k+1; i<zMatRowSizeNC(u); i++ ){ /* L column and elimination of the other rows */
      ui = zMatRowBufNC(u,i);
      zMatSetElemNC( l, zIndexElemNC(index_row,i), k, ( f = ui[cr] ) );
      ui[cr] = 0;
      for( j=k+1; j<zMatColSizeNC(u); j++ ){
        cj = zIndexElemNC(index_col,j);
        ui[cj] -= f * ur[cj];
      }
    }
  }
  if( k < zMatRowSizeNC(u) ){ /* automatically adjust sizes of factorized matrices */
    zMatColResize( l, k ); /* to be column-full-rank */
    zMatRowResize( u, k ); /* to be row-full-rank */
  }
  return k; /* rank of the original matrix */
}

/* LU decomposition of row-long rectangular matrices with complete pivoting (destructive).
 * l has to be initially set for the original matrix. Columns of l are physically swapped, while rows
 * are not, so that l keeps the correspondence of rows to the original matrix. */
static int _zMatDecompLUCPDST_L(zMat l, zMat u, zIndex index_row, zIndex index_col, double tol)
{
  int i, j, k, p, q;
  double pivot, f, *lk, *li, *uk;

  zMatZero( u );
  zIndexOrder( index_row, 0 );
  zIndexOrder( index_col, 0 );
  for( k=0; k<zMatColSizeNC(l); k++ ){
    if( !( _zMatDecompLUCPPivot( l, index_row, NULL, k, &p, &q ) > tol ) ) break; /* remaining columns are dependent */
    zIndexSwapNC( index_row, k, p );
    if( q != k ){
      zMatSwapColNC( l, k, q );
      zIndexSwapNC( index_col, k, q );
    }
    lk = zMatRowBufNC(l,zIndexElemNC(index_row,k));
    uk = zMatRowBufNC(u,k);
    pivot = lk[k];
    uk[zIndexElemNC(index_col,k)] = 1;
    for( j=k+1; j<zMatColSizeNC(l); j++ ) uk[zIndexElemNC(index_col,j)] = lk[j] / pivot; /* U row */
    zRawVecZero( lk + k + 1, zMatColSizeNC(l) - k - 1 ); /* the pivot remains as a component of L */
    for( i=k+1; i<zMatRowSizeNC(l); i++ ){ /* elimination of the other rows, the k-th components of which remain as L */
      li = zMatRowBufNC(l,zIndexElemNC(index_row,i));
      f = li[k];
      for( j=k+1; j<zMatColSizeNC(l); j++ ) li[j] -= f * uk[zIndexElemNC(index_col,j)];
    }
  }
  if( k < zMatColSizeNC(l) ){ /* automatically adjust sizes of factorized matrices */
    zMatColResize( l, k ); /* to be column-full-rank */
    zMatRowResize( u, k ); /* to be row-full-rank */
  }
  return k; /* rank of the original matrix */
}

/* LU decomposition of a matrix with complete pivoting. */
int zMatDecompLUCP(const zMat m, zMat l, zMat u, zIndex index_row, zIndex index_col)
{
  double tol;

  if( zMatRowSize(m) != zIndexSizeNC(index_row) || zMatColSize(m) != zIndexSizeNC(index_col) ){
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
    return _zMatDecompLUCPDST_L( l, u, index_row, index_col, tol );
  }
  /* column-long rectangular matrix */
  if( !zMatSizeEqual( m, u ) || !zMatRowSizeEqual( m, l ) || !zMatIsSqr( l ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return -1;
  }
  zMatCopyNC( m, u );
  return _zMatDecompLUCPDST_U( l, u, index_row, index_col, tol );
}

/* LU decomposition with complete pivoting with an automatic matrix allocation. */
int zMatDecompLUCPAlloc(const zMat m, zMat *l, zMat *u, zIndex *index_row, zIndex *index_col)
{
  int minsize;

  minsize = zMin( zMatRowSizeNC(m), zMatColSizeNC(m) );
  *l = zMatAlloc( zMatRowSizeNC(m), minsize );
  *u = zMatAlloc( minsize, zMatColSizeNC(m) );
  *index_row = zIndexCreate( zMatRowSizeNC(m) );
  *index_col = zIndexCreate( zMatColSizeNC(m) );
  if( !*l || !*u || !*index_row || !*index_col ){
    zMatFree( *l );
    zMatFree( *u );
    zIndexFree( *index_row );
    zIndexFree( *index_col );
    *l = *u = NULL;
    *index_row = *index_col = NULL;
    return -1;
  }
  return zMatDecompLUCP( m, *l, *u, *index_row, *index_col );
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
      return -1; /* not positive semi-definite */
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
