/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_mat_decomp_lq - matrix decomposition: LQ/QR decomposition.
 */

#include <zm/zm_mat_decomp.h>

/* internal functions for Householder reflections and scaling. */

/* exactly scale a matrix by a power of two to avoid overflow and underflow. */
double _zMatScalePow2(zMat m)
{
  double scale = 0, x;
  int i;

  for( i=0; i<zMatRowSizeNC(m); i++ ){
    x = zRawVecAbsMaxPow2( zMatRowBufNC(m,i), zMatColSizeNC(m) );
    if( !_zIsFinite( x ) ) return 1;
    if( x > scale ) scale = x;
  }
  if( scale == 0 ) return 1;
  zMatDivDRC( m, scale );
  return scale;
}

/* unit normal vector nu of the Householder reflection. */
bool _zMatHouseholderUnitNormal(const double *x, int n, double *nu, double *s)
{
  double scale, r = 0, d;
  int i;

  if( ( scale = zRawVecAbsMaxPow2( x, n ) ) == 0 ) return false;
  for( i=0; i<n; i++ ){ /* exactly scaled by a power of two to avoid overflow and underflow */
    nu[i] = x[i] / scale;
    r += nu[i] * nu[i];
  }
  r = sqrt( r );
  *s = x[0] >= 0 ? -r : r;
  nu[0] -= *s; /* |nu|^2 = 2 r ( r + |x0| ) */
  d = sqrt( 2 * r * ( r + fabs( x[0] / scale ) ) );
  for( i=0; i<n; i++ ) nu[i] /= d;
  *s *= scale;
  return true;
}

/* partially apply a Householder reflection to rows from the left side. */
void _zMatHouseholderL(zMat m, const double *nu, int n, int r0, int c0, int c1, double *w)
{
  int i, j;
  double a, *row;

  zRawVecZero( w + c0, c1 - c0 + 1 );
  for( i=0; i<n; i++ ){ /* w = nu^T m */
    row = zMatRowBufNC( m, r0+i );
    for( a=nu[i], j=c0; j<=c1; j++ ) w[j] += a * row[j];
  }
  for( i=0; i<n; i++ ){ /* m -= 2 nu w */
    row = zMatRowBufNC( m, r0+i );
    for( a=2*nu[i], j=c0; j<=c1; j++ ) row[j] -= a * w[j];
  }
}

/* partially apply a Householder reflection to columns from the right side. */
void _zMatHouseholderR(zMat m, const double *nu, int n, int c0, int r0, int r1)
{
  int i;
  double *p;

  for( i=r0; i<=r1; i++ ){
    p = &zMatElemNC(m,i,c0);
    zRawVecCatDRC( p, -2*zRawVecInnerProd( p, nu, n ), nu, n );
  }
}

/* the maximum norm of row vectors of a matrix. */
static double _zMatRowNormMax(const zMat m)
{
  int i;
  double norm, norm_max = 0;

  for( i=0; i<zMatRowSizeNC(m); i++ )
    if( ( norm = zRawVecNorm( zMatRowBufNC(m,i), zMatColSizeNC(m) ) ) > norm_max ) norm_max = norm;
  return norm_max;
}

/* workspace for the row pivoting of the Gram=Schmidt process: the squared norm of the residual of each
 * remaining row, which is downdated at each orthogonalization, and its base, namely, the squared norm exactly
 * computed last time, from which it has been downdated. they are stored in the last two columns of the space
 * given for the left factor, which are not used for the factor until the last two steps. ws is the whole space,
 * and row is the index of a row of the original matrix. */
#define _zMatDecompLQResSqrNorm(ws,row)     zMatElemNC(ws,row,zMatColSizeNC(ws)-1)
#define _zMatDecompLQResSqrNormBase(ws,row) zMatElemNC(ws,row,zMatColSizeNC(ws)-2)
/* the squared norm is computed again if it is downdated to be less than this ratio to the base, where the
 * relative error of the downdated one is in the order of the machine epsilon divided by this ratio. */
#define ZM_LQ_DOWNDATE_TOL ( 1.0e-8 )

/* LQ decomposition of column-long rectangular matrices based on modified Gram=Schmidt's method with row
 * pivoting. (destructive)
 * if reorth is true, the pivot row is orthogonalized again against the preceding orthonormal rows before
 * normalized, so that the orthogonality of q is kept within the order of the machine epsilon regardless
 * of the condition number of the matrix (twice is enough). the rank is detected in the same way in both
 * cases, namely, with the norms of the residuals before the reorthogonalization.
 * the norms of the residuals for pivoting are downdated instead of computed at each step while three or
 * more rows remain, where the workspace is borrowed from l. the norm of the pivot row is exactly computed. */
static int _zMatDecompLQDST_Q(zMat l, zMat q, zIndex index, bool reorth) /* q has to be initially set for the original matrix. */
{
  int i, p, rank, row;
  double *qk, *qi, r, norm, tol, s;
  zMat ws = l; /* workspace for pivoting borrowed from l */
  bool use_ws;

  zMatZero( l );
  zIndexOrder( index, 0 );
  tol = zTOL * _zMatRowNormMax( q );
  if( zMatRowSizeNC(q) > 2 )
    for( i=0; i<zMatRowSizeNC(q); i++ )
      _zMatDecompLQResSqrNormBase(ws,i) = _zMatDecompLQResSqrNorm(ws,i) = zRawVecSqrNorm( zMatRowBufNC(q,i), zMatColSizeNC(q) );
  for( rank=0; rank<zMatRowSizeNC(q); rank++ ){
    /* pivoting: the row with the largest residual is orthonormalized next (rows are physically swapped). */
    if( ( use_ws = zMatRowSizeNC(q) - rank > 2 ) ){
      for( norm=-1, p=i=rank; i<zMatRowSizeNC(q); i++ ){
        row = zIndexElemNC(index,i);
        if( _zMatDecompLQResSqrNorm(ws,row) < ZM_LQ_DOWNDATE_TOL * _zMatDecompLQResSqrNormBase(ws,row) ) /* cancellation */
          _zMatDecompLQResSqrNormBase(ws,row) = _zMatDecompLQResSqrNorm(ws,row) = zRawVecSqrNorm( zMatRowBufNC(q,i), zMatColSizeNC(q) );
        if( _zMatDecompLQResSqrNorm(ws,row) > norm ){
          norm = _zMatDecompLQResSqrNorm(ws,row);
          p = i;
        }
      }
      r = zRawVecNorm( zMatRowBufNC(q,p), zMatColSizeNC(q) );
    } else{ /* the workspace is not available in the last two steps. */
      for( r=0, p=i=rank; i<zMatRowSizeNC(q); i++ )
        if( ( norm = zRawVecNorm( zMatRowBufNC(q,i), zMatColSizeNC(q) ) ) > r ){
          r = norm;
          p = i;
        }
    }
    if( !( r > tol ) ) break; /* the residuals of all the remaining rows are negligible. */
    if( p != rank ){
      zMatSwapRowNC( q, rank, p );
      zIndexSwapNC( index, rank, p );
    }
    if( use_ws ) /* the workspace of the pivot row is released to be a part of the upper triangle of l. */
      _zMatDecompLQResSqrNormBase(ws,zIndexElemNC(index,rank)) = _zMatDecompLQResSqrNorm(ws,zIndexElemNC(index,rank)) = 0;
    qk = zMatRowBufNC(q,rank);
    if( reorth ){ /* the coefficients of the reorthogonalization are added to l */
      for( i=0; i<rank; i++ ){
        qi = zMatRowBufNC(q,i);
        s = zRawVecInnerProd( qk, qi, zMatColSizeNC(q) );
        zRawVecCatDRC( qk, -s, qi, zMatColSizeNC(q) );
        zMatElemNC(l,zIndexElemNC(index,rank),i) += s;
      }
      if( !( ( r = zRawVecNorm( qk, zMatColSizeNC(q) ) ) > 0 ) ) break;
    }
    zRawVecDivDRC( qk, r, zMatColSizeNC(q) );
    zMatSetElemNC( l, zIndexElemNC(index,rank), rank, r );
    for( i=rank+1; i<zMatRowSizeNC(q); i++ ){ /* orthogonalization of the remaining rows */
      qi = zMatRowBufNC(q,i);
      r = zRawVecInnerProd( qi, qk, zMatColSizeNC(q) );
      zRawVecCatDRC( qi, -r, qk, zMatColSizeNC(q) );
      zMatSetElemNC( l, zIndexElemNC(index,i), rank, r );
      if( zMatRowSizeNC(q) - rank > 3 ) /* downdated for the next step. */
        _zMatDecompLQResSqrNorm(ws,zIndexElemNC(index,i)) -= r*r;
      else
      if( zMatRowSizeNC(q) - rank == 3 ) /* the workspace is released since it is not used anymore. */
        _zMatDecompLQResSqrNormBase(ws,zIndexElemNC(index,i)) = _zMatDecompLQResSqrNorm(ws,zIndexElemNC(index,i)) = 0;
    }
  }
  if( rank < zMatRowSizeNC(q) - 2 ) /* the workspace of the remaining rows is released. */
    for( i=rank; i<zMatRowSizeNC(q); i++ )
      _zMatDecompLQResSqrNormBase(ws,zIndexElemNC(index,i)) = _zMatDecompLQResSqrNorm(ws,zIndexElemNC(index,i)) = 0;
  if( rank < zMatRowSizeNC(q) ){ /* automatically adjust sizes of factorized matrices */
    zMatColResize( l, rank ); /* to be column-full-rank */
    zMatRowResize( q, rank ); /* to be row-full-rank */
  }
  return rank;
}

/* full-sized LQ decomposition based on Householder method with row pivoting. (destructive) */
static int _zMatDecompLQFullDST(zMat m, zMat q, zIndex index) /* m has to be initially set for the original matrix. */
{
  double scale, s, norm, norm_max, tol;
  double *u, *nu;
  int i, p, colsize, rank;

  scale = _zMatScalePow2( m );
  tol = zTOL * _zMatRowNormMax( m );
  zIndexOrder( index, 0 );
  /* forward reflections: the k-th reflection is I - 2 nu nu^T, where the unit vector nu is stored in
   * the k-th row of q. Rows of m are not physically swapped to keep the correspondence to the original. */
  zMatZero( q );
  for( rank=0; rank<zMatRowSizeNC(m) && rank<zMatColSizeNC(m); rank++ ){
    colsize = zMatColSizeNC(m) - rank;
    /* pivoting: the row with the largest remaining components is reflected next. */
    for( norm_max=0, p=i=rank; i<zMatRowSizeNC(m); i++ )
      if( ( norm = zRawVecNorm( &zMatElemNC(m,zIndexElemNC(index,i),rank), colsize ) ) > norm_max ){
        norm_max = norm;
        p = i;
      }
    if( !( norm_max > tol ) ) break; /* the remaining components of all the remaining rows are negligible. */
    zIndexSwapNC( index, rank, p );
    u = &zMatElemNC(m,zIndexElemNC(index,rank),rank);
    nu = &zMatElemNC(q,rank,rank);
    _zMatHouseholderUnitNormal( u, colsize, nu, &s );
    for( i=rank+1; i<zMatRowSizeNC(m); i++ )
      _zMatHouseholderR( m, nu, colsize, rank, zIndexElemNC(index,i), zIndexElemNC(index,i) );
    *u = s;
    zRawVecZero( u + 1, colsize - 1 );
  }
  for( i=rank; i<zMatRowSizeNC(m); i++ ) /* the remaining components of linearly dependent rows are discarded. */
    zRawVecZero( &zMatElemNC(m,zIndexElemNC(index,i),rank), zMatColSizeNC(m) - rank );
  /* backward accumulation: q = H_{rank-1} ... H_0 is computed by multiplying H_{rank-1}, ..., H_0 from
   * the right side of the identity matrix, where only the trailing block of q is modified at each step. */
  for( i=rank; i<zMatColSizeNC(q); i++ )
    zMatSetElemNC( q, i, i, 1 );
  for( i=rank-1; i>=0; i-- ){
    nu = &zMatElemNC(q,i,i);
    colsize = zMatColSizeNC(q) - i;
    /* the i-th component of the j-th row (j>i) is zero. */
    if( i < zMatRowSizeNC(q) - 1 )
      _zMatHouseholderR( q, nu, colsize, i, i+1, zMatRowSizeNC(q)-1 );
    s = -2 * *nu; /* the i-th row is e_i - 2 nu_0 nu */
    zRawVecMulDRC( nu, s, colsize );
    *nu += 1.0;
  }
  zMatMulDRC( m, scale );
  return rank;
}

/* LQ decomposition, where the rows of q of a column-long matrix are reorthogonalized if reorth is true. */
static int _zMatDecompLQ(const zMat m, zMat l, zMat q, zIndex index, bool reorth)
{
  int rank;
  double scale;

  if( zMatRowSizeNC(m) != zIndexSizeNC(index) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return -1;
  }
  if( zMatRowSizeNC(m) > zMatColSizeNC(m) ){ /* row-long rectangular matrix */
    if( !zMatSizeEqual( m, l ) || !zMatRowColSizeEqual( q, l ) || !zMatIsSqr( q ) ){
      ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
      return -1;
    }
    zMatCopyNC( m, l );
    rank = _zMatDecompLQFullDST( l, q, index );
    if( rank < zMatColSizeNC(l) ){ /* automatically adjust sizes of factorized matrices */
      zMatColResize( l, rank ); /* to be column-full-rank */
      zMatRowResize( q, rank ); /* to be row-full-rank */
    }
  } else{ /* column-long rectangular matrix */
    if( !zMatSizeEqual( m, q ) || !zMatRowSizeEqual( m, l ) || !zMatIsSqr( l ) ){
      ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
      return -1;
    }
    zMatCopyNC( m, q );
    scale = _zMatScalePow2( q );
    rank = _zMatDecompLQDST_Q( l, q, index, reorth );
    zMatMulDRC( l, scale ); /* the scaling by a power of two to avoid overflow and underflow is undone */
  }
  return rank;
}

/* LQ decomposition. */
int zMatDecompLQ(const zMat m, zMat l, zMat q, zIndex index)
{
  return _zMatDecompLQ( m, l, q, index, false );
}

/* LQ decomposition with reorthogonalization. */
int zMatDecompLQReorth(const zMat m, zMat l, zMat q, zIndex index)
{
  return _zMatDecompLQ( m, l, q, index, true );
}

/* LQ decomposition with an automatic matrix allocation. */
int zMatDecompLQAlloc(const zMat m, zMat *l, zMat *q, zIndex *index)
{
  int minsize;

  minsize = zMin( zMatRowSizeNC(m), zMatColSizeNC(m) );
  *l = zMatAlloc( zMatRowSizeNC(m), minsize );
  *q = zMatAlloc( minsize, zMatColSizeNC(m) );
  *index = zIndexCreate( zMatRowSizeNC(m) );
  if( !*l || !*q || !*index ){
    zMatFree( *l );
    zMatFree( *q );
    zIndexFree( *index );
    *l = *q = NULL;
    *index = NULL;
    return -1;
  }
  return zMatDecompLQ( m, *l, *q, *index );
}

/* full-sized LQ decomposition based on Householder method. */
int zMatDecompLQFull(const zMat m, zMat l, zMat q, zIndex index)
{
  if( !zMatSizeEqual( m, l ) || !zMatRowColSizeEqual( q, m ) || !zMatIsSqr( q ) ||
      zMatRowSizeNC(m) != zIndexSizeNC(index) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return -1;
  }
  zMatCopyNC( m, l );
  return _zMatDecompLQFullDST( l, q, index );
}

/* full-sized LQ decomposition with an automatic matrix allocation based on Householder method. */
int zMatDecompLQFullAlloc(const zMat m, zMat *l, zMat *q, zIndex *index)
{
  *l = zMatAlloc( zMatRowSizeNC(m), zMatColSizeNC(m) );
  *q = zMatAllocSqr( zMatColSizeNC(m) );
  *index = zIndexCreate( zMatRowSizeNC(m) );
  if( !*l || !*q || !*index ){
    zMatFree( *l );
    zMatFree( *q );
    zIndexFree( *index );
    *l = *q = NULL;
    *index = NULL;
    return -1;
  }
  return zMatDecompLQFull( m, *l, *q, *index );
}

/* LQ decomposition based on Householder method. */
int zMatDecompLQNull(const zMat m, zMat l, zMat q, zMat qnull, zIndex index)
{
  int rank;

  if( ( rank = zMatDecompLQFull( m, l, q, index ) ) >= 0 && qnull ){
    zMatSetColSizeNC( qnull, zMatColSizeNC(m) - rank );
    zMatTGetNC( q, rank, 0, qnull );
    zMatSetColSizeNC( l, rank );
    zMatSetRowSizeNC( q, rank );
  }
  return rank;
}

/* QR decomposition. */
int zMatDecompQR(const zMat m, zMat q, zMat r, zIndex index)
{
  zMat mcp, qcp, rcp;
  int rank = -1;

  mcp = zMatAlloc( zMatColSizeNC(m), zMatRowSizeNC(m) );
  qcp = zMatAlloc( zMatColSizeNC(q), zMatRowSizeNC(q) );
  rcp = zMatAlloc( zMatColSizeNC(r), zMatRowSizeNC(r) );
  if( !mcp || !qcp || !rcp ) goto TERMINATE;
  zMatTNC( m, mcp );
  if( ( rank = zMatDecompLQ( mcp, rcp, qcp, index ) ) < 0 ) goto TERMINATE;
  zMatColResize( q, rank ); /* to be column-full-rank */
  zMatRowResize( r, rank ); /* to be row-full-rank */
  zMatTNC( qcp, q );
  zMatTNC( rcp, r );

 TERMINATE:
  zMatFreeAtOnce( 3, mcp, qcp, rcp );
  return rank;
}
