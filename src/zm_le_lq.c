/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_lq - linear equation: LQ/QR decomposition.
 */

#include <zm/zm_le.h>

/* internal functions for Householder reflections and scaling. */

/* exactly scale a matrix by a power of two to avoid overflow and underflow. */
double _zMatScalePow2(zMat m)
{
  double scale;
  int e;

  if( ( scale = zDataAbsMax( zMatBuf(m), zMatRowSizeNC(m)*zMatColSizeNC(m), NULL ) ) == 0 ) return 1;
  frexp( scale, &e );
  scale = ldexp( 1.0, e - 1 );
  zMatDivDRC( m, scale );
  return scale;
}

/* unit normal vector nu of the Householder reflection. */
bool _zMatHouseholderUnitNormal(const double *x, int n, double *nu, double *s)
{
  double scale = 0, r = 0, d, scale_tmp;
  int i;

  for( i=0; i<n; i++ )
    if( ( scale_tmp = fabs( x[i] ) ) > scale ) scale = scale_tmp;
  if( scale == 0 ) return false;
  for( i=0; i<n; i++ ){ /* scaled to avoid overflow and underflow */
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

/* check if the residual of a row vector after orthogonalization is negligible relative to its original norm. */
static bool _zMatDecompLQIsDependent(double residual, double norm, double tol)
{
  return residual <= tol * norm;
}

/* LQ decomposition of column-long rectangular matrices based on modified Gram=Schmidt's method. (destructive) */
static int _zMatDecompLQDST_Q(zMat l, zMat q) /* q has to be initially set for the original matrix. */
{
  int i, j, rank;
  double *qp, r, norm;

  zMatZero( l );
  for( rank=0, i=0; i<zMatRowSizeNC(q); i++ ){
    qp = zMatRowBufNC(q,i);
    norm = zRawVecNorm( qp, zMatColSizeNC(q) );
    for( j=0; j<rank; j++ ){
      r = zRawVecInnerProd( qp, zMatRowBufNC(q,j), zMatColSizeNC(q) );
      zRawVecCatDRC( qp, -r, zMatRowBufNC(q,j), zMatColSizeNC(q) );
      zMatSetElemNC( l, i, j, r );
    }
    if( _zMatDecompLQIsDependent( ( r = zRawVecNorm(qp,zMatColSizeNC(q)) ), norm, ZM_LQ_DECOMP_TOL ) ) continue;
    zRawVecDiv( qp, r, zMatRowBufNC(q,rank), zMatColSizeNC(q) );
    zMatSetElemNC( l, i, rank++, r );
  }
  if( rank == 0 ) ZRUNERROR( ZM_ERR_MAT_CANNOTDECOMPOSEZEROMAT );
  if( rank < zMatRowSizeNC(q) ){ /* automatically adjust sizes of factorized matrices */
    zMatColResize( l, rank ); /* to be column-full-rank */
    zMatRowResize( q, rank ); /* to be row-full-rank */
  }
  return rank;
}

/* LQ decomposition of row-long rectangular matrices based on modified Gram=Schmidt's method. (destructive) */
static int _zMatDecompLQDST_L(zMat l, zMat q) /* l has to be initially set for the original matrix. */
{
  int i, j, rank = 0;
  double r0, r, norm, *ws = NULL, *v;

  for( i=0; i<zMatRowSizeNC(l); i++ ){
    if( ( r0 = zRawVecNorm( zMatRowBufNC(l,i), zMatColSizeNC(l) ) ) > 0 ){
      ws = zMatRowBufNC(l,i);
      break;
    }
    zRawVecZero( zMatRowBufNC(l,i), zMatColSizeNC(l) ); /* leading null rows */
  }
  if( !ws ){
    ZRUNERROR( ZM_ERR_MAT_CANNOTDECOMPOSEZEROMAT );
    goto ZMAT_DECOMP_LQ_L_RESIZE;
  }
  zRawVecDiv( ws, r0, zMatRowBufNC(q,0), zMatColSizeNC(l) );
  for( rank=1, i++; i<zMatRowSizeNC(l); i++ ){
    /* the row vector is orthogonalized in the next row of q, or in ws if q is already full. */
    v = rank < zMatColSizeNC(l) ? zMatRowBufNC(q,rank) : ws;
    zRawVecCopy( zMatRowBufNC(l,i), v, zMatColSizeNC(l) );
    norm = zRawVecNorm( v, zMatColSizeNC(l) );
    for( j=0; j<rank; j++ ){
      r = zRawVecInnerProd( v, zMatRowBufNC(q,j), zMatColSizeNC(q) );
      zRawVecCatDRC( v, -r, zMatRowBufNC(q,j), zMatColSizeNC(q) );
      zMatSetElemNC( l, i, j, r );
    }
    if( rank >= zMatColSizeNC(l) ) continue;
    zRawVecZero( zMatRowBufNC(l,i) + rank, zMatColSizeNC(l) - rank );
    if( _zMatDecompLQIsDependent( ( r = zRawVecNorm( v, zMatColSizeNC(q) ) ), norm, ZM_LQ_DECOMP_TOL ) ) continue;
    zRawVecDivDRC( v, r, zMatColSizeNC(q) );
    zMatSetElemNC( l, i, rank++, r );
  }
  *ws = r0;
  zRawVecZero( ws + 1, zMatColSizeNC(l) - 1 );
 ZMAT_DECOMP_LQ_L_RESIZE:
  if( rank < zMatRowSizeNC(q) ){ /* automatically adjust sizes of factorized matrices */
    zMatColResize( l, rank ); /* to be column-full-rank */
    zMatRowResize( q, rank ); /* to be row-full-rank */
  }
  return rank;
}

/* LQ decomposition based on modified Gram=Schmidt's method. */
int zMatDecompLQ(const zMat m, zMat l, zMat q)
{
  int rank;
  double scale;

  if( zMatRowSizeNC(m) > zMatColSizeNC(m) ){ /* row-long rectangular matrix */
    if( !zMatSizeEqual( m, l ) || !zMatRowColSizeEqual( q, l ) || !zMatIsSqr( q ) ){
      ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
      return -1;
    }
    zMatCopyNC( m, l );
    scale = _zMatScalePow2( l );
    rank = _zMatDecompLQDST_L( l, q );
  } else{ /* column-long rectangular matrix */
    if( !zMatSizeEqual( m, q ) || !zMatRowSizeEqual( m, l ) || !zMatIsSqr( l ) ){
      ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
      return -1;
    }
    zMatCopyNC( m, q );
    scale = _zMatScalePow2( q );
    rank = _zMatDecompLQDST_Q( l, q );
  }
  zMatMulDRC( l, scale ); /* the scaling by a power of two to avoid overflow and underflow is undone */
  return rank;
}

/* LQ decomposition with an automatic matrix allocation. */
int zMatDecompLQAlloc(const zMat m, zMat *l, zMat *q)
{
  int minsize;

  minsize = zMin( zMatRowSizeNC(m), zMatColSizeNC(m) );
  *l = zMatAlloc( zMatRowSizeNC(m), minsize );
  *q = zMatAlloc( minsize, zMatColSizeNC(m) );
  if( !*l || !*q ){
    zMatFree( *l );
    zMatFree( *q );
    *l = *q = NULL;
    return -1;
  }
  return zMatDecompLQ( m, *l, *q );
}

/* full-sized LQ decomposition based on Householder method. (destructive) */
static int _zMatDecompLQFullDST(zMat m, zMat q)
{
  double scale, s, norm;
  double *u, *nu;
  int i, colsize, rank;

  scale = _zMatScalePow2( m );
  /* forward reflections: the k-th reflection is I - 2 nu nu^T, where the unit vector nu is stored in
   * the k-th row of q. */
  zMatZero( q );
  for( rank=0, i=0; i<zMatRowSizeNC(m) && rank<zMatColSizeNC(m); i++ ){
    /* the norm of a row vector is preserved through the orthogonal transformations. */
    norm = zRawVecNorm( zMatRowBufNC(m,i), zMatColSizeNC(m) );
    u = &zMatElemNC(m,i,rank);
    colsize = zMatColSizeNC(m) - rank;
    nu = &zMatElemNC(q,rank,rank);
    if( !_zMatHouseholderUnitNormal( u, colsize, nu, &s ) ||
        _zMatDecompLQIsDependent( fabs( s ), norm, ZM_LQ_DECOMP_TOL ) ){
      zRawVecZero( nu, colsize );
      zRawVecZero( u, colsize ); /* linearly dependent on the preceding rows */
      continue;
    }
    if( i < zMatRowSizeNC(m) - 1 )
      _zMatHouseholderR( m, nu, colsize, rank, i+1, zMatRowSizeNC(m)-1 );
    *u = s;
    zRawVecZero( u + 1, colsize - 1 );
    rank++;
  }
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

/* full-sized LQ decomposition based on Householder method. */
int zMatDecompLQFull(const zMat m, zMat l, zMat q)
{
  if( !zMatSizeEqual( m, l ) || !zMatRowColSizeEqual( q, m ) || !zMatIsSqr( q ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return -1;
  }
  zMatCopyNC( m, l );
  return _zMatDecompLQFullDST( l, q );
}

/* full-sized LQ decomposition with an automatic matrix allocation based on Householder method. */
int zMatDecompLQFullAlloc(const zMat m, zMat *l, zMat *q)
{
  *l = zMatAlloc( zMatRowSizeNC(m), zMatColSizeNC(m) );
  *q = zMatAllocSqr( zMatColSizeNC(m) );
  if( !*l || !*q ){
    zMatFree( *l );
    zMatFree( *q );
    *l = *q = NULL;
    return -1;
  }
  return zMatDecompLQFull( m, *l, *q );
}

/* LQ decomposition based on Householder method. */
int zMatDecompLQNull(const zMat m, zMat l, zMat q, zMat qnull)
{
  int rank;

  if( ( rank = zMatDecompLQFull( m, l, q ) ) >= 0 && qnull ){
    zMatSetColSizeNC( qnull, zMatColSizeNC(m) - rank );
    zMatTGetNC( q, rank, 0, qnull );
    zMatSetColSizeNC( l, rank );
    zMatSetRowSizeNC( q, rank );
  }
  return rank;
}

/* QR decomposition. */
int zMatDecompQR(const zMat m, zMat q, zMat r)
{
  zMat mcp, qcp, rcp;
  int rank = -1;

  mcp = zMatAlloc( zMatColSizeNC(m), zMatRowSizeNC(m) );
  qcp = zMatAlloc( zMatColSizeNC(q), zMatRowSizeNC(q) );
  rcp = zMatAlloc( zMatColSizeNC(r), zMatRowSizeNC(r) );
  if( !mcp || !qcp || !rcp ){
    ZALLOCERROR();
    goto TERMINATE;
  }
  zMatTNC( m, mcp );
  rank = zMatDecompLQ( mcp, rcp, qcp );
  zMatTNC( qcp, q );
  zMatTNC( rcp, r );

 TERMINATE:
  zMatFreeAtOnce( 3, mcp, qcp, rcp );
  return rank;
}
