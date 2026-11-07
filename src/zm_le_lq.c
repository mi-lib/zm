/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_lq - linear equation: LQ/QR decomposition.
 */

#include <zm/zm_le.h>

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

  if( zMatRowSizeNC(m) > zMatColSizeNC(m) ){ /* row-long rectangular matrix */
    if( !zMatSizeEqual( m, l ) || !zMatRowColSizeEqual( q, l ) || !zMatIsSqr( q ) ){
      ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
      return -1;
    }
    zMatCopyNC( m, l );
    rank = _zMatDecompLQDST_L( l, q );
  } else{ /* column-long rectangular matrix */
    if( !zMatSizeEqual( m, q ) || !zMatRowSizeEqual( m, l ) || !zMatIsSqr( l ) ){
      ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
      return -1;
    }
    zMatCopyNC( m, q );
    rank = _zMatDecompLQDST_Q( l, q );
  }
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
    return -1;
  }
  return zMatDecompLQ( m, *l, *q );
}

/* full-sized LQ decomposition based on Householder method. (destructive) */
static int _zMatDecompLQFullDST(zMat m, zMat q)
{
  double s, r, norm, d;
  double *u, *v;
  int i, j, colsize, rank;

  /* forward reflections: the k-th reflection is I - v v^T, where v (|v|^2=2) is stored in the k-th row of q. */
  zMatZero( q );
  for( rank=0, i=0; i<zMatRowSizeNC(m) && rank<zMatColSizeNC(m); i++ ){
    /* the norm of a row vector is preserved through the orthogonal transformations. */
    norm = zRawVecNorm( zMatRowBufNC(m,i), zMatColSizeNC(m) );
    u = &zMatElemNC(m,i,rank);
    colsize = zMatColSizeNC(m) - rank;
    if( _zMatDecompLQIsDependent( ( r = zRawVecNorm( u, colsize ) ), norm, ZM_LQ_DECOMP_TOL ) ){
      zRawVecZero( u, colsize ); /* linearly dependent on the preceding rows */
      continue;
    }
    s = -zSgn(*u) * r;
    v = &zMatElemNC(q,rank,rank);
    zRawVecCopy( u, v, colsize );
    *v -= s; /* reflection vector u - s e1 */
    zRawVecDivDRC( v, sqrt( s * ( s - *u ) ), colsize ); /* s ( s - u0 ) = r ( r + |u0| ) > 0 */
    for( j=i+1; j<zMatRowSizeNC(m); j++ ){
      d = zRawVecInnerProd( &zMatElemNC(m,j,rank), v, colsize );
      zRawVecCatDRC( &zMatElemNC(m,j,rank), -d, v, colsize );
    }
    *u = s;
    zRawVecZero( u + 1, colsize - 1 );
    rank++;
  }
  /* backward accumulation: q = H_{rank-1} ... H_0 is computed by multiplying H_{rank-1}, ..., H_0 from
   * the right side of the identity matrix, where only the trailing block of q is modified at each step. */
  for( i=rank; i<zMatColSizeNC(q); i++ )
    zMatSetElemNC( q, i, i, 1 );
  for( i=rank-1; i>=0; i-- ){
    v = &zMatElemNC(q,i,i);
    colsize = zMatColSizeNC(q) - i;
    for( j=i+1; j<zMatRowSizeNC(q); j++ ){ /* the i-th component of the j-th row is zero. */
      d = zRawVecInnerProd( &zMatElemNC(q,j,i), v, colsize );
      zRawVecCatDRC( &zMatElemNC(q,j,i), -d, v, colsize );
    }
    d = -*v; /* the i-th row is e_i - v_0 v */
    zRawVecMulDRC( v, d, colsize );
    *v += 1.0;
  }
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
