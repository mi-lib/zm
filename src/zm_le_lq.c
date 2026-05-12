/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_lq - linear equation: LQ/QR decomposition.
 */

#include <zm/zm_le.h>

/* LQ decomposition of column-long rectangular matrices based on Gram=Schmidt's method. (destructive) */
static int _zMatDecompLQDST_Q(zMat l, zMat q) /* the original matrix has to be stored in q. */
{
  int i, j, rank;
  double *qp, r;

  zMatZero( l );
  for( rank=0, i=0; i<zMatRowSizeNC(q); i++ ){
    qp = zMatRowBufNC(q,i);
    for( j=0; j<rank; j++ ){
      r = zRawVecInnerProd( qp, zMatRowBufNC(q,j), zMatColSizeNC(q) );
      zRawVecCatDRC( qp, -r, zMatRowBufNC(q,j), zMatColSizeNC(q) );
      zMatSetElemNC( l, i, j, r );
    }
    if( zIsTol( ( r = zRawVecNorm(qp,zMatColSizeNC(q)) ) / zDataAbsMax(zMatRowBufNC(l,i),rank-1,NULL), ZM_LQ_DECOMP_GRAMSCHMIDT_TOL ) ) continue;
    zRawVecDivDRC( qp, r, zMatColSizeNC(q) );
    zMatSetElemNC( l, i, rank, r );
    if( rank < zMatColSizeNC(q) ) rank++;
  }
  if( rank < zMatRowSizeNC(q) ){ /* automatically adjust sizes of factorized matrices */
    zMatColResize( l, rank ); /* to be column-full-rank */
    zMatRowResize( q, rank ); /* to be row-full-rank */
  }
  return rank;
}

/* LQ decomposition of row-long rectangular matrices based on Gram=Schmidt's method. (destructive) */
static int _zMatDecompLQDST_L(zMat l, zMat q) /* the original matrix has to be stored in l. */
{
  int i, j, rank;
  double r0, r, *ws = NULL;

  for( i=0; i<zMatRowSizeNC(l); i++ ){
    if( !zIsTol( ( r0 = zRawVecNorm( zMatRowBufNC(l,i), zMatColSizeNC(l) ) ), ZM_LQ_DECOMP_GRAMSCHMIDT_TOL ) ){
      ws = zMatRowBufNC(l,i);
      break;
    }
  }
  if( !ws ){
    ZRUNERROR( ZM_ERR_MAT_CANNOTDECOMPOSEZEROMAT );
    return 0;
  }
  zRawVecDiv( ws, r0, zMatRowBufNC(q,0), zMatColSizeNC(l) ); /* temporary workspace */
  for( rank=1, i++; i<zMatRowSizeNC(l); i++ ){
    zRawVecCopy( zMatRowBufNC(l,i), ws, zMatColSizeNC(l) ); /* temporarily store the raw vector */
    for( j=0; j<rank; j++ )
      zMatSetElemNC( l, i, j, zRawVecInnerProd( ws, zMatRowBufNC(q,j), zMatColSizeNC(q) ) );
    if( rank >= zMatColSizeNC(l) ) continue;
    zRawVecCopy( ws, zMatRowBufNC(q,rank), zMatColSizeNC(q) ); /* temporarily store the raw vector */
    for( j=0; j<rank; j++ )
      zRawVecCatDRC( zMatRowBufNC(q,rank), -zMatElemNC(l,i,j), zMatRowBufNC(q,j), zMatColSizeNC(q) );
    zMatSetElemNC( l, i, rank, ( r = zRawVecNorm( zMatRowBufNC(q,rank), zMatColSizeNC(q) ) ) );
    zRawVecZero( zMatRowBufNC(l,i) + rank + 1, zMatColSizeNC(l) - rank - 1 );
    if( zIsTol( r / zDataAbsMax(zMatRowBufNC(l,i),rank-1,NULL), ZM_LQ_DECOMP_GRAMSCHMIDT_TOL ) ) continue;
    zRawVecDivDRC( zMatRowBufNC(q,rank), r, zMatColSizeNC(q) );
    rank++;
  }
  *ws = r0;
  zRawVecZero( ws + 1, zMatColSizeNC(l) - 1 );
  if( rank < zMatRowSizeNC(q) ){ /* automatically adjust sizes of factorized matrices */
    zMatColResize( l, rank ); /* to be column-full-rank */
    zMatRowResize( q, rank ); /* to be row-full-rank */
  }
  return rank;
}

/* LQ decomposition based on Gram=Schmidt's method. */
int zMatDecompLQ(const zMat m, zMat l, zMat q)
{
  int rank;

  if( zMatRowSizeNC(m) > zMatColSizeNC(m) ){ /* column-long rectangular matrix */
    if( !zMatSizeEqual( m, l ) || !zMatRowColSizeEqual( q, l ) || !zMatIsSqr( q ) ){
      ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
      return -1;
    }
    zMatCopyNC( m, l );
    rank = _zMatDecompLQDST_L( l, q );
  } else{ /* row-long rectangular matrix */
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
  double s, ds, norm_inv, reflection;
  double *u;
  int i, j, size, colsize, rank;

  zMatIdentNC( q );
  size = zMatMinSize( m );
  for( rank=0, i=0; i<size; i++ ){
    u = &zMatElemNC(m,i,i);
    colsize = zMatColSizeNC(m) - i;
    if( zIsTiny( ( s = -zSgn( zMatElemNC(m,i,i) ) * zRawVecNorm( u, colsize ) ) ) ||
        zIsTiny( ( ds = s - zMatElemNC(m,i,i) ) ) ) continue;
    norm_inv = 1.0 / ( s * ds );
    *u = -ds; /* to use u temporarily as a reflection vector */
    for( j=0; j<zMatColSizeNC(q); j++ ){
      reflection = -norm_inv * zRawMatColInnerProd( zMatRowBufNC(q,i), zMatColCapacity(q), u, zMatRowSizeNC(q)-i, zMatColSizeNC(q), j );
      zRawMatColCatDRC( zMatRowBufNC(q,i), zMatColCapacity(q), reflection, u, zMatRowSizeNC(q)-i, zMatColSizeNC(q), j );
    }
    for( j=zMatRowSizeNC(m)-1; j>i; j-- ){
      reflection = -norm_inv * zRawVecInnerProd( &zMatElemNC(m,j,i), u, colsize );
      zRawVecCatDRC( &zMatElemNC(m,j,i), reflection, u, colsize );
    }
    *u = s;
    zRawVecZero( u + 1, colsize - 1 );
    rank++;
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

  if( ( rank = zMatDecompLQFull( m, l, q ) ) > 0 && qnull ){
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
