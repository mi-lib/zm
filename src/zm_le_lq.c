/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_lq - linear equation: LQ/QR decomposition.
 */

#include <zm/zm_le.h>

/* LQ decomposition based on Gram=Schmidt's method. (destructive) */
int zMatDecompLQDST(zMat m, zMat l, zMat q)
{
  int i, j, rank;
  double *mp, r;

  zMatZero( l );
  for( rank=0, i=0; i<zMatRowSizeNC(m); i++ ){
    mp = zMatRowBuf(m,i);
    for( j=0; j<rank; j++ ){
      r = zRawVecInnerProd( mp, zMatRowBuf(q,j), zMatColSizeNC(m) );
      zRawVecCatDRC( mp, -r, zMatRowBuf(q,j), zMatColSizeNC(q) );
      zMatSetElemNC( l, i, j, r );
    }
    if( zIsTol( ( r = zRawVecNorm(mp,zMatColSizeNC(m)) ), ZM_LQ_DECOMP_GRAMSCHMIDT_TOL ) ) continue;
    zRawVecDiv( mp, r, zMatRowBuf(q,rank), zMatColSizeNC(m) );
    zMatSetElemNC( l, i, rank, r );
    if( rank < zMatColSizeNC(m) ) rank++;
  }
  if( rank < zMatRowSizeNC(m) ){ /* automatically adjust sizes of factorized matrices */
    zMatColResize( l, rank ); /* to be column-full-rank */
    zMatRowResize( q, rank ); /* to be row-full-rank */
  }
  return rank;
}

/* LQ decomposition based on Gram=Schmidt's method. */
int zMatDecompLQ(const zMat m, zMat l, zMat q)
{
  zMat mcp;
  int rank;

  if( !( mcp = zMatClone( m ) ) ){
    ZALLOCERROR();
    return -1;
  }
  rank = zMatDecompLQDST( mcp, l, q );
  zMatFree( mcp );
  return rank;
}

/* LQ decomposition with an automatic matrix allocation. */
int zMatDecompLQAlloc(const zMat m, zMat *l, zMat *q)
{
  *l = zMatAllocSqr( zMatRowSizeNC(m) );
  *q = zMatAlloc( zMatRowSizeNC(m), zMatColSizeNC(m) );
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
  if( !zMatSizeEqual( m, l ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return -1;
  }
  zMatCopyNC( m, l );
  return _zMatDecompLQFullDST( l, q );
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
  rank = zMatDecompLQDST( mcp, rcp, qcp );
  zMatTNC( qcp, q );
  zMatTNC( rcp, r );

 TERMINATE:
  zMatFreeAtOnce( 3, mcp, qcp, rcp );
  return rank;
}
