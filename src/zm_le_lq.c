/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_lq - linear equation: LQ/QR decomposition.
 */

#include <zm/zm_le.h>

/* LQ decomposition based on Gram=Schmidt's method. (destructive) */
int zMatDecompLQDST(zMat m, zMat l, zMat q, zIndex idx)
{
  int i, j, rank;
  double *mp, r;

  zMatZero( l );
  zIndexOrder( idx, 0 );
  for( rank=0, i=0; i<zMatRowSizeNC(m); i++ ){
    mp = zMatRowBuf(m,i);
    for( j=0; j<rank; j++ ){
      r = zRawVecInnerProd( mp, zMatRowBuf(q,j), zMatColSizeNC(m) );
      zRawVecCatDRC( mp, -r, zMatRowBuf(q,j), zMatColSizeNC(q) );
      zMatElemNC(l,i,j) = r;
    }
    if( zIsTiny( r = zRawVecNorm(mp,zMatColSizeNC(m)) ) ){
      zIndexMove( idx, rank, zArraySize(idx)-1 );
      continue;
    }
    zRawVecDiv( mp, r, zMatRowBuf(q,rank), zMatColSizeNC(m) );
    zMatElemNC(l,i,rank) = r;
    if( rank < zMatColSizeNC(m) ) rank++;
  }
  return rank;
}

/* LQ decomposition based on Gram=Schmidt's method. */
int zMatDecompLQ(zMat m, zMat l, zMat q, zIndex idx)
{
  zMat mcp;
  int rank;

  if( !( mcp = zMatClone( m ) ) ){
    ZALLOCERROR();
    return -1;
  }
  rank = zMatDecompLQDST( mcp, l, q, idx );
  zMatFree( mcp );
  return rank;
}

/* LQ decomposition and regression. */
int zMatDecompLQReg(zMat m, zMat l, zMat q, zIndex idx)
{
  int rank;

  if( ( rank = zMatDecompLQ( m, l, q, idx ) ) < 0 ) return -1;
  if( rank < (int)zMatRowSizeNC(m) ){
    zMatColReg( l, rank );
    zMatRowReg( q, rank );
  }
  return rank;
}

/* LQ decomposition with an automatic matrix allocation and resize. */
int zMatDecompLQAlloc(zMat m, zMat *l, zMat *q, zIndex *idx)
{
  *l = zMatAllocSqr( zMatRowSizeNC(m) );
  *q = zMatAlloc( zMatRowSizeNC(m), zMatColSizeNC(m) );
  *idx = zIndexCreate( zMatRowSizeNC(m) );
  if( !*l || !*q || !*idx ){
    zMatFree( *l );
    zMatFree( *q );
    zIndexFree( *idx );
    return -1;
  }
  return zMatDecompLQReg( m, *l, *q, *idx );
}

/* QR decomposition based on Gram=Schmidt's method. */
int zMatDecompQR(zMat m, zMat q, zMat r, zIndex idx)
{
  zMat mcp, qcp, rcp;
  int rank = -1;

  mcp = zMatAlloc( zMatColSizeNC(m), zMatRowSizeNC(m) );
  qcp = zMatAlloc( zMatColSizeNC(q), zMatRowSizeNC(q) );
  rcp = r ? zMatAlloc( zMatColSizeNC(r), zMatRowSizeNC(r) ) : NULL;
  if( !mcp || !qcp || ( r && !rcp ) ){
    ZALLOCERROR();
    goto TERMINATE;
  }
  zMatTNC( m, mcp );
  rank = zMatDecompLQDST( mcp, rcp, qcp, idx );
  zMatTNC( qcp, q );
  if( r ) zMatTNC( rcp, r );

 TERMINATE:
  zMatFreeAO( 3, mcp, qcp, rcp );
  return rank;
}
