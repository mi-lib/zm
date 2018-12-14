/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_lq - linear equation: LQ/QR decomposition.
 */

#include <zm/zm_le.h>

/* zLQDecompDST
 * - LQ decomposition based on Gram=Schmidt's method.
 *   (destructive)
 */
int zLQDecompDST(zMat m, zMat l, zMat q, zIndex idx)
{
  register int i, j, rank;
  double *mp, r;

  zMatClear( l );
  zIndexOrder( idx, 0 );
  for( rank=0, i=0; i<zMatRowSizeNC(m); i++ ){
    mp = zMatRowBuf(m,i);
    for( j=0; j<rank; j++ ){
      r = zRawVecInnerProd( mp, zMatRowBuf(q,j), zMatColSizeNC(m) );
      zRawVecCatDRC( mp, -r, zMatRowBuf(q,j), zMatColSizeNC(q) );
      zMatElem(l,i,j) = r;
    }
    if( zIsTiny( r = zRawVecNorm(mp,zMatColSizeNC(m)) ) ){
      zIndexMove( idx, rank, zArrayNum(idx)-1 );
      continue;
    }
    zRawVecDiv( mp, r, zMatRowBuf(q,rank), zMatColSizeNC(m) );
    zMatElem(l,i,rank) = r;
    if( rank < zMatColSizeNC(m) ) rank++;
  }
  return rank;
}

/* zLQDecomp
 * - LQ decomposition based on Gram=Schmidt's method.
 */
int zLQDecomp(zMat m, zMat l, zMat q, zIndex idx)
{
  zMat mcp;
  int rank;

  if( !( mcp = zMatClone( m ) ) ){
    ZALLOCERROR();
    return 0;
  }
  rank = zLQDecompDST( mcp, l, q, idx );
  zMatFree( mcp );
  return rank;
}

/* zLQDecompReg
 * - LQ decomposition and regression.
 */
int zLQDecompReg(zMat m, zMat l, zMat q, zIndex idx)
{
  int rank;

  if( ( rank = zLQDecomp( m, l, q, idx ) ) < zMatRowSizeNC(m) ){
    zMatColReg( l, rank );
    zMatRowReg( q, rank );
  }
  return rank;
}

/* zLQDecompAlloc
 * - LQ decomposition with an automatic matrix allocation and resize.
 */
int zLQDecompAlloc(zMat m, zMat *l, zMat *q, zIndex *idx)
{
  *l = zMatAllocSqr( zMatRowSizeNC(m) );
  *q = zMatAlloc( zMatRowSizeNC(m), zMatColSizeNC(m) );
  *idx = zIndexCreate( zMatRowSizeNC(m) );
  if( !*l || !*q || !*idx ){
    zMatFree( *l );
    zMatFree( *q );
    zIndexFree( *idx );
    return 0;
  }
  return zLQDecompReg( m, *l, *q, *idx );
}

/* zQRDecomp
 * - QR decomposition based on Gram=Schmidt's method.
 */
int zQRDecomp(zMat m, zMat q, zMat r, zIndex idx)
{
  zMat mcp, qcp, rcp;
  int rank = 0;

  mcp = zMatAlloc( zMatColSizeNC(m), zMatRowSizeNC(m) );
  qcp = zMatAlloc( zMatColSizeNC(q), zMatRowSizeNC(q) );
  rcp = r ? zMatAlloc( zMatColSizeNC(r), zMatRowSizeNC(r) ) : NULL;
  if( !mcp || !qcp || ( r && !rcp ) ){
    ZALLOCERROR();
    goto TERMINATE;
  }
  zMatTNC( m, mcp );
  rank = zLQDecompDST( mcp, rcp, qcp, idx );
  zMatTNC( qcp, q );
  if( r ) zMatTNC( rcp, r );

 TERMINATE:
  zMatFreeAO( 3, mcp, qcp, rcp );
  return rank;
}
