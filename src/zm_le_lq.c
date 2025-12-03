/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_lq - linear equation: LQ/QR decomposition.
 */

#include <zm/zm_le.h>

/* LQ decomposition based on Gram=Schmidt's method. (destructive) */
int zMatDecompLQ_GramSchmidt_DST(zMat m, zMat l, zMat q)
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
  return rank;
}

/* LQ decomposition based on Gram=Schmidt's method. */
int zMatDecompLQ_GramSchmidt(const zMat m, zMat l, zMat q)
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

/* LQ decomposition based on Householder method. (destructive) */
static int _zMatDecompLQ_Householder_DST(zMat m, zMat q)
{
  double s, ds, norm_inv, reflection;
  double *u;
  int i, j, size, colsize;

  zMatIdentNC( q );
  size = zMatMinSize( m );
  for( i=0; i<size; i++ ){
    u = &zMatElemNC(m,i,i);
    colsize = zMatColSizeNC(m) - i;
    if( zIsTiny( ( s = zRawVecNorm( u, colsize ) ) ) ) break;
    if( zSgn( zMatElemNC(m,i,i) ) < 0 ) s = -s;
    if( zIsTiny( ( ds = s - zMatElemNC(m,i,i) ) ) ) continue;
    norm_inv = 1.0 / ( s * ds );
    *u = -ds;
    for( j=0; j<zMatColSizeNC(q); j++ ){
      reflection = -norm_inv * zRawMatColInnerProd( zMatRowBufNC(q,i), u, zMatRowSizeNC(q)-i, zMatColSizeNC(q), j );
      zRawMatColCatDRC( zMatRowBufNC(q,i), reflection, u, zMatRowSizeNC(q)-i, zMatColSizeNC(q), j );
    }
    for( j=zMatRowSizeNC(m)-1; j>i; j-- ){
      reflection = -norm_inv * zRawVecInnerProd( &zMatElemNC(m,j,i), u, colsize );
      zRawVecCatDRC( &zMatElemNC(m,j,i), reflection, u, colsize );
    }
    *u = s;
    zRawVecZero( u + 1, colsize - 1 );
  }
  return i; /* rank */
}

/* LQ decomposition based on Householder method. */
int zMatDecompLQ_Householder(const zMat m, zMat l, zMat q)
{
  zMat ltmp, qtmp;
  int rank = -1;

  ltmp = zMatClone( m );
  qtmp = zMatAllocSqr( zMatColSizeNC(m) );
  if( !ltmp || !qtmp ){
    ZALLOCERROR();
    goto TERMINATE;
  }
  if( ( rank = _zMatDecompLQ_Householder_DST( ltmp, qtmp ) ) <= 0 ){
    ZRUNERROR( ZM_ERR_MAT_CANNOTDECOMPOSEZEROMAT );
    goto TERMINATE;
  }
  zMatColResize( l, rank );
  zMatRowResize( q, rank );
  if( l )
    zMatGet( ltmp, 0, 0, l );
  zMatGet( qtmp, 0, 0, q );

 TERMINATE:
  zMatFreeAtOnce( 2, ltmp, qtmp );
  return rank;
}

/* LQ decomposition and resizing of a matrix. */
int zMatDecompLQAndResize(const zMat m, zMat l, zMat q)
{
  int rank;

  if( ( rank = zMatDecompLQ( m, l, q ) ) < 0 ) return -1;
  if( rank < zMatRowSizeNC(m) ){
    zMatColResize( l, rank );
    zMatRowResize( q, rank );
  }
  return rank;
}

/* LQ decomposition with an automatic matrix allocation and resizing. */
int zMatDecompLQAlloc(const zMat m, zMat *l, zMat *q)
{
  *l = zMatAllocSqr( zMatRowSizeNC(m) );
  *q = zMatAlloc( zMatRowSizeNC(m), zMatColSizeNC(m) );
  if( !*l || !*q ){
    zMatFree( *l );
    zMatFree( *q );
    return -1;
  }
  return zMatDecompLQAndResize( m, *l, *q );
}

/* QR decomposition. */
int zMatDecompQR(const zMat m, zMat q, zMat r)
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
  rank = zMatDecompLQDST( mcp, rcp, qcp );
  zMatTNC( qcp, q );
  if( r ) zMatTNC( rcp, r );

 TERMINATE:
  zMatFreeAtOnce( 3, mcp, qcp, rcp );
  return rank;
}
