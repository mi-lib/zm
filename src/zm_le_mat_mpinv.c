/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_mpinv - linear equation: Moore-Penrose inverse matrix.
 */

#include <zm/zm_le.h>

/* initialization: LQ decomposition */
static int _zMatMPInvAllocWork1(const zMat m, zMat *l, zMat *q, zIndex *idx)
{
  int rank;

  if( ( rank = zMatDecompLQAlloc( m, l, q, idx ) ) < 0 ){
    zMatFreeAtOnce( 2, *l, *q );
    zIndexFree( *idx );
    return -1;
  }
  return rank;
}

/* initialization: workspace for matrix computation */
static int _zMatMPInvAllocWork2(int rank, int row, zMat *tmp1, zMat *tmp2, zMat *tmp3)
{
  *tmp1 = zMatAlloc( rank, row );
  *tmp2 = zMatAllocSqr( rank );
  *tmp3 = zMatAlloc( rank, row );
  if( !*tmp1 || !*tmp2 || !*tmp3 ){
    zMatFreeAtOnce( 3, *tmp1, *tmp2, *tmp3 );
    return -1;
  }
  return rank;
}

/* Moore-Penrose's inverse matrix. */
int zMatMPInv(const zMat m, zMat mp)
{
  int rank;
  zMat l, q, tmp1, tmp2, tmp3;
  zIndex idx;

  if( ( rank = _zMatMPInvAllocWork1( m, &l, &q, &idx ) ) < 0 ) return -1;
  if( _zMatMPInvAllocWork2( rank, zMatRowSizeNC(l), &tmp1, &tmp2, &tmp3 ) < 0 ) return -1;
  if( zMatIsSqr( l ) )
    zMatInv( l, tmp3 );
  else{
    zMatT( l, tmp1 );
    zMulMatTMat( l, l, tmp2 );
    zMulInvMatMat( tmp2, tmp1, tmp3 );
  }
  zMulMatTMat( q, tmp3, mp );

  zMatFreeAtOnce( 5, l, q, tmp1, tmp2, tmp3 );
  zIndexFree( idx );
  return rank;
}

/* Moore-Penrose's inverse matrix with its null space. */
int zMatMPInvNull(const zMat m, zMat mp, zMat mn)
{
  int i, rank;
  zMat l, q, tmp1, tmp2, tmp3;
  zIndex idx;

  if( ( rank = _zMatMPInvAllocWork1( m, &l, &q, &idx ) ) < 0 ) return -1;
  if( _zMatMPInvAllocWork2( rank, zMatRowSizeNC(l), &tmp1, &tmp2, &tmp3 ) < 0 ) return -1;
  if( zMatIsSqr( l ) ){
    zMatInv( l, tmp3 );
  } else{
    zMatT( l, tmp1 );
    zMulMatTMat( l, l, tmp2 );
    zMulInvMatMat( tmp2, tmp1, tmp3 );
  }
  zMulMatTMat( q, tmp3, mp );
  /* null space */
  zMulMatTMat( q, q, mn );
  for( i=0; i<zMatRowSizeNC(mn); i++ ) zMatElemNC(mn,i,i) -= 1.0;

  zMatFreeAtOnce( 5, l, q, tmp1, tmp2, tmp3 );
  zIndexFree( idx );
  return rank;
}

/* multiply Moore-Penrose's inverse matrix from the left side to another matrix. */
zMat zMulMPInvMatMat(const zMat m1, const zMat m2, zMat m)
{
  int rank;
  zMat l, q, tmp1, tmp2, tmp3;
  zIndex idx;

  if( ( rank = _zMatMPInvAllocWork1( m1, &l, &q, &idx ) ) < 0 ) return NULL;
  if( _zMatMPInvAllocWork2( rank, zMatColSizeNC(m2), &tmp1, &tmp2, &tmp3 ) < 0 ) return NULL;
  if( zMatIsSqr( l ) )
    zMulInvMatMat( l, m2, tmp3 );
  else{
    zMulMatTMat( l, m2, tmp1 );
    zMulMatTMat( l, l, tmp2 );
    zMulInvMatMat( tmp2, tmp1, tmp3 );
  }
  zMulMatTMat( q, tmp3, m );

  zMatFreeAtOnce( 5, l, q, tmp1, tmp2, tmp3 );
  zIndexFree( idx );
  return m;
}

/* Penrose's iterative algorithm */

/* internal call of Moore-Penrose inverse matrix based on Penrose's iterative algorithm. */
static int _zMatMPInvPenrose(const zMat m, zMat mp)
{
  zMat b, c, c2, cb;
  double trace = 0;
  int rank = 0;

  b = zMatAllocSqr( zMatColSizeNC(m) );
  c = zMatAllocSqr( zMatColSizeNC(m) );
  c2 = zMatAllocSqr( zMatColSizeNC(m) );
  cb = zMatAllocSqr( zMatColSizeNC(m) );
  if( !b || !c || !c2 || !cb ) goto TERMINATE;

  zMulMatTMatNC( m, m, b );
  zMatIdentNC( c2 );
  zMatCopyNC( b, cb );
  for( rank=1; !zMatIsTiny(cb); rank++ ){
    zMatCopyNC( c2, c );
    trace = zMatTraceNC( cb );
    zMatIdentNC( c2 );
    zMatMulNCDRC( c2, trace / rank );
    zMatSubNCDRC( c2, cb );
    zMulMatMatNC( c2, b, cb );
  }
  zMatMulNCDRC( c, --rank / trace );
  zMulMatMatTNC( c, m, mp );
 TERMINATE:
  zMatFreeAtOnce( 4, b, c, c2, cb );
  return rank;
}

/* Moore-Penrose inverse matrix based on Penrose's iterative algorithm. */
int zMatMPInvPenrose(const zMat m, zMat mp)
{
  int rank;

  if( zMatColSizeNC(m) > zMatRowSizeNC(m) ){
    zMatTDRC( m );
    zMatTDRC( mp );
    rank = _zMatMPInvPenrose( m, mp );
    zMatTDRC( m );
    zMatTDRC( mp );
  } else
    rank = _zMatMPInvPenrose( m, mp );
  return rank;
}
