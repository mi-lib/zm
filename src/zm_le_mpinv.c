/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_mpinv - linear equation: Moore=Penrose inverse matrix.
 */

#include <zm/zm_le.h>

static int _zMPInvAllocWork1(zMat m, zMat *l, zMat *q, zIndex *idx);
static int _zMPInvAllocWork2(int rank, int row, zMat *tmp1, zMat *tmp2, zMat *tmp3);

/* (static)
 * _zMPInvAllocWork1
 * - initialization: LQ decomposition
 */
int _zMPInvAllocWork1(zMat m, zMat *l, zMat *q, zIndex *idx)
{
  int rank;

  if( ( rank = zLQDecompAlloc( m, l, q, idx ) ) == -1 ){
    zMatFreeAO( 2, *l, *q );
    zIndexFree( *idx );
    return -1;
  }
  return rank;
}

/* (static)
 * _zMPInvAllocWork2
 * - initialization: workspace for matrix computation
 */
int _zMPInvAllocWork2(int rank, int row, zMat *tmp1, zMat *tmp2, zMat *tmp3)
{
  *tmp1 = zMatAlloc( rank, row );
  *tmp2 = zMatAllocSqr( rank );
  *tmp3 = zMatAlloc( rank, row );
  if( !*tmp1 || !*tmp2 || !*tmp3 ){
    zMatFreeAO( 3, *tmp1, *tmp2, *tmp3 );
    return -1;
  }
  return rank;
}

/* zMPInv
 * - Moore=Penrose's inverse matrix.
 */
int zMPInv(zMat m, zMat mp)
{
  int rank;
  zMat l, q, tmp1, tmp2, tmp3;
  zIndex idx;

  if( ( rank = _zMPInvAllocWork1( m, &l, &q, &idx ) ) == -1 )
    return -1;
  if( _zMPInvAllocWork2( rank, zMatRowSizeNC(l), &tmp1, &tmp2, &tmp3 ) == -1 )
    return -1;

  if( zMatIsSqr(l) )
    zMatInv( l, tmp3 );
  else{
    zMatT( l, tmp1 );
    zMulMatTMat( l, l, tmp2 );
    zMulInvMatMat( tmp2, tmp1, tmp3 );
  }
  zMulMatTMat( q, tmp3, mp );

  zMatFreeAO( 5, l, q, tmp1, tmp2, tmp3 );
  zIndexFree( idx );
  return rank;
}

/* zMPInvNull
 * - Moore=Penrose's inverse matrix with its null space.
 */
int zMPInvNull(zMat m, zMat mp, zMat mn)
{
  int rank;
  zMat l, q, tmp1, tmp2, tmp3;
  zIndex idx;
  register int i;

  if( ( rank = _zMPInvAllocWork1( m, &l, &q, &idx ) ) == -1 )
    return -1;
  if( _zMPInvAllocWork2( rank, zMatRowSizeNC(l), &tmp1, &tmp2, &tmp3 ) == -1 )
    return -1;

  if( zMatIsSqr(l) )
    zMatInv( l, tmp3 );
  else{
    zMatT( l, tmp1 );
    zMulMatTMat( l, l, tmp2 );
    zMulInvMatMat( tmp2, tmp1, tmp3 );
  }
  zMulMatTMat( q, tmp3, mp );
  /* null space */
  zMulMatTMat( q, q, mn );
  for( i=0; i<zMatRowSizeNC(mn); i++ ) zMatElem(mn,i,i) -= 1.0;

  zMatFreeAO( 5, l, q, tmp1, tmp2, tmp3 );
  zIndexFree( idx );
  return rank;
}

/* zMulMPInvMatMat
 * - multiply Moore=Penrose's inverse matrix from the left side to
 *   another matrix.
 */
zMat zMulMPInvMatMat(zMat m1, zMat m2, zMat m)
{
  int rank;
  zMat l, q, tmp1, tmp2, tmp3;
  zIndex idx;

  if( ( rank = _zMPInvAllocWork1( m1, &l, &q, &idx ) ) == -1 )
    return NULL;
  if( _zMPInvAllocWork2( rank, zMatColSizeNC(m2), &tmp1, &tmp2, &tmp3 ) == -1 )
    return NULL;

  if( zMatIsSqr(l) )
    zMulInvMatMat( l, m2, tmp3 );
  else{
    zMulMatTMat( l, m2, tmp1 );
    zMulMatTMat( l, l, tmp2 );
    zMulInvMatMat( tmp2, tmp1, tmp3 );
  }
  zMulMatTMat( q, tmp3, m );

  zMatFreeAO( 5, l, q, tmp1, tmp2, tmp3 );
  zIndexFree( idx );
  return m;
}

/* Penrose's iterative algorithm */

static int _zMPInvPenrose(zMat m, zMat mp);

/* _zMPInvPenrose
 * - internal call of Moore=Penrose inverse matrix based on Penrose's
 *   iterative algorithm.
 */
int _zMPInvPenrose(zMat m, zMat mp)
{
  zMat b, c, c2, cb;
  double trace = 0;
  int rank = 0;

  b = zMatAllocSqr( zMatColSizeNC(m) );
  c = zMatAllocSqr( zMatColSizeNC(m) );
  c2 = zMatAllocSqr( zMatColSizeNC(m) );
  cb = zMatAllocSqr( zMatColSizeNC(m) );
  if( !b || !c || !c2 || !cb ){
    ZALLOCERROR();
    goto TERMINATE;
  }
  zMulMatTMatNC( m, m, b );
  zMatIdentNC( c2 );
  zMatCopyNC( b, cb );
  for( rank=1; !zMatIsTiny(cb); rank++ ){
    zMatCopyNC( c2, c );
    trace = zMatTrNC( cb );
    zMatIdentNC( c2 );
    zMatMulNCDRC( c2, trace / rank );
    zMatSubNCDRC( c2, cb );
    zMulMatMatNC( c2, b, cb );
  }
  zMatMulNCDRC( c, --rank / trace );
  zMulMatMatTNC( c, m, mp );
 TERMINATE:
  zMatFreeAO( 4, b, c, c2, cb );
  return rank;
}

/* zMPInvPenrose
 * - Moore=Penrose inverse matrix based on Penrose's iterative algorithm.
 */
int zMPInvPenrose(zMat m, zMat mp)
{
  int rank;

  if( zMatColSizeNC(m) > zMatRowSizeNC(m) ){
    zMatTDST( m );
    zMatTDST( mp );
    rank = _zMPInvPenrose( m, mp );
    zMatTDST( m );
    zMatTDST( mp );
  } else
    rank = _zMPInvPenrose( m, mp );
  return rank;
}
