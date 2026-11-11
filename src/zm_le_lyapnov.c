/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_lyapnov - linear equation: Lyapnov equation.
 */

#include <zm/zm_le.h>

/* Lyapnov equation solver. */
zMat zLELyapnovSolve(const zMat a, const zMat b, zMat ans)
{
  zVecStruct bvec;
  zVec s, x;
  zMat aw, bt;
  zIndex idx;
  int n, nn, i, j, k;

  if( !zMatIsSqr(a) || !zMatIsSqr(b) || !zMatIsSqr(ans) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return NULL;
  }
  if( !zMatSizeEqual(a,ans) || !zMatSizeEqual(b,ans) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  n = zMatRowSizeNC(ans);
  nn = n * n;
  bt = zMatAllocSqr( n );
  aw = zMatAllocSqr( nn );
  idx = zIndexCreate( nn );
  s = zVecAlloc( nn );
  x = zVecAlloc( nn );
  if( !bt || !aw || !idx || !s || !x ){
    ans = NULL;
    goto TERMINATE;
  }
  zMatTNC( b, bt );
  zVecAssignArray( &bvec, nn, zMatBuf(bt) ); /* bt is newly allocated, so that its rows are contiguous. */
  for( i=0; i<n; i++ )
    for( j=0; j<n; j++ )
      for( k=0; k<n; k++ ){
        zMatElemNC(aw,n*i+k,n*k+j) += zMatElemNC(a,j,i);
        zMatElemNC(aw,n*k+i,n*j+k) += zMatElemNC(a,j,i);
      }
  if( !zLESolveGaussDST( aw, &bvec, x, idx, s ) ){
    ZRUNERROR( ZM_ERR_MAT_SINGULAR );
    ans = NULL;
  } else /* the rows of ans are aligned with its column capacity, which might be larger than n. */
    for( i=0; i<n; i++ )
      zRawVecCopy( zVecBufNC(x) + n*i, zMatRowBufNC(ans,i), n );
 TERMINATE:
  zMatFree( bt );
  zMatFree( aw );
  zIndexFree( idx );
  zVecFreeAtOnce( 2, s, x );
  return ans;
}
