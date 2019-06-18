/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_lyapnov - linear equation: Lyapnov equation.
 */

#include <zm/zm_le.h>

/* zLyapnovSolve
 * - Lyapnov equation solver.
 */
zMat zLyapnovSolve(zMat a, zMat b, zMat ans)
{
  zVecStruct bvec, ansvec;
  zVec s;
  zMat aw, bt;
  zIndex idx;
  register int n, nn, i, j, k;

  if( !zMatIsSqr(a) || !zMatIsSqr(b) || !zMatIsSqr(ans) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return NULL;
  }
  if( !zMatSizeIsEqual(a,ans) || !zMatSizeIsEqual(b,ans) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  n = zMatRowSizeNC(ans);
  nn = n * n;
  bt = zMatAllocSqr( n );
  aw = zMatAllocSqr( nn );
  idx = zIndexCreate( nn );
  s = zVecAlloc( nn );
  if( !bt || !aw || !idx || !s ){
    ZALLOCERROR();
    ans = NULL;
    goto TERMINATE;
  }
  zMatTNC( b, bt );
  zVecSetSizeNC( &bvec, nn );
  zVecBufNC(&bvec) = zMatBuf(bt);
  zVecSetSizeNC( &ansvec, nn );
  zVecBufNC(&ansvec) = zMatBuf(ans);
  for( i=0; i<n; i++ )
    for( j=0; j<n; j++ )
      for( k=0; k<n; k++ ){
        zMatElemNC(aw,n*i+k,n*k+j) += zMatElemNC(a,j,i);
        zMatElemNC(aw,n*k+i,n*j+k) += zMatElemNC(a,j,i);
      }
  if( !zLESolveGaussDST( aw, &bvec, &ansvec, idx, s ) ){
    ZRUNERROR( ZM_ERR_LE_SINGULAR );
    ans = NULL;
  }
 TERMINATE:
  zMatFree( bt );
  zMatFree( aw );
  zIndexFree( idx );
  zVecFree( s );
  return ans;
}
