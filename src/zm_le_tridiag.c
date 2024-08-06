/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_tridiag - linear equation: tridiagonal equation.
 */

#include <zm/zm_le.h>

/* tridiagonal equation solver (destructive). */
zVec zTridiagSolveDST(zVec a, zVec b, zVec c, zVec d, zVec ans)
{
  int i;
  int n;

  n = zVecSize( a );
  for( i=1; i<n; i++ ){
    if( zVecElemNC(b,i-1) == 0 ){
      ZRUNERROR( ZM_ERR_LE_SINGULAR );
      return NULL;
    }
    zVecElemNC(b,i) -= zVecElemNC(c,i-1) / zVecElemNC(b,i-1) * zVecElemNC(a,i);
    zVecElemNC(d,i) -= zVecElemNC(d,i-1) / zVecElemNC(b,i-1) * zVecElemNC(a,i);
    zVecSetElemNC( a, i, 0 );
  }
  zVecSetElemNC( ans, n-1, zVecElemNC(d,n-1)/zVecElemNC(b,n-1) );
  for( i=n-2; i>=0; i-- )
    zVecSetElemNC( ans, i,
      (zVecElemNC(d,i)-zVecElemNC(c,i)*zVecElemNC(ans,i+1)) / zVecElemNC(b,i) );
  return ans;
}

/* tridiagonal equation solver. */
zVec zTridiagSolve(zVec a, zVec b, zVec c, zVec d, zVec ans)
{
  zVec acp, bcp, dcp;

  if( !zVecSizeIsEqual(a,b) || !zVecSizeIsEqual(b,c) ||
      !zVecSizeIsEqual(c,d) || !zVecSizeIsEqual(d,ans) ){
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );
    return NULL;
  }
  acp = zVecClone( a );
  bcp = zVecClone( b );
  dcp = zVecClone( d );
  if( acp && bcp && dcp )
    ans = zTridiagSolveDST( acp, bcp, c, dcp, ans );
  zVecFreeAtOnce( 3, acp, bcp, dcp );
  return ans;
}
