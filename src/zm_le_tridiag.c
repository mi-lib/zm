/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_tridiag - linear equation: tridiagonal equation.
 */

#include <zm/zm_le.h>

/* zTridiagSolveDST
 * - tridiagonal equation solver (destructive).
 */
zVec zTridiagSolveDST(zVec a, zVec b, zVec c, zVec d, zVec ans)
{
  register int i;
  int n;

  n = zVecSize( a );
  for( i=1; i<n; i++ ){
    if( zVecElem(b,i-1) == 0 ){
      ZRUNERROR( ZM_ERR_LE_SINGULAR );
      return NULL;
    }
    zVecElem(b,i) -= zVecElem(c,i-1) / zVecElem(b,i-1) * zVecElem(a,i);
    zVecElem(d,i) -= zVecElem(d,i-1) / zVecElem(b,i-1) * zVecElem(a,i);
    zVecSetElem( a, i, 0 );
  }
  zVecSetElem( ans, n-1, zVecElem(d,n-1)/zVecElem(b,n-1) );
  for( i=n-2; i>=0; i-- )
    zVecSetElem( ans, i,
      (zVecElem(d,i)-zVecElem(c,i)*zVecElem(ans,i+1)) / zVecElem(b,i) );
  return ans;
}

/* zTridiagSolve
 * - tridiagonal equation solver.
 */
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
  zVecFreeAO( 3, acp, bcp, dcp );
  return ans;
}
