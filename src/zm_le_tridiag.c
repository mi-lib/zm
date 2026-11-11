/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_tridiag - linear equation: tridiagonal equation.
 */

#include <zm/zm_le.h>

/* tridiagonal equation solver based on Gaussian elimination with partial pivoting (destructive).
 * At the i-th step, the i-th row (with the components at the columns i and i+1) and the (i+1)-th row
 * (with the components at the columns i, i+1 and i+2) are compared, and the one with the larger
 * absolute component at the column i is chosen as the pivot row. The resulting upper triangular matrix
 * has two superdiagonals due to the row exchanges, whose components are stored in a (the first) and
 * ans (the second), while the diagonal components are stored in b. c is not modified. */
zVec zLETridiagSolveDST(zVec a, zVec b, zVec c, zVec d, zVec ans)
{
  int i, n;
  double p0, p1, p2, pd, q0, q1, q2, qd, c_cur, x;

  if( ( n = zVecSizeNC(b) ) == 0 ) return ans;
  c_cur = n > 1 ? zVecElemNC(c,0) : 0; /* superdiagonal component of the current row */
  for( i=0; i<n-1; i++ ){
    if( fabs( zVecElemNC(a,i+1) ) > fabs( zVecElemNC(b,i) ) ){ /* exchange the rows */
      p0 = zVecElemNC(a,i+1); p1 = zVecElemNC(b,i+1); p2 = i < n-2 ? zVecElemNC(c,i+1) : 0; pd = zVecElemNC(d,i+1);
      q0 = zVecElemNC(b,i);   q1 = c_cur;             q2 = 0;                               qd = zVecElemNC(d,i);
    } else{
      p0 = zVecElemNC(b,i);   p1 = c_cur;             p2 = 0;                               pd = zVecElemNC(d,i);
      q0 = zVecElemNC(a,i+1); q1 = zVecElemNC(b,i+1); q2 = i < n-2 ? zVecElemNC(c,i+1) : 0; qd = zVecElemNC(d,i+1);
    }
    if( p0 == 0 ){
      ZRUNERROR( ZM_ERR_MAT_SINGULAR );
      return NULL;
    }
    /* the i-th row of the upper triangular matrix (a(i) is not used any more) */
    zVecSetElemNC( b, i, p0 );
    zVecSetElemNC( a, i, p1 );
    zVecSetElemNC( ans, i, p2 );
    zVecSetElemNC( d, i, pd );
    /* elimination of the other row, which becomes the (i+1)-th row */
    zVecSetElemNC( b, i+1, q1 - p1 / p0 * q0 );
    c_cur = q2 - p2 / p0 * q0;
    zVecSetElemNC( d, i+1, qd - pd / p0 * q0 );
  }
  if( zVecElemNC(b,n-1) == 0 ){
    ZRUNERROR( ZM_ERR_MAT_SINGULAR );
    return NULL;
  }
  /* backward substitution, where the i-th component of ans stores the second superdiagonal component
   * of the i-th row until it is overwritten by the answer. */
  zVecSetElemNC( ans, n-1, zVecElemNC(d,n-1) / zVecElemNC(b,n-1) );
  for( i=n-2; i>=0; i-- ){
    x = zVecElemNC(d,i) - zVecElemNC(a,i) * zVecElemNC(ans,i+1);
    if( i < n-2 && zVecElemNC(ans,i) != 0 )
      x -= zVecElemNC(ans,i) * zVecElemNC(ans,i+2);
    zVecSetElemNC( ans, i, x / zVecElemNC(b,i) );
  }
  return ans;
}

/* tridiagonal equation solver. */
zVec zLETridiagSolve(const zVec a, const zVec b, const zVec c, const zVec d, zVec ans)
{
  zVec acp, bcp, dcp;

  if( !zVecSizeEqual(a,b) || !zVecSizeEqual(b,c) ||
      !zVecSizeEqual(c,d) || !zVecSizeEqual(d,ans) ){
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );
    return NULL;
  }
  acp = zVecClone( a );
  bcp = zVecClone( b );
  dcp = zVecClone( d );
  if( acp && bcp && dcp )
    ans = zLETridiagSolveDST( acp, bcp, c, dcp, ans );
  zVecFreeAtOnce( 3, acp, bcp, dcp );
  return ans;
}
