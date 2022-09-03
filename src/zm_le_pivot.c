/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_pivot - linear equation: pivoting.
 */

#include <zm/zm_le.h>

/* matrix pivoting. */
int zPivoting(zMat m, zIndex idx, int r, int c)
{
  int i, pi;
  double tmp, max;

  pi = r;
  max = fabs( zMatElemNC( m, zIndexElemNC(idx,pi), c ) );
  for( i=r+1; i<zArraySize(idx); i++ ){
    tmp = fabs( zMatElemNC( m, zIndexElemNC(idx,i), c ) );
    if( tmp > max ){
      max = tmp;
      pi = i;
    }
  }
  return zIndexSwap( idx, pi, r );
}

/* matrix pivoting in diagonal values. */
int zPivotingDiag(zMat m, zIndex idx, int i)
{
  int j, pi;
  double tmp, max;

  pi = i;
  max = fabs( zMatElemNC( m, zIndexElemNC(idx,i), zIndexElemNC(idx,i) ) );
  for( j=i+1; j<zArraySize(idx); j++ ){
    tmp = fabs( zMatElemNC( m, zIndexElemNC(idx,j), zIndexElemNC(idx,j) ) );
    if( tmp > max ){
      max = tmp;
      pi = j;
    }
  }
  return zIndexSwap( idx, pi, i );
}

/* sweep out matrix column. */
double zSweepOutMat(zMat m1, zMat m2, int r, int c)
{
  int i, j;
  double value, d, ratio;

  value = zMatElemNC( m1, r, c );
  if( value == 0 ){
    ZRUNWARN( ZM_WARN_LE_ZEROPIVOT );
    return 0;
  }
  d = 1.0 / value;
  /* measurement */
  for( j=0; j<zMatColSizeNC(m1); j++ )
    j != c ? ( zMatElemNC(m1,r,j) *= d ) : zMatSetElemNC( m1, r, c, 1 );
  for( j=0; j<zMatColSizeNC(m2); j++ )
    zMatElemNC(m2,r,j) *= d;
  /* sweep out */
  for( i=0; i<zMatRowSizeNC(m1); i++ ){
    if( i == r ) continue;
    if( !zIsTiny( ( ratio = zMatElemNC( m1, i, c ) ) ) ){
      for( j=0; j<zMatColSizeNC(m1); j++ )
        if( j != c )
          zMatElemNC(m1,i,j) -= zMatElemNC(m1,r,j) * ratio;
      for( j=0; j<zMatColSizeNC(m2); j++ )
        zMatElemNC(m2,i,j) -= zMatElemNC(m2,r,j) * ratio;
    }
    zMatSetElemNC( m1, i, c, 0 );
  }
  return value;
}

/* sweep out vector. */
double zSweepOutVec(zMat m, zVec v, int r, int c)
{
  int i, j;
  double value, d, ratio;

  value = zMatElemNC( m, r, c );
  if( value == 0 ){
    ZRUNWARN( ZM_WARN_LE_ZEROPIVOT );
    return 0;
  }
  d = 1.0 / value;
  /* measurement */
  for( j=0; j<zMatColSizeNC(m); j++ )
    j != c ? ( zMatElemNC(m,r,j) *= d ) : zMatSetElemNC( m, r, c, 1 );
  zVecElemNC(v,r) *= d;
  /* sweep out */
  for( i=0; i<zMatRowSizeNC(m); i++ ){
    if( i == r ) continue;
    if( !zIsTiny( ( ratio = zMatElemNC( m, i, c ) ) ) ){
      for( j=0; j<zMatColSizeNC(m); j++ )
        if( j != c )
          zMatElemNC(m,i,j) -= zMatElemNC(m,r,j) * ratio;
      zVecElemNC(v,i) -= zVecElemNC(v,r) * ratio;
    }
    zMatSetElemNC( m, i, c, 0 );
  }
  return value;
}
