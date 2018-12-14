/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_pivot - linear equation: pivoting.
 */

#include <zm/zm_le.h>

/* zPivoting
 * - matrix pivoting.
 */
int zPivoting(zMat m, zIndex idx, int r, int c)
{
  register int i;
  int pi;
  double tmp, max;

  pi = r;
  max = fabs( zMatElem( m, zIndexElem(idx,pi), c ) );
  for( i=r+1; i<zArrayNum(idx); i++ ){
    tmp = fabs( zMatElem( m, zIndexElem(idx,i), c ) );
    if( tmp > max ){
      max = tmp;
      pi = i;
    }
  }
  return zIndexSwap( idx, pi, r );
}

/* zPivotingDiag
 * - matrix pivoting in diagonal values.
 */
int zPivotingDiag(zMat m, zIndex idx, int i)
{
  register int j;
  int pi;
  double tmp, max;

  pi = i;
  max = fabs( zMatElem( m, zIndexElem(idx,i), zIndexElem(idx,i) ) );
  for( j=i+1; j<zArrayNum(idx); j++ ){
    tmp = fabs( zMatElem( m, zIndexElem(idx,j), zIndexElem(idx,j) ) );
    if( tmp > max ){
      max = tmp;
      pi = j;
    }
  }
  return zIndexSwap( idx, pi, i );
}

/* zSweepOutMat
 * - sweep out matrix column.
 */
double zSweepOutMat(zMat m1, zMat m2, int r, int c)
{
  register int i, j;
  double value, d, ratio;

  value = zMatElem( m1, r, c );
  if( value == 0 ){
    ZRUNWARN( ZM_WARN_LE_ZEROPIVOT );
    return 0;
  }
  d = 1.0 / value;
  /* measurement */
  for( j=0; j<zMatColSizeNC(m1); j++ )
    j != c ? ( zMatElem(m1,r,j) *= d ) : zMatSetElem( m1, r, c, 1 );
  for( j=0; j<zMatColSizeNC(m2); j++ )
    zMatElem(m2,r,j) *= d;
  /* sweep out */
  for( i=0; i<zMatRowSizeNC(m1); i++ ){
    if( i == r ) continue;
    if( !zIsTiny( ( ratio = zMatElem( m1, i, c ) ) ) ){
      for( j=0; j<zMatColSizeNC(m1); j++ )
        if( j != c )
          zMatElem(m1,i,j) -= zMatElem(m1,r,j) * ratio;
      for( j=0; j<zMatColSizeNC(m2); j++ )
        zMatElem(m2,i,j) -= zMatElem(m2,r,j) * ratio;
    }
    zMatSetElem( m1, i, c, 0 );
  }
  return value;
}

/* zSweepOutVec
 * - sweep out vector.
 */
double zSweepOutVec(zMat m, zVec v, int r, int c)
{
  register int i, j;
  double value, d, ratio;

  value = zMatElem( m, r, c );
  if( value == 0 ){
    ZRUNWARN( ZM_WARN_LE_ZEROPIVOT );
    return 0;
  }
  d = 1.0 / value;
  /* measurement */
  for( j=0; j<zMatColSizeNC(m); j++ )
    j != c ? ( zMatElem(m,r,j) *= d ) : zMatSetElem( m, r, c, 1 );
  zVecElem(v,r) *= d;
  /* sweep out */
  for( i=0; i<zMatRowSizeNC(m); i++ ){
    if( i == r ) continue;
    if( !zIsTiny( ( ratio = zMatElem( m, i, c ) ) ) ){
      for( j=0; j<zMatColSizeNC(m); j++ )
        if( j != c )
          zMatElem(m,i,j) -= zMatElem(m,r,j) * ratio;
      zVecElem(v,i) -= zVecElem(v,r) * ratio;
    }
    zMatSetElem( m, i, c, 0 );
  }
  return value;
}
