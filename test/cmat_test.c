#include <zm/zm.h>

void assert_get_put(void)
{
  const int rowsize = 12;
  const int colsize = 10;
  zCMat mat_test1, mat_test2;
  zCVec vec_test1, vec_test2;
  int i;
  bool result;

  mat_test1 = zCMatAlloc( rowsize, colsize );
  mat_test2 = zCMatAlloc( rowsize, colsize );
  vec_test1 = zCVecAlloc( colsize );
  vec_test2 = zCVecAlloc( rowsize );
  zCMatRandUniform( mat_test1, -10, 10, -5, 5 );
  zCMatRandUniform( mat_test2, -10, 10, -5, 5 );

  zCMatGetRow( mat_test1, 3, vec_test1 );
  for( result=true, i=0; i<colsize; i++ )
    if( !zComplexMatch( zCMatElemNC(mat_test1,3,i), zCVecElemNC(vec_test1,i) ) ) result = false;
  zAssert( zCMatGetRow, result );
  zCMatPutRow( mat_test2, 2, vec_test1 );
  for( result=true, i=0; i<colsize; i++ )
    if( !zComplexMatch( zCMatElemNC(mat_test2,2,i), zCVecElemNC(vec_test1,i) ) ) result = false;
  zAssert( zCMatPutRow, result );
  zCMatGetCol( mat_test1, 2, vec_test2 );
  for( result=true, i=0; i<rowsize; i++ )
    if( !zComplexMatch( zCMatElemNC(mat_test1,i,2), zCVecElemNC(vec_test2,i) ) ) result = false;
  zAssert( zCMatGetCol, result );
  zCMatPutCol( mat_test2, 3, vec_test2 );
  for( result=true, i=0; i<rowsize; i++ )
    if( !zComplexMatch( zCMatElemNC(mat_test2,i,3), zCVecElemNC(vec_test2,i) ) ) result = false;
  zAssert( zCMatPutCol, result );
  zCMatCopy( mat_test1, mat_test2 );
  zCMatFree( mat_test1 );
  zCMatFree( mat_test2 );
  zCVecFree( vec_test1 );
  zCVecFree( vec_test2 );
}

void assert_cmul_cdiv_drc(void)
{
  zCMat m1, m2, m3;
  zComplex z;

  m1 = zCMatAlloc( 10, 20 );
  m2 = zCMatAlloc( 10, 20 );
  zCMatRandUniform( m1, -10, 10, -10, 10 );
  zComplexCreate( &z, zRandF(1,10), zRandF(-10,10) );
  m3 = zCMatClone( m1 );
  zCMatCMul( m1, &z, m2 );
  zCMatCMulDRC( m3, &z );
  zCMatSubDRC( m3, m2 );
  zAssert( zCMatCMulDRC, zCMatIsTol( m3, zTOL ) );
  zCMatCopy( m1, m3 );
  zCMatCDiv( m1, &z, m2 );
  zCMatCDivDRC( m3, &z );
  zCMatSubDRC( m3, m2 );
  zAssert( zCMatCDivDRC, zCMatIsTol( m3, zTOL ) );
  zCMatFree( m1 );
  zCMatFree( m2 );
  zCMatFree( m3 );
}

/* division by a complex number with a very small or large absolute value */
void assert_cdiv_scale(void)
{
  zCMat m1, m2;
  zComplex z, r;
  double scale[] = { 1.0e-7, 1.0e-200, 1.0e200 };
  int i, j, k;
  bool result = true, result_zero;

  m1 = zCMatAlloc( 10, 20 );
  m2 = zCMatAlloc( 10, 20 );
  zCMatRandUniform( m1, -10, 10, -10, 10 );
  for( k=0; k<sizeof(scale)/sizeof(double); k++ ){
    zComplexCreate( &z, zRandF(-10,10)*scale[k], zRandF(-10,10)*scale[k] );
    if( !zCMatCDiv( m1, &z, m2 ) ){
      result = false;
      continue;
    }
    for( i=0; i<zCMatRowSizeNC(m1); i++ )
      for( j=0; j<zCMatColSizeNC(m1); j++ ){ /* m2 z = m1 */
        zComplexCMul( zCMatElemNC(m2,i,j), &z, &r );
        zComplexSubDRC( &r, zCMatElemNC(m1,i,j) );
        if( zComplexAbs( &r ) > 10 * DBL_EPSILON * zComplexAbs( zCMatElemNC(m1,i,j) ) ) result = false;
      }
  }
  zComplexCreate( &z, 0, 0 );
  result_zero = zCMatCDiv( m1, &z, m2 ) == NULL;
  zCMatFree( m1 );
  zCMatFree( m2 );
  zAssert( zCMatCDiv (very small or large divisor), result );
  zAssert( zCMatCDiv (division by zero), result_zero );
}

/* division by a tiny real number */
void assert_div_tiny(void)
{
  zCMat m1, m2;
  int i, j;
  bool result = true, result_zero;

  m1 = zCMatAlloc( 10, 20 );
  m2 = zCMatAlloc( 10, 20 );
  zCMatRandUniform( m1, -1, 1, -1, 1 );
  if( !zCMatDiv( m1, 1.0e-13, m2 ) ) result = false;
  for( i=0; i<zCMatRowSizeNC(m1); i++ )
    for( j=0; j<zCMatColSizeNC(m1); j++ )
      if( fabs( zCMatElemNC(m2,i,j)->re - zCMatElemNC(m1,i,j)->re/1.0e-13 ) > DBL_EPSILON * fabs( zCMatElemNC(m1,i,j)->re/1.0e-13 ) ||
          fabs( zCMatElemNC(m2,i,j)->im - zCMatElemNC(m1,i,j)->im/1.0e-13 ) > DBL_EPSILON * fabs( zCMatElemNC(m1,i,j)->im/1.0e-13 ) )
        result = false;
  result_zero = zCMatDiv( m1, 0, m2 ) == NULL;
  zCMatFree( m1 );
  zCMatFree( m2 );
  zAssert( zCMatDiv (tiny divisor), result );
  zAssert( zCMatDiv (division by zero), result_zero );
}

int main(void)
{
  zRandInit();
  assert_get_put();
  assert_cmul_cdiv_drc();
  assert_cdiv_scale();
  assert_div_tiny();
  return EXIT_SUCCESS;
}
