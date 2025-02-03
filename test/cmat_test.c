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

int main(void)
{
  zRandInit();
  assert_get_put();
  return EXIT_SUCCESS;
}
