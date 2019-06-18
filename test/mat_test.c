#include <zm/zm.h>

#define MAT_ROW_SIZE 12
#define MAT_COL_SIZE 10

void assert_get_put(void)
{
  const int rowsize = MAT_ROW_SIZE;
  const int colsize = MAT_COL_SIZE;
  zMat mat_test1, mat_test2, mat_test3;
  zVec vec_test1, vec_test2;
  register int i, j;
  bool result;

  mat_test1 = zMatAlloc( rowsize, colsize );
  mat_test2 = zMatAlloc( rowsize, colsize );
  mat_test3 = zMatAlloc( 3, 2 );
  vec_test1 = zVecAlloc( colsize );
  vec_test2 = zVecAlloc( rowsize );
  zMatRandUniform( mat_test1, -10, 10 );
  zMatRandUniform( mat_test2, -10, 10 );
  zMatGet( mat_test1, 4, 3, mat_test3 );
  for( result=true, i=0; i<3; i++ )
    for( j=0; j<2; j++ )
      if( zMatElemNC(mat_test1,i+4,j+3) != zMatElemNC(mat_test3,i,j) ) result = false;
  zAssert( zMatGet, result );
  zMatPut( mat_test2, 5, 6, mat_test3 );
  for( result=true, i=0; i<3; i++ )
    for( j=0; j<2; j++ )
      if( zMatElemNC(mat_test2,i+5,j+6) != zMatElemNC(mat_test3,i,j) ) result = false;
  zAssert( zMatPut, result );
  zMatTGet( mat_test1, 4, 3, mat_test3 );
  for( result=true, i=0; i<3; i++ )
    for( j=0; j<2; j++ )
      if( zMatElemNC(mat_test1,4+j,i+3) != zMatElemNC(mat_test3,i,j) ) result = false;
  zAssert( zMatTGet, result );
  zMatTPut( mat_test2, 5, 6, mat_test3 );
  for( result=true, i=0; i<3; i++ )
    for( j=0; j<2; j++ )
      if( zMatElemNC(mat_test2,j+5,i+6) != zMatElemNC(mat_test3,i,j) ) result = false;
  zAssert( zMatTPut, result );

  zMatGetRow( mat_test1, 3, vec_test1 );
  zAssert( zMatGetRow, memcmp( zMatRowBufNC(mat_test1,3), zVecBufNC(vec_test1), sizeof(double)*colsize ) == 0 );
  zMatPutRow( mat_test2, 2, vec_test1 );
  zAssert( zMatPutRow, memcmp( zMatRowBufNC(mat_test2,2), zMatRowBufNC(mat_test1,3), sizeof(double)*colsize ) == 0 );
  zMatGetCol( mat_test1, 2, vec_test2 );
  for( result=true, i=0; i<rowsize; i++ )
    if( zMatElemNC(mat_test1,i,2) != zVecElemNC(vec_test2,i) ) result = false;
  zAssert( zMatGetCol, result );
  zMatPutCol( mat_test2, 3, vec_test2 );
  for( result=true, i=0; i<rowsize; i++ )
    if( zMatElemNC(mat_test2,i,3) != zMatElemNC(mat_test1,i,2) ) result = false;
  zAssert( zMatPutCol, result );
  zMatCopy( mat_test1, mat_test2 );
  zMatSwapRow( mat_test2, 3, 7 );
  zAssert( zMatSwapRow,
    memcmp( zMatRowBufNC(mat_test1,3), zMatRowBufNC(mat_test2,7), sizeof(double)*colsize ) == 0 &&
    memcmp( zMatRowBufNC(mat_test1,7), zMatRowBufNC(mat_test2,3), sizeof(double)*colsize ) == 0 );
  zMatCopy( mat_test1, mat_test2 );
  zMatSwapCol( mat_test2, 3, 7 );
  for( result=true, i=0; i<rowsize; i++ )
    if( zMatElemNC(mat_test1,i,3) != zMatElemNC(mat_test2,i,7) ||
        zMatElemNC(mat_test1,i,7) != zMatElemNC(mat_test2,i,3) ) result = false;
  zAssert( zMatSwapCol, result );
  zMatFreeAO( 3, mat_test1, mat_test2, mat_test3 );
  zVecFreeAO( 2, vec_test1, vec_test2 );
}

void assert_arith(void)
{
  const int rowsize = MAT_ROW_SIZE;
  const int colsize = MAT_COL_SIZE;
  zMat mat_test1, mat_test2, mat_test3;
  double k;
  register int i, j;
  bool result;

  mat_test1 = zMatAlloc( rowsize, colsize );
  mat_test2 = zMatAlloc( rowsize, colsize );
  mat_test3 = zMatAlloc( rowsize, colsize );
  zMatRandUniform( mat_test1, -10, 10 );
  zMatRandUniform( mat_test2, -10, 10 );
  k = zRandF( -10, 10 );
  zMatAdd( mat_test1, mat_test2, mat_test3 );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( !zIsTiny( zMatElemNC(mat_test1,i,j)+zMatElemNC(mat_test2,i,j)-zMatElemNC(mat_test3,i,j) ) ) result = false;
  zAssert( zMatAdd, result );
  zMatSub( mat_test1, mat_test2, mat_test3 );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( !zIsTiny( zMatElemNC(mat_test1,i,j)-zMatElemNC(mat_test2,i,j)-zMatElemNC(mat_test3,i,j) ) ) result = false;
  zAssert( zMatSub, result );
  zMatRev( mat_test1, mat_test3 );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( !zIsTiny( zMatElemNC(mat_test1,i,j)+zMatElemNC(mat_test3,i,j) ) ) result = false;
  zAssert( zMatRev, result );
  zMatMul( mat_test1, k, mat_test3 );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( !zIsTiny( zMatElemNC(mat_test1,i,j)*k-zMatElemNC(mat_test3,i,j) ) ) result = false;
  zAssert( zMatMul, result );
  zMatDiv( mat_test1, k, mat_test3 );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( !zIsTiny( zMatElemNC(mat_test1,i,j)/k-zMatElemNC(mat_test3,i,j) ) ) result = false;
  zAssert( zMatDiv, result );
  zMatCat( mat_test1, k, mat_test2, mat_test3 );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( !zIsTiny( zMatElemNC(mat_test1,i,j)+k*zMatElemNC(mat_test2,i,j)-zMatElemNC(mat_test3,i,j) ) ) result = false;
  zAssert( zMatCat, result );

  zMatCopy( mat_test1, mat_test3 );
  zMatAddDRC( mat_test3, mat_test2 );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( !zIsTiny( zMatElemNC(mat_test1,i,j)+zMatElemNC(mat_test2,i,j)-zMatElemNC(mat_test3,i,j) ) ) result = false;
  zAssert( zMatAddDRC, result );
  zMatCopy( mat_test1, mat_test3 );
  zMatSubDRC( mat_test3, mat_test2 );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( !zIsTiny( zMatElemNC(mat_test1,i,j)-zMatElemNC(mat_test2,i,j)-zMatElemNC(mat_test3,i,j) ) ) result = false;
  zAssert( zMatSubDRC, result );
  zMatCopy( mat_test1, mat_test3 );
  zMatRevDRC( mat_test3 );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( !zIsTiny( zMatElemNC(mat_test1,i,j)+zMatElemNC(mat_test3,i,j) ) ) result = false;
  zAssert( zMatRevDRC, result );
  zMatCopy( mat_test1, mat_test3 );
  zMatMulDRC( mat_test3, k );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( !zIsTiny( zMatElemNC(mat_test1,i,j)*k-zMatElemNC(mat_test3,i,j) ) ) result = false;
  zAssert( zMatMulDRC, result );
  zMatCopy( mat_test1, mat_test3 );
  zMatDivDRC( mat_test3, k );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( !zIsTiny( zMatElemNC(mat_test1,i,j)/k-zMatElemNC(mat_test3,i,j) ) ) result = false;
  zAssert( zMatDivDRC, result );
  zMatCopy( mat_test1, mat_test3 );
  zMatCatDRC( mat_test3, k, mat_test2 );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( !zIsTiny( zMatElemNC(mat_test1,i,j)+k*zMatElemNC(mat_test2,i,j)-zMatElemNC(mat_test3,i,j) ) ) result = false;
  zAssert( zMatCatDRC, result );
  zMatFreeAO( 3, mat_test1, mat_test2, mat_test3 );
}

void assert_transpose(void)
{
  const int rowsize = MAT_ROW_SIZE;
  const int colsize = MAT_COL_SIZE;
  zMat mat_test1, mat_test2;
  double tr;
  register int i, j;
  bool result;

  mat_test1 = zMatAlloc( colsize, rowsize );
  mat_test2 = zMatAlloc( rowsize, colsize );
  zMatRandUniform( mat_test1, -10, 10 );
  zMatT( mat_test1, mat_test2 );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( zMatElemNC(mat_test1,j,i) != zMatElemNC(mat_test2,i,j) ) result = false;
  zAssert( zMatT, result );
  zMatSetSize( mat_test2, colsize, rowsize );
  zMatCopy( mat_test1, mat_test2 );
  zMatTDRC( mat_test2 );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( zMatElemNC(mat_test1,j,i) != zMatElemNC(mat_test2,i,j) ) result = false;
  zAssert( zMatTDRC, result );
  for( tr=0, i=0; i<zMin(rowsize,colsize); i++ )
    tr += zMatElem(mat_test2,i,i);
  zAssert( zMatTrNC, zIsTiny( tr - zMatTrNC( mat_test2 ) ) );
  zMatFreeAO( 2, mat_test1, mat_test2 );
}

void assert_mul_mat_vec(void)
{
  zMat mat_test1, mat_test2, mat_test3, ans3, ans4, ans5, mat_error;
  zVec vec_test1, vec_test2, vec_test3, ans1, ans2, vec_error;

  mat_test1 = zMatCreateList( 3, 4,
    2.0, 1.0,-1.0, 2.0,
    1.0, 3.0, 2.0,-1.0,
   -1.0, 2.0,-3.0, 1.0 );
  mat_test2 = zMatCreateList( 4, 3,
    1.0, 1.0, 2.0,
    2.0,-2.0,-1.0,
   -1.0, 3.0, 2.0,
   -2.0, 1.0,-3.0 );
  vec_test1 = zVecCreateList( 4,
    1.0, 3.0, 4.0,-2.0 );
  vec_test2 = zVecCreateList( 3,
    2.0,-1.0, 3.0 );
  ans1 = zVecCreateList( 3,
   -3.0, 20.0, -9.0 );
  ans2 = zVecCreateList( 4, 0.0, 5.0,-13.0, 8.0 );
  ans3 = zMatCreateList( 3, 3,
    1.0, -1.0, -5.0,
    7.0,  0.0,  6.0,
    4.0,-13.0,-13.0 );
  ans4 = zMatCreateList( 3, 3,
    5.0,  2.0, -5.0,
    6.0,-10.0,  1.0,
   -3.0,  6.0,-12.0 );
  ans5 = zMatCreateList( 4, 4,
    -2.0, 3.0, 2.0, 10.0,
    -1.0,-6.0, 1.0,  5.0,
   -11.0, 3.0,-7.0, 13.0,
     6.0, 1.0, 6.0, -2.0 );

  vec_test3 = zVecAlloc( 3 );
  vec_error = zVecAlloc( 3 );
  zMulMatVec( mat_test1, vec_test1, vec_test3 );
  zVecSub( vec_test3, ans1, vec_error );
  zAssert( zMulMatVec, zVecIsTiny( vec_error ) );
  zVecFreeAO( 2, vec_test3, vec_error );
  vec_test3 = zVecAlloc( 4 );
  vec_error = zVecAlloc( 4 );
  zMulMatTVec( mat_test1, vec_test2, vec_test3 );
  zVecSub( vec_test3, ans2, vec_error );
  zAssert( zMulMatTVec, zVecIsTiny( vec_error ) );
  mat_test3 = zMatAllocSqr( 3 );
  mat_error = zMatAllocSqr( 3 );
  zMulMatMat( mat_test1, mat_test2, mat_test3 );
  zMatSub( mat_test3, ans3, mat_error );
  zAssert( zMulMatMat, zMatIsTiny( mat_error ) );
  zMatSetSize( mat_test2, 3, 4 );
  zMulMatMatT( mat_test1, mat_test2, mat_test3 );
  zMatSub( mat_test3, ans4, mat_error );
  zAssert( zMulMatMatT, zMatIsTiny( mat_error ) );
  zMatFreeAO( 2, mat_test3, mat_error );
  mat_test3 = zMatAllocSqr( 4 );
  mat_error = zMatAllocSqr( 4 );
  zMulMatTMat( mat_test1, mat_test2, mat_test3 );
  zMatSub( mat_test3, ans5, mat_error );
  zAssert( zMulMatTMat, zMatIsTiny( mat_error ) );
  zMatFreeAO( 7, mat_test1, mat_test2, mat_test3, ans3, ans4, ans5, mat_error );
  zVecFreeAO( 6, vec_test1, vec_test2, vec_test3, ans1, ans2, vec_error );
}

void assert_dyad(void)
{
  zVec vec_test1, vec_test2;
  zMat mat_test1, mat_test2, ans1, ans2, ans3, ans4, error;

  vec_test1 = zVecCreateList( 4, 1.0, 2.0, 3.0, 4.0 );
  vec_test2 = zVecCreateList( 3, 3.0, 2.0, 1.0 );
  mat_test1 = zMatCreateList( 4, 3,
    3.0, 1.0,-5.0,
    2.0,-4.0, 2.0,
    1.0,-3.0, 2.0,
   -2.0, 4.0,-1.0 );
  mat_test2 = zMatAlloc( 4, 3 );
  ans1 = zMatCreateList( 4, 3,
    3.0, 2.0, 1.0,
    6.0, 4.0, 2.0,
    9.0, 6.0, 3.0,
   12.0, 8.0, 4.0 );
  ans2 = zMatCreateList( 4, 3,
    6.0, 3.0,-4.0,
    8.0, 0.0, 4.0,
   10.0, 3.0, 5.0,
   10.0,12.0, 3.0 );
  ans3 = zMatCreateList( 4, 3,
    0.0,-1.0,-6.0,
   -4.0,-8.0, 0.0,
   -8.0,-9.0,-1.0,
  -14.0,-4.0,-5.0 );
  ans4 = zMatCreateList( 4, 3,
   -3.0, -3.0,-7.0,
  -10.0,-12.0,-2.0,
  -17.0,-15.0,-4.0,
  -26.0,-12.0,-9.0 );
  error = zMatAlloc( 4, 3 );

  zVecDyad( vec_test1, vec_test2, mat_test2 );
  zMatSub( mat_test2, ans1, error );
  zAssert( zVecDyad, zMatIsTiny( error ) );
  zMatCopy( mat_test1, mat_test2 );
  zMatAddDyad( mat_test2, vec_test1, vec_test2 );
  zMatSub( mat_test2, ans2, error );
  zAssert( zMatAddDyad, zMatIsTiny( error ) );
  zMatCopy( mat_test1, mat_test2 );
  zMatSubDyad( mat_test2, vec_test1, vec_test2 );
  zMatSub( mat_test2, ans3, error );
  zAssert( zMatSubDyad, zMatIsTiny( error ) );
  zMatCopy( mat_test1, mat_test2 );
  zMatCatDyad( mat_test2,-2, vec_test1, vec_test2 );
  zMatSub( mat_test2, ans4, error );
  zAssert( zMatCatDyad, zMatIsTiny( error ) );
  zVecFreeAO( 2, vec_test1, vec_test2 );
  zMatFreeAO( 7, mat_test1, mat_test2, ans1, ans2, ans3, ans4, error );
}

void assert_quad(void)
{
  const int rowsize = MAT_ROW_SIZE;
  const int colsize = MAT_COL_SIZE;
  zMat a, w1, w2, tmp1, tmp2, q1, q2, qt1, qt2, e1, e2;
  zVec wv1, wv2;

  a = zMatAlloc( rowsize, colsize );
  wv1 = zVecAlloc( colsize );
  wv2 = zVecAlloc( rowsize );
  w1 = zMatAllocSqr( colsize );
  w2 = zMatAllocSqr( rowsize );
  tmp1 = zMatAlloc( colsize, rowsize );
  tmp2 = zMatAlloc( rowsize, colsize );
  q1 = zMatAllocSqr( rowsize );
  q2 = zMatAllocSqr( colsize );
  qt1 = zMatAllocSqr( rowsize );
  qt2 = zMatAllocSqr( colsize );
  e1 = zMatAllocSqr( rowsize );
  e2 = zMatAllocSqr( colsize );

  zMatRandUniform( a, -10, 10 );
  zVecRandUniform( wv1, 0, 10 );
  zVecRandUniform( wv2, 0, 10 );
  zMatDiag( w1, wv1 );
  zMatDiag( w2, wv2 );

  zMatQuad( a, wv1, q1 );
  zMulMatMatT( w1, a, tmp1 );
  zMulMatMat( a, tmp1, qt1 );
  zMatSub( q1, qt1, e1 );
  zAssert( zMatQuad, zMatIsTiny( e1 ) );

  zMatTQuad( a, wv2, q2 );
  zMulMatMat( w2, a, tmp2 );
  zMulMatTMat( a, tmp2, qt2 );
  zMatSub( q2, qt2, e2 );
  zAssert( zMatTQuad, zMatIsTiny( e2 ) );

  zMatFreeAO( 11, a, w1, w2, tmp1, tmp2, q1, q2, qt1, qt2, e1, e2 );
  zVecFreeAO( 2, wv1, wv2 );
}

int main(void)
{
  zRandInit();
  assert_get_put();
  assert_arith();
  assert_transpose();
  assert_mul_mat_vec();
  assert_dyad();
  assert_quad();

  return EXIT_SUCCESS;
}
