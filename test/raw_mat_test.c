#include <zm/zm.h>

#define MAT_ROW_SIZE 12
#define MAT_COL_SIZE 10

void assert_get_put(void)
{
  const int rowsize = MAT_ROW_SIZE;
  const int colsize = MAT_COL_SIZE;
  double mat_test1[MAT_ROW_SIZE * MAT_COL_SIZE];
  double mat_test2[MAT_ROW_SIZE * MAT_COL_SIZE];
  double mat_test3[MAT_ROW_SIZE * MAT_COL_SIZE];
  register int i, j;
  bool result;

  zRawMatRandUniform( mat_test1, rowsize, colsize, -10, 10 );
  zRawMatRandUniform( mat_test2, rowsize, colsize, -10, 10 );
  zRawMatGet( mat_test1, rowsize, colsize, 4, 3, mat_test3, 3, 2 );
  for( result=true, i=0; i<3; i++ )
    for( j=0; j<2; j++ )
      if( mat_test1[(4+i)*colsize+j+3] != mat_test3[i*2+j] ) result = false;
  zAssert( zRawMatGet, result );
  zRawMatPut( mat_test2, rowsize, colsize, 5, 6, mat_test3, 3, 2 );
  for( result=true, i=0; i<3; i++ )
    for( j=0; j<2; j++ )
      if( mat_test2[(5+i)*colsize+j+6] != mat_test3[i*2+j] ) result = false;
  zAssert( zRawMatPut, result );
  zRawMatTGet( mat_test1, rowsize, colsize, 4, 3, mat_test3, 3, 2 );
  for( result=true, i=0; i<3; i++ )
    for( j=0; j<2; j++ )
      if( mat_test1[(4+j)*colsize+i+3] != mat_test3[i*2+j] ) result = false;
  zAssert( zRawMatTGet, result );
  zRawMatTPut( mat_test2, rowsize, colsize, 5, 6, mat_test3, 3, 2 );
  for( result=true, i=0; i<3; i++ )
    for( j=0; j<2; j++ )
      if( mat_test2[(5+j)*colsize+i+6] != mat_test3[i*2+j] ) result = false;
  zAssert( zRawMatTPut, result );
  zRawMatGetRow( mat_test1, rowsize, colsize, 3, mat_test3 );
  zAssert( zRawMatGetRow, memcmp( mat_test1+colsize*3, mat_test3, sizeof(double)*colsize ) == 0 );
  zRawMatPutRow( mat_test2, rowsize, colsize, 3, mat_test3 );
  zAssert( zRawMatPutRow, memcmp( mat_test2+colsize*3, mat_test3, sizeof(double)*colsize ) == 0 );
  zRawMatGetCol( mat_test1, rowsize, colsize, 3, mat_test3 );
  for( result=true, i=0; i<rowsize; i++ )
    if( mat_test1[i*colsize+3] != mat_test3[i] ) result = false;
  zAssert( zRawMatGetCol, result );
  zRawMatPutCol( mat_test2, rowsize, colsize, 3, mat_test3 );
  for( result=true, i=0; i<rowsize; i++ )
    if( mat_test2[i*colsize+3] != mat_test3[i] ) result = false;
  zAssert( zRawMatPutCol, result );
  zRawMatCopy( mat_test1, mat_test2, rowsize, colsize );
  zRawMatSwapRow( mat_test2, rowsize, colsize, 3, 7 );
  zAssert( zRawMatSwapRow,
    memcmp( mat_test1+colsize*3, mat_test2+colsize*7, sizeof(double)*colsize ) == 0 &&
    memcmp( mat_test1+colsize*7, mat_test2+colsize*3, sizeof(double)*colsize ) == 0 );
  zRawMatCopy( mat_test1, mat_test2, rowsize, colsize );
  zRawMatSwapCol( mat_test2, rowsize, colsize, 3, 7 );
  for( result=true, i=0; i<rowsize; i++ )
    if( mat_test1[i*colsize+3] != mat_test2[i*colsize+7] ||
        mat_test1[i*colsize+7] != mat_test2[i*colsize+3] ) result = false;
  zAssert( zRawMatSwapCol, result );
}

void assert_arith(void)
{
  const int rowsize = MAT_ROW_SIZE;
  const int colsize = MAT_COL_SIZE;
  double mat_test1[MAT_ROW_SIZE * MAT_COL_SIZE];
  double mat_test2[MAT_ROW_SIZE * MAT_COL_SIZE];
  double mat_test3[MAT_ROW_SIZE * MAT_COL_SIZE];
  double k;
  register int i, j;
  bool result;

  zRawMatRandUniform( mat_test1, rowsize, colsize, -10, 10 );
  zRawMatRandUniform( mat_test2, rowsize, colsize, -10, 10 );
  k = zRandF( -10, 10 );
  zRawMatAdd( mat_test1, mat_test2, mat_test3, rowsize, colsize );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( !zIsTiny( mat_test1[i*colsize+j]+mat_test2[i*colsize+j]-mat_test3[i*colsize+j] ) ) result = false;
  zAssert( zRawMatAdd, result );
  zRawMatSub( mat_test1, mat_test2, mat_test3, rowsize, colsize );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( !zIsTiny( mat_test1[i*colsize+j]-mat_test2[i*colsize+j]-mat_test3[i*colsize+j] ) ) result = false;
  zAssert( zRawMatSub, result );
  zRawMatRev( mat_test1, mat_test3, rowsize, colsize );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( !zIsTiny( mat_test1[i*colsize+j]+mat_test3[i*colsize+j] ) ) result = false;
  zAssert( zRawMatRev, result );
  zRawMatMul( mat_test1, k, mat_test3, rowsize, colsize );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( !zIsTiny( mat_test1[i*colsize+j]*k-mat_test3[i*colsize+j] ) ) result = false;
  zAssert( zRawMatMul, result );
  zRawMatDiv( mat_test1, k, mat_test3, rowsize, colsize );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( !zIsTiny( mat_test1[i*colsize+j]/k-mat_test3[i*colsize+j] ) ) result = false;
  zAssert( zRawMatDiv, result );
  zRawMatCat( mat_test1, k, mat_test2, mat_test3, rowsize, colsize );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( !zIsTiny( mat_test1[i*colsize+j]+k*mat_test2[i*colsize+j]-mat_test3[i*colsize+j] ) ) result = false;
  zAssert( zRawMatCat, result );

  zRawMatCopy( mat_test1, mat_test3, rowsize, colsize );
  zRawMatAddDRC( mat_test3, mat_test2, rowsize, colsize );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( !zIsTiny( mat_test1[i*colsize+j]+mat_test2[i*colsize+j]-mat_test3[i*colsize+j] ) ) result = false;
  zAssert( zRawMatAddDRC, result );
  zRawMatCopy( mat_test1, mat_test3, rowsize, colsize );
  zRawMatSubDRC( mat_test3, mat_test2, rowsize, colsize );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( !zIsTiny( mat_test1[i*colsize+j]-mat_test2[i*colsize+j]-mat_test3[i*colsize+j] ) ) result = false;
  zAssert( zRawMatSubDRC, result );
  zRawMatCopy( mat_test1, mat_test3, rowsize, colsize );
  zRawMatRevDRC( mat_test3, rowsize, colsize );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( !zIsTiny( mat_test1[i*colsize+j]+mat_test3[i*colsize+j] ) ) result = false;
  zAssert( zRawMatRevDRC, result );
  zRawMatCopy( mat_test1, mat_test3, rowsize, colsize );
  zRawMatMulDRC( mat_test3, k, rowsize, colsize );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( !zIsTiny( mat_test1[i*colsize+j]*k-mat_test3[i*colsize+j] ) ) result = false;
  zAssert( zRawMatMulDRC, result );
  zRawMatCopy( mat_test1, mat_test3, rowsize, colsize );
  zRawMatDivDRC( mat_test3, k, rowsize, colsize );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( !zIsTiny( mat_test1[i*colsize+j]/k-mat_test3[i*colsize+j] ) ) result = false;
  zAssert( zRawMatDivDRC, result );
  zRawMatCopy( mat_test1, mat_test3, rowsize, colsize );
  zRawMatCatDRC( mat_test3, k, mat_test2, rowsize, colsize );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( !zIsTiny( mat_test1[i*colsize+j]+k*mat_test2[i*colsize+j]-mat_test3[i*colsize+j] ) ) result = false;
  zAssert( zRawMatCatDRC, result );
}

void assert_transpose(void)
{
  const int rowsize = MAT_ROW_SIZE;
  const int colsize = MAT_COL_SIZE;
  double mat_test1[MAT_COL_SIZE * MAT_ROW_SIZE];
  double mat_test2[MAT_ROW_SIZE * MAT_COL_SIZE];
  double tr;
  register int i, j;
  bool result;

  zRawMatRandUniform( mat_test1, colsize, rowsize, -10, 10 );
  zRawMatT( mat_test1, mat_test2, rowsize, colsize );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( mat_test1[j*rowsize+i] != mat_test2[i*colsize+j] ) result = false;
  zAssert( zRawMatT, result );
  zRawMatCopy( mat_test1, mat_test2, colsize, rowsize );
  zRawMatTDRC( mat_test2, colsize, rowsize );
  for( result=true, i=0; i<rowsize; i++ )
    for( j=0; j<colsize; j++ )
      if( mat_test1[j*rowsize+i] != mat_test2[i*colsize+j] ) result = false;
  zAssert( zRawMatTDRC, result );
  for( tr=0, i=0; i<zMin(rowsize,colsize); i++ )
    tr += mat_test2[i*colsize+i];
  zAssert( zRawMatTrace, zIsTiny( tr - zRawMatTrace( mat_test2, rowsize, colsize ) ) );
}

void assert_mul_mat_vec(void)
{
  double mat_test1[] = {
    2, 1,-1, 2,
    1, 3, 2,-1,
   -1, 2,-3, 1,
/* transpose:
  2, 1,-1,
  1, 3, 2,
 -1, 2,-3,
  2,-1, 1,
 */
  };
  double mat_test2[] = {
    1, 1, 2,
    2,-2,-1,
   -1, 3, 2,
   -2, 1,-3,
/* if 3x4:
    1, 1, 2, 2,
   -2,-1,-1, 3,
    2,-2, 1,-3,
 */
  };
  double vec_test1[] = {
    1, 3, 4,-2,
  };
  double vec_test2[] = {
    2,-1, 3,
  };
  double mat_test3[16];
  double vec_test3[4];
  double ans1[] = { -3, 20, -9 };
  double ans2[] = { 0, 5,-13, 8 };
  double ans3[] = {
    1, -1, -5,
    7,  0,  6,
    4,-13,-13,
  };
  double ans4[] = {
    5,  2, -5,
    6,-10,  1,
   -3,  6,-12,
  };
  double ans5[] = {
    -2, 3, 2, 10,
    -1,-6, 1,  5,
   -11, 3,-7, 13,
     6, 1, 6, -2,
  };
  double error[16];

  zRawMulMatVec( mat_test1, vec_test1, 3, 4, vec_test3 );
  zRawVecSub( vec_test3, ans1, error, 3 );
  zAssert( zRawMulMatVec, zRawVecIsTiny( error, 3 ) );
  zRawMulMatTVec( mat_test1, vec_test2, 3, 4, vec_test3 );
  zRawVecSub( vec_test3, ans2, error, 4 );
  zAssert( zRawMulMatTVec, zRawVecIsTiny( error, 4 ) );
  zRawMulMatMat( mat_test1, 3, 4, mat_test2, 4, 3, mat_test3 );
  zRawMatSub( mat_test3, ans3, error, 3, 3 );
  zAssert( zRawMulMatMat, zRawMatIsTiny( error, 3, 3 ) );
  zRawMulMatMatT( mat_test1, 3, 4, mat_test2, 3, 4, mat_test3 );
  zRawMatSub( mat_test3, ans4, error, 3, 3 );
  zAssert( zRawMulMatMatT, zRawMatIsTiny( error, 3, 3 ) );
  zRawMulMatTMat( mat_test1, 3, 4, mat_test2, 3, 4, mat_test3 );
  zRawMatSub( mat_test3, ans5, error, 4, 4 );
  zAssert( zRawMulMatTMat, zRawMatIsTiny( error, 4, 4 ) );
}

void assert_dyad(void)
{
  const int rowsize = 4;
  const int colsize = 3;
  double vec_test1[] = { 1, 2, 3, 4 };
  double vec_test2[] = { 3, 2, 1 };
  double mat_test1[] = {
    3, 1,-5,
    2,-4, 2,
    1,-3, 2,
   -2, 4,-1,
  };
  double mat_test2[12];
  double ans1[] = {
    3, 2, 1,
    6, 4, 2,
    9, 6, 3,
   12, 8, 4,
  };
  double ans2[] = {
    6, 3,-4,
    8, 0, 4,
   10, 3, 5,
   10,12, 3,
  };
  double ans3[] = {
    0,-1,-6,
   -4,-8, 0,
   -8,-9,-1,
  -14,-4,-5,
  };
  double ans4[] = {
   -3, -3,-7,
  -10,-12,-2,
  -17,-15,-4,
  -26,-12,-9,
  };
  double error[12];

  zRawVecDyad( vec_test1, rowsize, vec_test2, colsize, mat_test2 );
  zRawMatSub( mat_test2, ans1, error, 3, 4 );
  zAssert( zRawVecDyad, zRawMatIsTiny( error, 3, 4 ) );
  zRawMatCopy( mat_test1, mat_test2, 4, 3 );
  zRawMatAddDyad( mat_test2, vec_test1, 4, vec_test2, 3 );
  zRawMatSub( mat_test2, ans2, error, 3, 4 );
  zAssert( zRawMatAddDyad, zRawMatIsTiny( error, 3, 4 ) );
  zRawMatCopy( mat_test1, mat_test2, 4, 3 );
  zRawMatSubDyad( mat_test2, vec_test1, 4, vec_test2, 3 );
  zRawMatSub( mat_test2, ans3, error, 3, 4 );
  zAssert( zRawMatSubDyad, zRawMatIsTiny( error, 3, 4 ) );
  zRawMatCopy( mat_test1, mat_test2, 4, 3 );
  zRawMatCatDyad( mat_test2,-2, vec_test1, 4, vec_test2, 3 );
  zRawMatSub( mat_test2, ans4, error, 3, 4 );
  zAssert( zRawMatCatDyad, zRawMatIsTiny( error, 3, 4 ) );
}

int main(void)
{
  zRandInit();
  assert_get_put();
  assert_arith();
  assert_transpose();
  assert_mul_mat_vec();
  assert_dyad();

  return EXIT_SUCCESS;
}
