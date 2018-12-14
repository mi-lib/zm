#include <zm/zm_vec.h>

int main(void)
{
  zVec v1, v2, v;
  double arr[] = { 1, 2, 3, 4 };
  FILE *fp;

  fp = fopen( "vec.dat", "r" );
  v1 = zVecFRead( fp );
  v2 = zVecFRead( fp );
  printf( "v1: " ); zVecWrite( v1 );
  printf( "v2: " ); zVecWrite( v2 );
  fclose( fp );

  printf( "v1 clone: " ); zVecWrite( ( v = zVecClone( v1 ) ) );
  printf( "v1 arr: " );   zVecWrite( zVecCopyArray( arr, 4, v1 ) );

  /* indirect test */
  printf( "indirect test\n" );
  printf( "v1: " ); zVecWrite( v1 );
  printf( "v2: " ); zVecWrite( v2 );
  printf( "v1 + v2: " ); zVecWrite( zVecAdd( v1, v2, v ) );
  printf( "v1 - v2: " ); zVecWrite( zVecSub( v1, v2, v ) );
  printf( "-v1    : " ); zVecWrite( zVecRev( v1, v ) );
  printf( "v1 * 2 : " ); zVecWrite( zVecMul( v1, 2, v ) );
  printf( "v1 / 2 : " ); zVecWrite( zVecDiv( v1, 2, v ) );
  printf( "v1 * v2: " ); zVecWrite( zVecAmp( v1, v2, v ) );
  printf( "v1 / v2: " ); zVecWrite( zVecDem( v1, v2, v ) );
  printf( "v1^T v2: %f\n", zVecInnerProd( v1, v2 ) );
  printf( "v = {PI,...}: " ); zVecWrite( zVecSetAll( v, zPI ) );
  printf( "v = {10,11,12,...}: " ); zVecWrite( zVecLinSpace( v, 10, 13 ) );

  /* direct test */
  printf( "direct test\n" );
  printf( "v1: " ); zVecWrite( v1 );
  printf( "v2: " ); zVecWrite( v2 );
  printf( "v1 + v2: " ); zVecWrite( zVecAddDRC( v1, v2 ) );
  printf( "v1 - v2: " ); zVecWrite( zVecSubDRC( v1, v2 ) );
  printf( "-v1    : " ); zVecWrite( zVecRevDRC( v1 ) );
  printf( "v1 * 2 : " ); zVecWrite( zVecMulDRC( v1, 2 ) );
  printf( "v1 / 2 : " ); zVecWrite( zVecDivDRC( v1, 2 ) );
  printf( "v1 * v2: " ); zVecWrite( zVecAmpDRC( v1, v2 ) );
  printf( "v1 / v2: " ); zVecWrite( zVecDemDRC( v1, v2 ) );

  zVecFree( v1 );
  zVecFree( v2 );
  zVecFree( v );
  return 0;
}
