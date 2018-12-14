#include <zm/zm_le.h>

void test1(void)
{
  zMat a;
  zVec b, w, x, y;

  printf( "+++ test 1 (norm min.) +++\n" );
  a = zMatCreateList( 1, 3, 1.0, 1.0, 1.0 );
  b = zVecCreateList( 1, 1.0 );
  w = zVecCreateList( 3, 0.1, 1.0, 1.0 );
  x = zVecAlloc( 3 );
  y = zVecAlloc( 1 );

  printf( "A: " ); zMatWrite( a );
  printf( "b: " ); zVecWrite( b );

  zLESolveNormMin( a, b, NULL, x );
  printf( "x(no-weight)= " ); zVecWrite( x );
  zMulMatVec( a, x, y );
  printf( "A x = " ); zVecWrite( y );

  zLESolveNormMin( a, b, w, x );
  printf( "x(with weight)= " ); zVecWrite( x );
  zMulMatVec( a, x, y );
  printf( "A x = " ); zVecWrite( y );
}

void test2(void)
{
  zMat a;
  zVec b, w, x, y;

  printf( "+++ test 2 (error min.) +++\n" );
  a = zMatCreateList( 3, 1, 1.0, 1.0, 1.0 );
  b = zVecCreateList( 3, 1.0, 2.0, 3.0 );
  w = zVecCreateList( 3, 0.1, 1.0, 1.0 );
  x = zVecAlloc( 1 );
  y = zVecAlloc( 3 );

  printf( "A: " ); zMatWrite( a );
  printf( "b: " ); zVecWrite( b );

  zLESolveErrorMin( a, b, NULL, x );
  printf( "x(no-weight)= " ); zVecWrite( x );
  zMulMatVec( a, x, y );
  printf( "A x = " ); zVecWrite( y );

  zLESolveErrorMin( a, b, w, x );
  printf( "x(with weight)= " ); zVecWrite( x );
  zMulMatVec( a, x, y );
  printf( "A x = " ); zVecWrite( y );
}

void test3(void)
{
  double a_arr[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1 };
  double b_arr[] = { 2, 3, 4 };
  double w_arr[] = { 1, 1, 1 };
  zMat a;
  zVec b, w, x;

  printf( "+++ test 3 +++\n" );
  a = zMatCloneArray( a_arr, 3, 3 );
  b = zVecCloneArray( b_arr, 3 );
  w = zVecCloneArray( w_arr, 3 );
  x = zVecAlloc( 3 );

  zMatWrite( a );
  zVecWrite( b );
  printf( "norm minimization\n" );
  zLESolveNormMin( a, b, NULL, x );
  zVecWrite( x );
  printf( "squared error minimization\n" );
  zLESolveErrorMin( a, b, NULL, x );
  zVecWrite( x );
}

int main(void)
{
  test1();
  test2();
  test3();
  return 0;
}
