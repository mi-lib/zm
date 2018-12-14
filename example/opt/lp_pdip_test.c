/* Primal-dual Interior-point method */
#include <zm/zm_opt.h>

#define TEST 7

int main(void)
{
#if TEST == 1
  double a_arr[] = {
    3.0, 1.0, 1.0, 0.0, 0.0,
    2.5, 2.0, 0.0, 1.0, 0.0,
    1.0, 2.0, 0.0, 0.0, 1.0,
  };
  double b_arr[] = {
    9.0, 12.5, 8.0,
  };
  double c_arr[] = {
    -3.0, -2.0, 0.0, 0.0, 0.0,
  };
  int row = 3, col = 5;
  /* answer: [2 3] */
#elif TEST == 2
  double a_arr[] = {
    1.0, 2.0, 1.0, 0.0, 0.0,
    3.0, 4.0, 0.0, 1.0, 0.0,
    3.0, 1.0, 0.0, 0.0, 1.0,
  };
  double b_arr[] = {
    800.0, 1800.0, 1500.0,
  };
  double c_arr[] = {
    -20.0, -30.0, 0.0, 0.0, 0.0,
  };
  int row = 3, col = 5;
  /* answer: [200 300] */
#elif TEST == 3
  double a_arr[] = {
   -1.0, 2.0,-2.0, 1.0, 0.0,
    2.0, 3.0, 1.0, 0.0, 1.0,
  };
  double b_arr[] = {
   -8.0, 5.0,
  };
  double c_arr[] = {
    5.0, 3.0,-1.0, 0.0, 0.0,
  };
  int row = 2, col = 5;
  /* answer: [0 0 5] */
#elif TEST == 4
  double a_arr[] = {
    1.0,-1.0, 1.0, 1.0, 0.0, 0.0,
    3.0, 2.0, 4.0, 0.0, 1.0, 0.0,
    3.0, 2.0, 0.0, 0.0, 0.0, 1.0,
  };
  double b_arr[] = {
    20.0, 42.0, 30.0,
  };
  double c_arr[] = {
    -5.0,-4.0,-6.0, 0.0, 0.0, 0.0,
  };
  int row = 3, col = 6;
  /* answer: [0 15 3] */
#elif TEST == 5
  double a_arr[] = {
    1.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 2.0, 0.0,-7.0, 0.0, 1.0, 0.0,
    0.0,-1.0, 1.0,-2.0, 0.0, 0.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
  };
  double b_arr[] = {
    740.0, 0.0, -0.5, 9.0,
  };
  double c_arr[] = {
   -1.0,-1.0,-3.0, 0.5, 0.0, 0.0, 0.0,
  };
  int row = 4, col = 7;
  /* answer: [0 3.325 4.725 0.95] */
#elif TEST == 6
  double a_arr[] = {
    1.0, 6.0,-1.0, 0.0,
    0.0,-3.0, 4.0, 1.0,
  };
  double b_arr[] = {
    2.0, 8.0,
  };
  double c_arr[] = {
    0.0,-2.0, 4.0, 0.0,
  };
  int row = 2, col = 4;
  /* answer: [0 1/3 0 9] */
#else
  double a_arr[] = {
    1.0, 4.0, 3.0, 1.0, 0.0,
   -1.0, 2.0,-3.0, 0.0, 1.0,
  };
  double b_arr[] = {
    12.0, 4.0,
  };
  double c_arr[] = {
    1.0, -3.0, 1.0, 0.0, 0.0,
  };
  int row = 2, col = 5;
  /* answer: [4/3 8/3] */
#endif

  zMat a;
  zVec b, c, x;
  double cost;

  a = zMatCloneArray( a_arr, row, col );
  b = zVecCloneArray( b_arr, row );
  c = zVecCloneArray( c_arr, col );
  x = zVecAlloc( col );

  printf( "Primal-dual Interior-point method\n" );
  printf( "minimize c^T x subject to Ax<=b and x>=0, where:\n" );
  printf( "c: " ); zVecWrite( c );
  printf( "A: " ); zMatWrite( a );
  printf( "b: " ); zVecWrite( b );
  printf( "<result>\n" );
  zLPSolvePDIP_PC( a, b, c, x, &cost );
  zVecWrite( x );
  printf( "%f\n", cost );
  return 0;
}
