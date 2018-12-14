#include <zm/zm_opt.h>

int main(void)
{
  double a_arr[] = {
    3.0, 1.0,
    2.5, 2.0,
    1.0, 2.0,
  };
  double b_arr[] = {
    9.0, 12.5, 8.0,
  };
  double c_arr[] = {
    -3.0, -2.0,
  };
  int row = 3, col = 2;
  /* answer: [2 3] */

  zMat a, as;
  zVec b, c, cs, x, xs;
  double cost;

  a = zMatCloneArray( a_arr, row, col );
  b = zVecCloneArray( b_arr, row );
  c = zVecCloneArray( c_arr, col );
  x = zVecAlloc( col );

  printf( "<original matrix/vector>\n" );
  printf( "Ax<=b, x>=0, c^T x -> min.\n" );
  printf( "A: " ); zMatWrite( a );
  printf( "b: " ); zVecWrite( b );
  printf( "c: " ); zVecWrite( c );
  zLPIneq2Std( a, c, x, &as, &cs, &xs );
  printf( "standardized form <original matrix/vector>\n" );
  printf( "Ax=b, x>=0, c^T x -> min.\n" );
  printf( "A: " ); zMatWrite( as );
  printf( "b: " ); zVecWrite( b );
  printf( "c: " ); zVecWrite( cs );

  printf( "<result (by simplex method)>\n" );
  if( !zLPSolveSimplex( as, b, cs, xs, &cost ) ){
    printf( "failed.\n" );
    return 0;
  }
  zVecGet( xs, 0, x );
  zVecWrite( x );
  printf( "cost=%f\n", cost );

  printf( "<result (by primal-dual interior-point method)>\n" );
  if( !zLPSolvePDIP_PC( as, b, cs, xs, &cost ) ){
    printf( "failed.\n" );
    return 0;
  }
  zVecGet( xs, 0, x );
  zVecWrite( x );
  printf( "cost=%f\n", cost );

  zMatFree( a );
  zVecFree( b );
  zVecFree( c );
  zVecFree( x );
  zMatFree( as );
  zVecFree( cs );
  zVecFree( xs );
  return 0;
}
