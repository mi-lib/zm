#include <zm/zm_opt.h>

int main(void)
{
  double a_arr[] = {
    3.0, 1.0,
    4.0, 2.0,
    1.0, 2.0,
  };
  double c_arr[] = {
   -3.0, -2.0,
  };
  int row = 3, col = 2;

  zMat a, as;
  zVec c, cs, x, xs;

  a = zMatCloneArray( a_arr, row, col );
  c = zVecCloneArray( c_arr, col );
  x = zVecAlloc( col );

  printf( "<original matrix/vector>\n" );
  printf( "A: " ); zMatWrite( a );
  printf( "c: " ); zVecWrite( c );

  printf( "<inequality case>\n" );
  zLPIneq2Std( a, c, x, &as, &cs, &xs );
  printf( "A: " ); zMatWrite( as );
  printf( "c: " ); zVecWrite( cs );
  zMatFree( as );
  zVecFree( cs );
  zVecFree( xs );

  printf( "<unbound case>\n" );
  zLPUnb2Std( a, c, x, &as, &cs, &xs );
  printf( "A: " ); zMatWrite( as );
  printf( "c: " ); zVecWrite( cs );
  zMatFree( as );
  zVecFree( cs );
  zVecFree( xs );

  zMatFree( a );
  zVecFree( c );
  zVecFree( x );
  return 0;
}
