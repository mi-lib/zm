#include <zm/zm_le.h>

int main(void)
{
  zMat m;
  zVec v, mv;

  m = zMatAlloc( 5, 6 );
  v = zVecAlloc( 4 );
  mv = zVecAlloc( 5 );

  zMatSetSize( m, 3, 3 );
  zVecSetSize( v, 3 );
  zVecSetSize( mv, 3 );
  zRandInit();
  zMatSetElemList( m,
    2.0, 1.0, 4.0,
   -1.0, 2.0,-3.0,
    3.0, 5.0,-2.0
  );
  zVecSetElemList( mv,
   -1.0, 1.0, 2.0
  );
  zMulMatVec( m, mv, v );
  zMatPrint( m );
  zVecPrint( v );
  printf( "The answer should be " );
  zVecPrint( mv );
  zLESolveGauss( m, v, mv );
  printf( "Answer " );
  zVecPrint( mv );

  zMatFree( m );
  zVecFree( v );
  zVecFree( mv );
  return 0;
}
