#include <zm/zm_opt.h>

#if 0
bool _zQPSolveASM(zMat q, zVec c, zMat a, zVec b, zVec ans, double *cost)
{

}
#endif

int main(void)
{
  double qarray[] = {
    1, 2,-4, 5
  };
  double aarray[] = {
    1, 2
  };
  zMat q, a;
  zVec c, b, x;
  double cost;

  q = zMatAllocSqr( 2 );
  a = zMatAlloc( 1, 2 );
  c = zVecAlloc( 2 );
  b = zVecAlloc( 1 );
  x = zVecAlloc( 2 );

  zMatCopyArray( qarray, 2, 2, q );
  zMatCopyArray( aarray, 1, 2, a );
  zVecSetElem( b, 0, 3 );
  zMatPrint( q );
  zVecPrint( c );
  zMatPrint( a );
  zVecPrint( b );
  _zQPSolveASM( q, c, a, b, x, &cost );
  zVecPrint( x );

  zMatFree( q );
  zMatFree( a );
  zVecFree( c );
  zVecFree( b );
  zVecFree( x );
  return 0;
}
