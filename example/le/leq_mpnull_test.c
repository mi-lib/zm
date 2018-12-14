#include <zm/zm_le.h>

#define ROW 200
#define COL 300

int main(int argc, char *argv[])
{
  zMat a, mn;
  zVec b, bn, y, yn, ans;

  a = zMatAlloc( ROW, COL );
  b = zVecAlloc( ROW );
  bn = zVecAlloc( ROW );
  y = zVecAlloc( COL );
  yn = zVecAlloc( COL );
  ans = zVecAlloc( COL );
  mn = zMatAllocSqr( COL );

  zRandInit();
  zMatRandUniform( a, -10, 10 );
  zVecRandUniform( b, -10, 10 );
  zVecRandUniform( y, -10, 10 );
  zLESolveMPNull( a, b, NULL, NULL, ans, mn );
  /* null-space */
  zMulMatVec( mn, y, yn );
  zVecAddDRC( ans, yn );
  zMulMatVec( a, ans, bn );
  printf( "%g %g\n", zVecNorm( yn ), zVecDist( b, bn ) );

  zMatFreeAO( 2, a, mn );
  zVecFreeAO( 5, b, bn, y, yn, ans );
  return 0;
}
