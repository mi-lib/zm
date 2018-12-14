#include <zm/zm_le.h>

int main(void)
{
  /* m1 ... r x c
     mp ... c x r
     m2 ... r x w
     ma ... c x w
     mb ... c x w
     v1 ... r x 1
     v2 ... c x 1
     v3 ... c x 1
   */
  zMat mp, m1, m2, ma, mb;
  zVec v1, v2, v3;
  int r, c, w;

  zRandInit();
  r = 400;
  c = 800;
  w = 100;
  mp = zMatAlloc( c, r );
  m1 = zMatAlloc( r, c );
  m2 = zMatAlloc( r, w );
  ma = zMatAlloc( c, w );
  mb = zMatAlloc( c, w );
  v1 = zVecAlloc( r );
  v2 = zVecAlloc( c );
  v3 = zVecAlloc( c );
  zMatRandUniform( m1, -10, 10 );
  zMatRandUniform( m2, -10, 10 );
  zMatGetCol( m2, 0, v1 );

  zMPInv( m1, mp );
  zMulMatMat( mp, m2, ma );
  zMulMPInvMatMat( m1, m2, mb );
  printf( "M1^# M2 vs M1^# M2 test: %s.\n", zMatIsEqual(ma,mb) ? "ok" : "error" );

  zLESolveMP( m1, v1, NULL, NULL, v2 );
  zMulMatVec( mp, v1, v3 );
  printf( "M1^# v  vs M1^# v  test: %s.\n", zVecIsEqual(v2,v3) ? "ok" : "error" );

  zMatGetCol( ma, 0, v3 );
  printf( "c(M1^# M2,1) vs M1^# v test: %s.\n", zVecIsEqual(v2,v3) ? "ok" : "error" );

  zMatFreeAO( 5, m1, m2, ma, mb, mp );
  zVecFreeAO( 3, v1, v2, v3 );
  return 0;
}
