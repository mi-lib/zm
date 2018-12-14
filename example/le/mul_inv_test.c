#include <zm/zm_le.h>

void test1(void)
{
  zMat m1, m2, m;
  double marray[] = {
    1, 1, 0, 1,-1,
    0, 1,-1, 3, 0,
   -2,-1, 1, 2, 3,
    1, 0, 2,-1,-1,
   -1,-2,-3,-7, 1
  };

  m1 = zMatCloneArray( marray, 5, 5 );
  m2 = zMatCloneArray( marray, 5, 5 );
  m = zMatAllocSqr( 5 );

  printf( "m1^-1 m2\n" );
  zMulInvMatMat( m1, m2, m );
  zMatWrite( m1 );
  zMatWrite( m2 );
  zMatWrite( m );

  printf( "m1 m2^-1\n" );
  zMulMatInvMat( m1, m2, m );
  zMatWrite( m1 );
  zMatWrite( m2 );
  zMatWrite( m );

  printf( "m2 = m1^-1\n" );
  zMatInv( m1, m2 );
  zMatWrite( m2 );
  zMulMatMat( m1, m2, m );
  zMatTouchup( m );
  zMatWrite( m );

  zMatFree( m1 );
  zMatFree( m2 );
  zMatFree( m );
}

void test2(void)
{
  zMat m1, m2, m3, m;

  zRandInit();
  m1 = zMatAlloc( 5, 2 );
  m2 = zMatAlloc( 5, 2 );
  m3 = zMatAlloc( 5, 2 );
  m = zMatAllocSqr( 5 );
  zMatRandUniform( m, -10, 10 );
  zMatRandUniform( m1,-10, 10 );

  zMulMatMat( m, m1, m2 );

  printf( "m^-1 m2\n" );
  zMulInvMatMat( m, m2, m3 );
  zMatWrite( m1 );
  zMatWrite( m3 );
  zMatSub( m1, m3, m2 );
  printf( "+++ error +++\n" );
  zMatWrite( m2 );

  zMatFreeAO( 4, m, m1, m2, m3 );
}

int main(void)
{
#if 0
  test1();
#else
  test2();
#endif
  return 0;
}
