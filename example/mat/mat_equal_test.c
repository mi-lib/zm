#include <zm/zm_le.h>

void check(zMat a, zMat b)
{
  printf( ">>>check<<<\n" );
  zMatWrite( a );
  zMatWrite( b );
  zMatIsEqual( a, b ) ? printf( "equal\n" ) : printf( "unequal\n" );
}

int main(void)
{
  zMat a, b, c;

  a = zMatCreateList( 2, 2, 1.0, 1.1,-1.2, 1.3 );
  b = zMatClone( a );
  c = zMatAllocSqr( 2 );

  check( a, b );
  zMatSubDRC( a, b );
  zMatAddDRC( a, b );
  check( a, b );

  zMulInvMatMat( a, b, c );
  zMulMatMatDRC( a, c );
  check( b, c );

  zMatFree( a );
  zMatFree( b );
  return 0;
}
