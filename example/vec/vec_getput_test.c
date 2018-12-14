#include <zm/zm_vec.h>

int main(void)
{
  zVec v, pv;
  double val[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  int s = 3, d = 8;

  v = zVecCloneArray( val, sizeof(val)/sizeof(double) );
  pv = zVecAlloc( 4 );

  zVecWrite( v );
  printf( "get %d-%d\n", s, s+zVecSize(pv) );
  zVecGet( v, s, pv );
  zVecWrite( pv );
  zVecPut( v, d, pv );
  printf( "put %d-%d\n", d, d+zVecSize(pv) );
  zVecWrite( v );
  zVecWrite( pv );
  return 0;
}
