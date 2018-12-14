#include <zm/zm_vec.h>

int main(void)
{
  zVec v;

  v = zVecAlloc( 10 );
  zVecWrite( v );
  printf( "vector %s NaN.\n", zVecIsNan(v) ? "includes" : "doesn't include" );
  zVecSetElem( v, 5, NAN );
  zVecWrite( v );
  printf( "vector %s NaN.\n", zVecIsNan(v) ? "includes" : "doesn't include" );
  zVecSetElem( v, 5, HUGE_VAL );
  zVecWrite( v );
  printf( "vector %s NaN.\n", zVecIsNan(v) ? "includes" : "doesn't include" );
  zVecFree( v );
  return 0;
}
