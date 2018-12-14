#include <zm/zm_vec.h>

#define DIM 3
#define N 5

int main(void)
{
  zVec s, k;
  zVecRing v;
  register int i;
  double val;

  s = zVecAlloc( DIM );
  k = zVecCreateList( N, 1.0/1.0, 1.0/2.0, 1.0/3.0, 1.0/4.0, 1.0/5.0 );
  zVecRingAlloc( &v, DIM, N );
  for( i=0; i<N; i++ ){
    val = i + 1;
    zVecSetElemList( *zRingElem(&v,i), val, val, val );
    printf( "v_%d: ", i ); zVecWrite( *zRingElem(&v,i) );
  }
  printf( "k: " ); zVecWrite( k );
  zVecRingLS( s, k, &v );
  zVecWrite( s );
  zVecRingFree( &v );
  return 0;
}
