#include <zm/zm_vec.h>

#define DIM 3
#define N 5

int main(void)
{
  zVec s, v[N];
  register int i;
  double val;

  s = zVecAlloc(DIM);
  for( i=0; i<N; i++ ){
    val = i + 1;
    v[i] = zVecCreateList( DIM, val, val, val );
    printf( "v_%d: ", i ); zVecWrite( v[i] );
  }
  zVecLS( s, N, 1.0/1.0, v[0], 1.0/2.0, v[1], 1.0/3.0, v[2], 1.0/4.0, v[3], 1.0/5.0, v[4] );
  zVecWrite( s );
  return 0;
}
