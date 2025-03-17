#include <zm/zm_vec.h>

#define N 20

int main(void)
{
  int i;
  zVec v;
  zIndex idx;

  zRandInit();
  v = zVecAlloc(N);
  idx = zIndexCreate(N);
  for( i=0; i<N; i++ )
    zVecSetElem( v, i, zRandF(-100,100) );
  zVecPrint( v );
  zIndexPrint( idx );
  zVecSort( v, idx );
  printf( "+++ sorted. +++\n" );
  zIndexPrint( idx );
  for( i=0; i<N; i++ )
    printf( "%g ", zVecElem(v,zIndexElem(idx,i)) );
  zEndl();
  zVecFree( v );
  zIndexFree( idx );
  return 0;
}
