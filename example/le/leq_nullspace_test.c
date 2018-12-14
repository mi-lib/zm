#include <zm/zm_le.h>

#define M 300
#define N 400

#define TURN 20

#define D_MAX 10.0

void check(zMat a, zVec b, zVec ans, zVec e, clock_t c1, clock_t c2, const char *key)
{
  zMulMatVec( a, ans, e );
  zVecSub( b, e, e );
  printf( "%s:%.1g %ld ", key, zVecNorm(e), (c2-c1)/1000 );
}

void test(zMat a, zVec b, zVec w, zVec w2, zVec x, zVec e, zVec aux)
{
  clock_t c1, c2;

  zVecRand( aux, -10, 10 );

  c1 = clock();
  zLESolveSR( a, b, w, w2, x );
  c2 = clock();
  check( a, b, x, e, c1, c2, "SR" );

  c1 = clock();
  zLESolveSRAux( a, b, w, w2, x, aux );
  c2 = clock();
  check( a, b, x, e, c1, c2, "SR(aux)" );

  c1 = clock();
  zLESolveMP( a, b, w, w2, x );
  c2 = clock();
  check( a, b, x, e, c1, c2, "MP" );

  c1 = clock();
  zLESolveMPAux( a, b, w, w2, x, aux );
  c2 = clock();
  check( a, b, x, e, c1, c2, "MP(aux)" );

  printf( "\n" );
}

int main(int argc, char *argv[])
{
  register int i, j;
  zMat a;
  zVec b, w, x, e, w2, aux;
  int turn;

  zRandInit();
  a = zMatAlloc( M, N );
  b = zVecAlloc( M );
  e= zVecAlloc( M );
  w = zVecAlloc( N );
  w2= zVecAlloc( M );
  x = zVecAlloc( N );
  aux = zVecAlloc( N );

  for( i=0; i<M; i++ ){
    for( j=0; j<N; j++ )
      zMatSetElem( a, i, j, zRandF(-D_MAX,D_MAX) );
    zVecSetElem( b, i, zRandF(-D_MAX,D_MAX) );
    zVecSetElem( w2, i, 1.0e5 );
  }
  for( j=0; j<N; j++ )
    zVecSetElem( w, j, 1.0 );

  turn = argc > 1 ? atoi(argv[1]) : TURN;
  while( turn-- > 0 ) test( a, b, w, w2, x, e, aux );

  zMatFree( a );
  zVecFreeAO( 6, b, e, w, w2, x, aux );
  return 0;
}
