#include <zm/zm_le.h>

#define M 300
#define N 400

#define TURN 20

#define D_MAX 10.0

void check(zMat a, zVec b, zVec ans, zVec e, const char *key)
{
  zMulMatVec( a, ans, e );
  zVecSub( b, e, e );
  printf( "%s:%.1g ", key, zVecNorm(e) );
}

void test(zMat a, zVec b, zVec w, zVec w2, zVec e, zVec aux)
{
  zMat m;
  zVec x1, x2, v, s, bb;
  zIndex idx;

  x1 = zVecAlloc( N );
  x2 = zVecAlloc( N );
  bb = zVecAlloc( zVecSizeNC(b) );
  zLEAllocWork( &m, &v, &s, &idx, N );
  zVecRand( aux, -10, 10 );

  zLESolveSRAux( a, b, w, w2, x1, aux );
  check( a, b, x1, e, "SR(aux)" );
  zLESolveSRAuxDST( a, b, w, w2, x2, aux, m, v, idx, s, bb );
  check( a, b, x2, e, "SR(aux,DST)" );
  zVecSub( x1, x2, s );
  printf( "error = %.10g\n", zVecNorm(s) );

  zLEFreeWork( m, v, s, idx );
}

int main(int argc, char *argv[])
{
  register int i, j;
  zMat a;
  zVec b, w, e, w2, aux;
  int turn;

  zRandInit();
  a = zMatAlloc( M, N );
  b = zVecAlloc( M );
  e= zVecAlloc( M );
  w = zVecAlloc( N );
  w2= zVecAlloc( M );
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
  while( turn-- > 0 ) test( a, b, w, w2, e, aux );

  zMatFree( a );
  zVecFreeAO( 5, b, e, w, w2, aux );
  return 0;
}
