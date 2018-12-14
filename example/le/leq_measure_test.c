#include <zm/zm_le.h>

#define M 150
#define N 180

#define TURN 20

#define D_MAX 10.0

void test(zMat a, zVec b, zVec w, zVec w2, zVec x, zVec _b)
{
  register int i, j;
  clock_t c1, c2;

#if 1
  for( i=0; i<M; i++ ){
    for( j=0; j<N; j++ )
      zMatSetElem( a, i, j, zRandF(-D_MAX,D_MAX) );
    zVecSetElem( b, i, zRandF(-D_MAX,D_MAX) );
    zVecSetElem( w2, i, 1000000.0 );
  }
#else /* ill-posed case */
  int n;
  double e;

  for( j=0; j<N; j++ )
    zMatSetElem( a, 0, j, zRandF(-D_MAX,D_MAX) );
  zVecSetElem( b, 0, zRandF(-D_MAX,D_MAX) );
  zVecSetElem( w2, 0, 1.0e5 );
  for( i=1; i<M; i++ ){
    zRawVecCopy( zMatRowArray(a,0), zMatRowArray(a,i), N );
    zVecSetElem( b, i, zVecElem(b,0) );
    zVecSetElem( w2, i, zVecElem(w2,0) );
  }
  n = 0.8 * zMin( M, N );
  for( i=0; i<n; i++ ){
    e = zRandF(-0.01,0.01);
    zMatElem(a,i,i) += e;
    zVecElem(b,i  ) += e;
  }
#endif

  for( j=0; j<N; j++ )
    zVecSetElem( w, j, 1.0 );

  /*
  c1 = clock();
  zLESolveNormMin( a, b, w, x );
  c2 = clock();
  zMulMatVec( a, x, _b );
  zVecSub( b, _b, _b );
  printf( "NM:%g %ld ", zVecNorm( _b ), (c2-c1)/1000 );
  */

  c1 = clock();
  zLESolveSR( a, b, w, w2, x );
  c2 = clock();
  zMulMatVec( a, x, _b );
  zVecSub( b, _b, _b );
  printf( "SR:%g %ld ", zVecNorm( _b ), (c2-c1)/1000 );

  c1 = clock();
  zLESolveMP( a, b, w, w2, x );
  c2 = clock();
  zMulMatVec( a, x, _b );
  zVecSub( b, _b, _b );
  printf( "MP(LQ):%g %ld ", zVecNorm( _b ), (c2-c1)/1000 );

  c1 = clock();
  zLESolveMP_LU( a, b, w, w2, x );
  c2 = clock();
  zMulMatVec( a, x, _b );
  zVecSub( b, _b, _b );
  printf( "MP(LU):%g %ld ", zVecNorm( _b ), (c2-c1)/1000 );

  c1 = clock();
  zLESolveMP_SVD( a, b, x );
  c2 = clock();
  zMulMatVec( a, x, _b );
  zVecSub( b, _b, _b );
  printf( "MP(SVD):%g %ld ", zVecNorm( _b ), (c2-c1)/1000 );

  printf( "\n" );
}

int main(int argc, char *argv[])
{
  zMat a;
  zVec b, w, x, _b, w2;
  int turn;

  zRandInit();
  a = zMatAlloc( M, N );
  b = zVecAlloc( M );
  _b= zVecAlloc( M );
  w = zVecAlloc( N );
  w2= zVecAlloc( M );
  x = zVecAlloc( N );

  turn = argc > 1 ? atoi(argv[1]) : TURN;
  while( turn-- > 0 ) test( a, b, w, w2, x, _b );

  zMatFree( a );
  zVecFree( b );
  zVecFree( _b );
  zVecFree( w );
  zVecFree( w2 );
  zVecFree( x );
  return 0;
}
