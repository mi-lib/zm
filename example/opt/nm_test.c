#include <zm/zm_opt.h>

#define N 5

double eval(zVec v, void *dummy)
{
  register int i;
  double result = 0;

  for( i=0; i<zVecSizeNC(v); i++ )
    result += zSqr( zVecElem(v,i) - 10*i );
  return result;
}

int main(void)
{
  zOptNM opt;
  zVec x;
  int i;

  zOptNMCreate( &opt, N, eval );
  x = zVecAlloc( N );
  printf( "initial value:  " ); zVecWrite( x );
  zOptNMSolve( &opt, x, NULL, zTOL, 0, NULL );
  printf( "solution value: " ); zVecWrite( x );

  for( i=0; i<N; i++ )
    zVecSetElem( x, i, -i );
  printf( "initial value:  " ); zVecWrite( x );
  zOptNMSolve( &opt, x, NULL, zTOL, 0, NULL );
  printf( "solution value: " ); zVecWrite( x );

  for( i=0; i<N; i++ )
    zVecSetElem( x, i, N-i );
  printf( "initial value:  " ); zVecWrite( x );
  zOptNMSolve( &opt, x, NULL, zTOL, 0, NULL );
  printf( "solution value: " ); zVecWrite( x );

  zVecRand( x, -100, 100 );
  printf( "initial value:  " ); zVecWrite( x );
  zOptNMSolve( &opt, x, NULL, zTOL, 0, NULL );
  printf( "solution value: " ); zVecWrite( x );

  zVecFree( x );
  zOptNMDestroy( &opt );
  return 0;
}
