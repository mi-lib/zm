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

void test(zOptNM *opt, zVec x)
{
  printf( "initial value:  " ); zVecPrint( x );
  zOptNMSolve( opt, x, NULL, zTOL, 0, NULL );
  zVecTouchup( x );
  printf( "solution value: " ); zVecPrint( x );
}

int main(void)
{
  zOptNM opt;
  zVec x;
  int i;

  zOptNMCreate( &opt, N, eval );
  x = zVecAlloc( N );
  test( &opt, x );

  for( i=0; i<N; i++ ) zVecSetElem( x, i, -i );
  test( &opt, x );

  for( i=0; i<N; i++ ) zVecSetElem( x, i, N-i );
  test( &opt, x );

  zVecRandUniform( x, -100, 100 );
  test( &opt, x );

  zVecFree( x );
  zOptNMDestroy( &opt );
  return 0;
}
