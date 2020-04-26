#include <zm/zm_opt.h>

#define N 5

double testfunc(zVec x, void *dummy)
{
  register int i;
  double result = 0;

  for( i=0; i<zVecSizeNC(x); i++ )
    result += zSqr( zVecElem(x,i) - 10*i );
  return result;
}

void test(zVec x)
{
  double val;

  printf( "initial value:  " ); zVecPrint( x );
  zOptSolveNM( testfunc, NULL, NULL, NULL, 0, zTOL, x, &val );
  zVecTouchup( x );
  printf( "solution value: " ); zVecPrint( x );
}

int main(void)
{
  zVec x;
  int i;

  x = zVecAlloc( N );
  test( x );

  for( i=0; i<N; i++ ) zVecSetElem( x, i, -i );
  test( x );

  for( i=0; i<N; i++ ) zVecSetElem( x, i, N-i );
  test( x );

  zVecRandUniform( x, -100, 100 );
  test( x );

  zVecFree( x );
  return 0;
}
