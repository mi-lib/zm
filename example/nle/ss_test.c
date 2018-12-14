#include <zm/zm_nle.h>

zVec f(zVec x, zVec y, void *util)
{
  zVecSetElem( y, 0, sqrt(zVecElem(x,0)) );
  zVecSetElem( y, 1, cos(zVecElem(x,0)) );
  zVecSetElem( y, 2, 2.0/(zSqr(zVecElem(x,0))+1.0)+1 );
  zVecSetElem( y, 3, 0.5*(zVecElem(x,0)+3) );
  zVecSetElem( y, 4, (1-zSqr(zVecElem(x,0)))/(1+sqrt(1+zSqr(zVecElem(x,0)))) );
  zVecSetElem( y, 5, 1.0/(1.0+zVecElem(x,0)) );
  return y;
}

int main(void)
{
  zVec x;

  x = zVecAlloc( 6 );
  zVecSetAll( x, 2.0 );
  zSSSolve( f, x, NULL, 0 );
  zVecWrite( x );
  printf( ">> assertion <<\n" );
  f( x, x, NULL );
  zVecWrite( x );
  zVecFree( x );
  return 0;
}
