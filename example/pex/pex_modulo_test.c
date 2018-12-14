#include <zm/zm_pex.h>

#define N 10

int main(void)
{
  double a;
  zPex p1, p2;
  zComplex arg, c;

  zRandInit();
  a = (double)zRandI(-10,10);
  p1 = zPexAlloc( N );
  zVecLinSpace( p1, 1, N+1 );
  p2 = zPexAlloc( N );
  zPexModulo( p1, a, p2 );
  zPexExprX( p1 );
  printf( "modulo: x%c%f\n", a>0?'+':'-', fabs(a) );
  zVecWrite( p2 );

  zComplexCreate( &arg, -a, 0 );
  zPexCVal( p1, &arg, &c );
  printf( "check: P(-a) = " );
  zComplexWrite( &c );
  printf( " = %g\n", zPexCoeff(p2,0) );

  zPexFree( p1 );
  zPexFree( p2 );
  return 0;
}
