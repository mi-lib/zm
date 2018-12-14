#include <zm/zm_complex.h>

void cardano(double a, double b, double c, double d)
{
  zComplex ans[3];

  printf( "%f x^3 + %f x^2 + %f x + %f = 0\n", a, b, c, d );
  if( !zCardano( a, b, c, d, ans ) ) return;
  printf( "answer->[0] " ); zComplexWrite( &ans[0] ); zEndl();
  printf( "        [1] " ); zComplexWrite( &ans[1] ); zEndl();
  printf( "        [2] " ); zComplexWrite( &ans[2] ); zEndl();
}

int main(void)
{
  cardano( 1, -6, 11, -6 );
  cardano( 1, -2, -1,  2 );
  cardano( 1, 0, 0, 0 );
  cardano( 1, 0, 0, 1 );
  cardano( 1, 0, 0,-1 );
  cardano( 0, 0, 0, 0 );
  return 0;
}
