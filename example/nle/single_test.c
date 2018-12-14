#include <zm/zm_nle.h>

double f(double x, void *priv)
{
  /*
  return (x+3)*(x+1)*(x-1)*(x-3);
  */
  return exp(-x)-x;
}

int main(void)
{
  double x;

  printf( "binary section method\n" );
  x = zNLE_Bisec( f, -10, 10, NULL, zTOL, 0 );
  printf( "f(%f) = %g\n", x, f(x,NULL) );

  printf( "Secant method\n" );
  x = zNLE_Secant( f, -10, 10, NULL, zTOL, 0 );
  printf( "f(%f) = %g\n", x, f(x,NULL) );

  printf( "Regula-Falsi method\n" );
  x = zNLE_RF( f, -10, 10, NULL, zTOL, 0 );
  printf( "f(%f) = %g\n", x, f(x,NULL) );

  printf( "Van Wijngaarden-Dekker-Brent method\n" );
  x = zNLE_VDB( f, -10, 10, NULL, zTOL, 0 );
  printf( "f(%f) = %g\n", x, f(x,NULL) );
  return 0;
}
