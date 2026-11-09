#include <zm/zm_nle.h>

double f(double x, void *util)
{
  /*
  return (x+3)*(x+1)*(x-1)*(x-3);
  */
  return exp(-x)-x;
}

void check(const char *label, int (* method)(double (*)(double,void*),double*,double,double,void*,double,double,int))
{
  double x;
  int ret;

  printf( "%s\n", label );
  ret = method( f, &x, -10, 10, NULL, zTOL, zTOL, 0 );
  printf( " returned %d, f(%f) = %g\n", ret, x, f(x,NULL) );
}

int main(void)
{
  check( "bisection method", zNLE1Bisec );
  check( "secant method", zNLE1Secant );
  check( "Regula-Falsi method", zNLE1RegulaFalsi );
  check( "Van Wijngaarden-Dekker-Brent method", zNLE1VDB );
  return 0;
}
