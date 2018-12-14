#include <zm/zm_opt.h>

double eval(double x, void *dummy){
  /* return cos(x); */
  return x*(x*(x*(x-12)+47)-60);       /* x=1.25 */
  return fabs( (x-1.6)*(x-10)*(x+8) ); /* x=1.6 */
  return (x*x+x+1)*(x-3)*(x-3);        /* x=3 */
  return x*x*(x-1)*(x-4)*(x-5)*(x-5);  /* symmetric about x=2.5 */
}

int main(void)
{
  double x;

  x = zOptLineGSEC( eval, -10, 20, NULL, 0 );
  printf( "golden section method : x* = %.16f ( evaluation = %.16f )\n", x, eval( x, NULL ) );
  x = zOptLineBisec( eval, -10, 20, NULL, 0 );
  printf( "bisection method :      x* = %.16f ( evaluation = %.16f )\n", x, eval( x, NULL ) );
  x = zOptLineBrent( eval, -10, 20, NULL, 0 );
  printf( "Brent's method :        x* = %.16f ( evaluation = %.16f )\n", x, eval( x, NULL ) );
  return 0;
}
