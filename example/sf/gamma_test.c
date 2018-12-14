#include <zm/zm_sf.h>

int main(void)
{
  register int i;
  double x;

  for( i=-5000; i<=5000; i++ ){
    x = ((double)i) / 1000.0;
    printf( "%.10f %.10f %.10f %.10f\n", x, exp(zLnGamma(x)), zGamma(x), tgamma(x) );
  }
  return 0;
}
