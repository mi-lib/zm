#include <zm/zm_complex.h>

#define N 1000
#define R 3.2

int main(void)
{
  zComplex c1, c2, c3;
  double theta;
  register int i;

  zComplexCopy( Z_ZEROCOMPLEX, &c3 );
  for( i=0; i<=N; i++ ){
    theta = 4 * zPI * i / N;
    zComplexPolar( &c1, 1.0, theta );
    zComplexPow( &c1, R, &c2 );
    zComplexPowRef( &c1, R, &c3, &c3 );
    printf( "%f %f %f %f %f %f %f\n", theta, c1.re, c1.im, c2.re, c2.im, c3.re, c3.im );
  }
  return 0;
}
