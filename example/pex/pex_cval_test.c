#include <zm/zm_pex.h>
#include <zm/zm_complex.h>

#define N 10

int main(void)
{
  register int i;
  double t;
  zPex p;
  zComplex arg, c;

  p = zPexAlloc(N);
  zPexSetCoeff( p, N, 1.0 );
  zPexSetCoeff( p, 0,-1.0 );
  for( i=0; i<=N; i++ ){
    t = 2 * zPI * i / N;
    zComplexPolar( &arg, 1.0, t );
    zPexCVal( p, &arg, &c );
    printf( "%.10g %.10g ", zComplexRe(&arg), zComplexIm(&arg) );
    printf( "%.10g %.10g\n", zComplexRe(&c), zComplexIm(&c) );
  }
  zPexFree( p );
  return 0;
}
