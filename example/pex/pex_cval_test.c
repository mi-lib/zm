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
    zComplexCreatePolar( &arg, 1.0, t);
    zPexCVal( p, &arg, &c );
    printf( "%.10g %.10g ", arg.re, arg.im );
    printf( "%.10g %.10g\n", c.re, c.im );
  }
  zPexFree( p );
  return 0;
}
