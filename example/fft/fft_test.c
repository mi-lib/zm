#include <zm/zm_fft.h>

#define N 100

double f(double t)
{
  static double c[] = {
    10, 3, 4, 2, 1, 0.5, 0.3, 0.8, 1.2,
  };
  register int i;
  double val = 0;

  /* dummy */
  for( i=0; i<9; i++ )
    val += c[i]*(cos(2*zPI*i*t/N)-sin(2*i*zPI*t/N));
  return val;
}

int main(void)
{
  double data[N];
  zComplex res[N];
  zComplex inv[N];
  register int i;

  for( i=0; i<N; i++ )
    data[i] = f( i );

  zFFT( data, N, res );
  zFFTInv( res, N, inv );
  for( i=0; i<N; i++ ){
    printf( "%f %f %f %f %f\n", data[i], res[i].re, res[i].im, inv[i].re, inv[i].im );
  }
  return 0;
}
