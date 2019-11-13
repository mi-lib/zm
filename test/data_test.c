#include <zm/zm_data.h>

#define FFT_TEST_N 100
#define FFT_TEST_DIM 9

double func(double t, double sc[], double cc[], int dim)
{
  register int i;
  double s, c, val = 0;

  for( i=0; i<dim; i++ ){
    zSinCos( 2*zPI*i*t/FFT_TEST_N, &s, &c );
    val += sc[i]*s + cc[i]*c;
  }
  return val;
}

void assert_fft(void)
{
  double sc[FFT_TEST_DIM];
  double cc[FFT_TEST_DIM];
  double data[FFT_TEST_N];
  zComplex res[FFT_TEST_N];
  double inv[FFT_TEST_N];
  register int i;
  bool result = true;

  for( i=0; i<FFT_TEST_DIM; i++ ){
    sc[i] = zRandF(-10,10);
    cc[i] = zRandF(-10,10);
  }
  sc[0] = 0; /* because sin 0 = 0 */
  for( i=0; i<FFT_TEST_N; i++ ) data[i] = func( i, sc, cc, FFT_TEST_DIM );

  zFFT( data, FFT_TEST_N, res );
  zFFTInv( res, FFT_TEST_N, inv );
  for( i=0; i<FFT_TEST_N; i++ )
    if( !zIsEqual( data[i], inv[i], zTOL ) ) result = false;
  for( i=0; i<FFT_TEST_DIM; i++ )
    if( !zIsEqual( res[i].re, cc[i], zTOL ) ||
        !zIsEqual( res[i].im, sc[i], zTOL ) ) result = false;
  zAssert( zFFT + zFFTInv, result );
}

int main(int argc, char *argv[])
{
  zRandInit();
  assert_fft();
  return EXIT_SUCCESS;
}
