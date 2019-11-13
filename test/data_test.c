#include <zm/zm_data.h>
#include <zm/zm_fourier.h>

#define FFT_TEST_N 100
#define FFT_TEST_DIM 9

void assert_fft(void)
{
  zFourier f;
  double data[FFT_TEST_N];
  zComplex res[FFT_TEST_N];
  double inv[FFT_TEST_N];
  register int i;
  bool result = true;

  zFourierAlloc( &f, FFT_TEST_DIM );
  for( i=0; i<f.n; i++ ){
    zFourierSetSinCoeff( &f, i, zRandF(-10,10) );
    zFourierSetCosCoeff( &f, i, zRandF(-10,10) );
  }
  for( i=0; i<FFT_TEST_N; i++ ) data[i] = zFourierVal( &f, 2*zPI*i/FFT_TEST_N );

  zFFT( data, FFT_TEST_N, res );
  zFFTInv( res, FFT_TEST_N, inv );
  for( i=0; i<FFT_TEST_N; i++ )
    if( !zIsEqual( data[i], inv[i], zTOL ) ) result = false;
  for( i=0; i<FFT_TEST_DIM; i++ )
    if( !zIsEqual( res[i].re, f.cc[i], zTOL ) ||
        !zIsEqual( res[i].im, f.sc[i], zTOL ) ) result = false;
  zAssert( zFFT + zFFTInv, result );
  zFourierFree( &f );
}

int main(int argc, char *argv[])
{
  zRandInit();
  assert_fft();
  return EXIT_SUCCESS;
}
