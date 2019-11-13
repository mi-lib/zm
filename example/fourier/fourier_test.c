#include <zm/zm_fourier.h>

#define FFT_TEST_N 100
#define FFT_TEST_DIM 9

int main(int argc, char *argv[])
{
  zFourier f;
  double phase;
  register int i;

  zRandInit();
  zFourierAlloc( &f, FFT_TEST_DIM );
  for( i=0; i<f.n; i++ ){
    zFourierSetSinCoeff( &f, i, zRandF(-10,10) );
    zFourierSetCosCoeff( &f, i, zRandF(-10,10) );
  }
  for( i=0; i<FFT_TEST_N; i++ ){
    phase = 2*zPI*i/FFT_TEST_N;
    printf( "%.10f %.10f %.10f\n", zFourierVal(&f,phase), zFourierVel(&f,phase), zFourierAcc(&f,phase) );
  }
  zFourierFree( &f );
  return 0;
}
