#include <zm/zm_data.h>
#include <zm/zm_fourier.h>

#define XMAX 1.0
#define DIV  1000

#define SAMPLE_N (DIV*2+1)

#define FOURIER_DIM 5
#define DIM 5

int main(int argc, char *argv[])
{
  zFourier fourier;
  double sample[SAMPLE_N];
  double smooth[SAMPLE_N];
  int i;

  zRandInit();
  zFourierAlloc( &fourier, FOURIER_DIM );
  for( i=0; i<fourier.n; i++ ){
    zFourierSetSinCoeff( &fourier, i, zRandF(-10,10) );
    zFourierSetCosCoeff( &fourier, i, zRandF(-10,10) );
  }
  for( i=-DIV; i<=DIV; i++ ){
    sample[i+DIV] = zFourierVal(&fourier,XMAX*2*zPI*i/DIV) + zRandF(-3,3);
    smooth[i+DIV] = 0;
  }

  zDataSmoothSG( sample, SAMPLE_N, 15, 3, smooth );
  for( i=0; i<SAMPLE_N; i++ ){
    printf( "%.10f %.10f\n", sample[i], smooth[i] );
  }
  return 0;
}
