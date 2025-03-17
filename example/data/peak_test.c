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
  zIndex peakidx, peakidx_sg;
  int i;
  FILE *fp;

  zRandInit();
  zFourierAlloc( &fourier, FOURIER_DIM );
  for( i=0; i<fourier.n; i++ ){
    zFourierSetSinCoeff( &fourier, i, zRandF(-10,10) );
    zFourierSetCosCoeff( &fourier, i, zRandF(-10,10) );
  }
  for( i=-DIV; i<=DIV; i++ ){
    sample[i+DIV] = zFourierVal(&fourier,XMAX*2*zPI*i/DIV) + zRandF(-2,2);
    smooth[i+DIV] = 0;
  }

  zDataSmoothSG( sample, SAMPLE_N, 50, 5, smooth );
  peakidx = zDataPeak( sample, SAMPLE_N, 100 );
  peakidx_sg = zDataPeakSG( sample, SAMPLE_N, 100, 5 );
  zDataSortIndex( sample, SAMPLE_N, peakidx );
  zDataSortIndex( smooth, SAMPLE_N, peakidx_sg );
  fp = fopen( "sample.dat", "w" );
  for( i=0; i<SAMPLE_N; i++ )
    fprintf( fp, "%.10f %.10f\n", sample[i], smooth[i] );
  fclose( fp );
  fp = fopen( "peak.dat", "w" );
  for( i=0; i<zArraySize(peakidx); i++ )
    fprintf( fp, "%d %.10f\n", zIndexElemNC(peakidx,i), sample[zIndexElemNC(peakidx,i)] );
  fclose( fp );
  fp = fopen( "peak_sg.dat", "w" );
  for( i=0; i<zArraySize(peakidx_sg); i++ )
    fprintf( fp, "%d %.10f\n", zIndexElemNC(peakidx_sg,i), smooth[zIndexElemNC(peakidx_sg,i)] );
  fclose( fp );
  zIndexFree( peakidx );
  zIndexFree( peakidx_sg );
  return 0;
}
