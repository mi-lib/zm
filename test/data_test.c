#include <zm/zm_data.h>
#include <zm/zm_fourier.h>

void assert_sort(void)
{
  const int n = 1000;
  int i;
  double *data;
  bool result = true;

  data = zAlloc( double, n );
  for( i=0; i<n; i++ )
    data[i] = zRandF(-10,10);
  zDataSort( data, n );
  for( i=1; i<n; i++ )
    if( data[i-1] > data[i] ) result = false;
  free( data );
  zAssert( zDataSort, result );
}

void assert_sort_index(void)
{
  const int n = 1000;
  int i;
  double *data_org, *data, *data_sorted;
  zIndex index;
  bool result1 = true, result2;

  data_org = zAlloc( double, n );
  data = zAlloc( double, n );
  data_sorted = zAlloc( double, n );
  index = zIndexCreate( n );
  for( i=0; i<n; i++ )
    data_org[i] = data[i] = data_sorted[i] = zRandF(-10,10);
  zDataSort( data_sorted, n );
  zDataSortIndex( data, n, index );
  for( i=0; i<n; i++ ){
    if( data[zIndexElem(index,i)] != data_sorted[i] ) result1 = false;
  }
  result2 = memcmp( data, data_org, n ) == 0;
  free( data_org );
  free( data );
  free( data_sorted );
  zIndexFree( index );
  zAssert( zDataSortIndex, result1 && result2 );
}

void assert_data_select(void)
{
  const int n = 1000;
  int i;
  double *data_org, *data, *data_sorted, val;
  bool result1 = true, result2;

  data_org = zAlloc( double, n );
  data = zAlloc( double, n );
  data_sorted = zAlloc( double, n );
  for( i=0; i<n; i++ )
    data_org[i] = data[i] = data_sorted[i] = zRandF(-100,100);
  zDataSort( data_sorted, n );
  for( i=0; i<n; i++ ){
    val = zDataSelect( data, n, i );
    if( val != data_sorted[i] ) result1 = false;
  }
  result2 = memcmp( data, data_org, n ) == 0;
  free( data_org );
  free( data );
  free( data_sorted );
  zAssert( zDataSelect, result1 && result2 );
}

void assert_data_median(void)
{
  int i, n;
  double *data;
  bool result;

  data = zAlloc( double, ( n = 11 ) );
  for( i=0; i<n; i++ )
    data[i] = i;
  result = zDataMedian( data, n ) == 5;
  free( data );
  zAssert( zDataMedian (odd members case), result );

  data = zAlloc( double, ( n = 10 ) );
  for( i=0; i<n; i++ )
    data[i] = i;
  result = zDataMedian( data, n ) == 4.5;
  free( data );
  zAssert( zDataMedian (even members case), result );
}

#define FFT_TEST_N   100
#define FFT_TEST_DIM  10

void assert_fft(void)
{
  zFourier f;
  zVec data;
  zCVec res;
  zVec inv;
  int i;
  bool result1 = true;
  bool result2 = true;

  zFourierAlloc( &f, FFT_TEST_DIM );
  for( i=0; i<f.n; i++ ){
    zFourierSetSinCoeff( &f, i, zRandF(-10,10) );
    zFourierSetCosCoeff( &f, i, zRandF(-10,10) );
  }
  data = zVecAlloc( FFT_TEST_N );
  res = zCVecAlloc( FFT_TEST_N );
  inv = zVecAlloc( FFT_TEST_N );
  for( i=0; i<FFT_TEST_N; i++ )
    zVecSetElemNC( data, i, zFourierVal( &f, 2*zPI*i/FFT_TEST_N ) );

  zFFT( data, res );
  for( i=0; i<FFT_TEST_DIM; i++ )
    if( !zEqual( zCVecElemNC(res,i)->re, f.cc[i], zTOL ) ||
        !zEqual( zCVecElemNC(res,i)->im, f.sc[i], zTOL ) ) result1 = false;
  zFFTInv( res, inv );
  if( !zVecEqual( data, inv, zTOL ) ) result2 = false;
  zVecFree( data );
  zCVecFree( res );
  zVecFree( inv );
  zFourierFree( &f );
  zAssert( zFFT, result1 );
  zAssert( zFFTInv, result2 );
}

void assert_fft2(void)
{
  zFourier f1, f2;
  zMat data;
  zCMat res;
  zMat inv;
  zComplex c1, c2, c;
  int i, j;
  bool result1 = true;
  bool result2 = true;

  zFourierAlloc( &f1, FFT_TEST_DIM );
  zFourierAlloc( &f2, FFT_TEST_DIM );
  for( i=0; i<f1.n; i++ ){
    zFourierSetSinCoeff( &f1, i, zRandF(-10,10) );
    zFourierSetCosCoeff( &f1, i, zRandF(-10,10) );
    zFourierSetSinCoeff( &f2, i, zRandF(-10,10) );
    zFourierSetCosCoeff( &f2, i, zRandF(-10,10) );
  }
  data = zMatAlloc( FFT_TEST_N, FFT_TEST_N );
  res = zCMatAlloc( FFT_TEST_N, FFT_TEST_N );
  inv = zMatAlloc( FFT_TEST_N, FFT_TEST_N );
  for( i=0; i<FFT_TEST_N; i++ )
    for( j=0; j<FFT_TEST_N; j++ )
      zMatSetElemNC( data, i, j, zFourierVal( &f1, 2*zPI*j/FFT_TEST_N ) * zFourierVal( &f2, 2*zPI*i/FFT_TEST_N ) );

  zFFT2( data, res );
  for( i=0; i<FFT_TEST_DIM; i++ ){
    zComplexCreate( &c2, f2.cc[i], f2.sc[i] );
    for( j=0; j<FFT_TEST_DIM; j++ ){
      zComplexCreate( &c1, f1.cc[j], f1.sc[j] );
      zComplexCMul( &c1, &c2, &c );
      if( !zComplexEqual( zCMatElemNC(res,i,j), &c, zTOL ) ){
        eprintf( "(%d,%d)\n", i, j );
        zComplexFPrint( stderr, zCMatElemNC(res,i,j) ); zFEndl( stderr );
        zComplexFPrint( stderr, &c ); zFEndl( stderr );
        zFEndl( stderr );
        result1 = false;
      }
    }
  }
  zFFT2Inv( res, inv );
  if( !zMatEqual( data, inv, zTOL ) ) result2 = false;

  zMatFree( data );
  zCMatFree( res );
  zMatFree( inv );
  zFourierFree( &f1 );
  zFourierFree( &f2 );
  zAssert( zFFT2, result1 );
  zAssert( zFFT2Inv, result2 );
}

int main(int argc, char *argv[])
{
  zRandInit();
  assert_sort();
  assert_sort_index();
  assert_data_select();
  assert_data_median();

  assert_fft();
  assert_fft2();
  return EXIT_SUCCESS;
}
