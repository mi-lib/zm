/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_data_fft - data analysis: fast Fourier transformation.
 */

#include <zm/zm_data.h>

/* factorization of the number of samples for FFT. */
static uint _zFFTRadix(size_t n)
{
  uint radix;

  for( radix=2; radix*radix<=n; radix++ )
    if( n % radix == 0 ) break;
  if( n % radix != 0 ) radix = n;
  return radix;
}

/* internal sub-FFT procedure. */
static void _zFFTInner(zComplex data[], zComplex buf[], size_t n, double theta)
{
  uint i, j, k;
  uint radix, n_radix;
  double theta_j, theta_jk, theta_ji;
  zComplex x, w, tmp;

  if( n <= 1 ) return;
  n_radix = n / ( radix = _zFFTRadix( n ) );
  /* butterflies */
  for( i=0; i<n_radix; i++ )
    for( j=0; j<radix; j++ ){
      theta_j = theta * j;
      zComplexCopy( &data[i], &x );
      for( k=n_radix; k<n; k+=n_radix ){
        theta_jk = theta_j * k;
        zComplexCreate( &w, cos(theta_jk), sin(theta_jk) );
	zComplexCMul( &w, &data[k+i], &tmp );
        zComplexAdd( &x, &tmp, &x );
      }
      theta_ji = theta_j * i;
      zComplexCreate( &w, cos(theta_ji), sin(theta_ji) );
      zComplexCMul( &x, &w, &buf[n_radix*j+i] );
    }
  /* recursive sub-FFT */
  for( i=0; i<n; i+=n_radix )
    _zFFTInner( &buf[i], data, n_radix, theta*radix );
  /* rearrangement */
  for( i=0; i<n_radix; i++ )
    for( j=0; j<radix; j++ )
      zComplexCopy( &buf[n_radix*j+i], &data[radix*i+j] );
}

/* Fast Fourier Transformation. */
bool zFFT(zVec data, zCVec result)
{
  zComplex *buf;

  if( zVecSize(data) != zCVecSize(result) ){
    ZRUNERROR( ZM_ERR_FFT_SIZEMISMATCH_VEC, zVecSize(data), zCVecSize(result) );
    return false;
  }
  if( !( buf = zAlloc( zComplex, zVecSizeNC(data) ) ) ){
    ZALLOCERROR();
    return false;
  }
  zVecToCVec( data, result );
  /* internal FFT */
  _zFFTInner( zCVecBuf(result), buf, zVecSizeNC(data), zPIx2/zVecSizeNC(data) );
  /* amplitude correction */
  zComplexMulDRC( zCVecElemNC(result,0), 0.5 );
  zCVecMulDRC( result, 2.0/zCVecSizeNC(result) );
  zFree( buf );
  return true;
}

/* inverse fast Fourier transformation. */
bool zFFTInv(zCVec data, zVec result)
{
  zComplex *buf;
  zCVec cres;
  bool ret = true;

  if( zCVecSize(data) != zVecSize(result) ){
    ZRUNERROR( ZM_ERR_FFT_SIZEMISMATCH_VEC, zCVecSize(data), zVecSize(result) );
    return false;
  }
  buf = zAlloc( zComplex, zCVecSizeNC(data) );
  cres = zCVecAlloc( zCVecSizeNC(data) );
  if( !buf || !cres ){
    ZALLOCERROR();
    ret = false;
    goto TERMINATE;
  }
  /* amplitude correction */
  zCVecMulNC( data, 0.5, cres );
  zComplexMulDRC( zCVecElemNC(cres,0), 2 );
  /* internal inverse FFT */
  _zFFTInner( zCVecBuf(cres), buf, zCVecSizeNC(data), -zPIx2/zCVecSizeNC(data) );
  zCVecToReVec( cres, result );

 TERMINATE:
  zFree( buf );
  zCVecFree( cres );
  return ret;
}

/* 2-dimensional Fast Fourier Transformation. */
bool zFFT2(zMat data, zCMat result)
{
  zComplex *buf;
  zCVec tmp;
  bool ret = false;
  int i;

  if( zMatRowSize(data) != zCMatRowSize(result) || zMatColSize(data) != zCMatColSize(result) ){
    ZRUNERROR( ZM_ERR_FFT_SIZEMISMATCH_MAT, zMatRowSize(data), zMatColSize(data), zCMatRowSize(result), zCMatColSize(result) );
    return false;
  }
  buf = zAlloc( zComplex, zMax( zMatColSizeNC(data), zMatRowSizeNC(data) ) );
  tmp = zCVecAlloc( zMatRowSizeNC(data) );
  if( !buf || !tmp ){
    ZALLOCERROR();
    goto TERMINATE;
  }
  zMatToCMat( data, result );
  /* rowwise FFT */
  for( i=0; i<zMatRowSizeNC(data); i++ ){
    _zFFTInner( zCMatRowBuf(result,i), buf, zMatColSizeNC(data), zPIx2/zMatColSizeNC(data) );
    zComplexMulDRC( zCMatElemNC(result,i,0), 0.5 ); /* amplitude correction */
  }
  zCMatMulDRC( result, 2.0/zMatColSizeNC(data) ); /* amplitude correction */
  /* columnwise FFT */
  for( i=0; i<zMatColSizeNC(data); i++ ){
    zCMatGetColNC( result, i, tmp );
    _zFFTInner( zCVecBuf(tmp), buf, zVecSizeNC(tmp), zPIx2/zVecSizeNC(tmp) );
    zComplexMulDRC( zCVecElemNC(tmp,0), 0.5 ); /* amplitude correction */
    zCMatPutColNC( result, i, tmp );
  }
  zCMatMulDRC( result, 2.0/zMatRowSizeNC(data) ); /* amplitude correction */
  ret = true;

 TERMINATE:
  zCVecFree( tmp );
  free( buf );
  return ret;
}

/* 2-dimensional inverse fast Fourier transformation. */
bool zFFT2Inv(zCMat data, zMat result)
{
  zComplex *buf;
  zCMat cres;
  zCVec tmp;
  bool ret = false;
  int i, j;

  if( zCMatRowSize(data) != zMatRowSize(result) || zCMatColSize(data) != zMatColSize(result) ){
    ZRUNERROR( ZM_ERR_FFT_SIZEMISMATCH_MAT, zCMatRowSize(data), zCMatColSize(data), zMatRowSize(result), zMatColSize(result) );
    return false;
  }
  buf = zAlloc( zComplex, zMax( zMatColSizeNC(data), zMatRowSizeNC(data) ) );
  cres = zCMatAlloc( zMatRowSizeNC(data), zMatColSizeNC(data) );
  tmp = zCVecAlloc( zMax( zMatRowSizeNC(data), zMatColSizeNC(data) ) );
  if( !buf || !cres || !tmp ){
    ZALLOCERROR();
    goto TERMINATE;
  }
  /* columnwise inverse FFT */
  zCVecSetSize( tmp, zMatRowSizeNC(data) );
  for( j=0; j<zMatColSizeNC(data); j++ ){
    zCMatGetColNC( data, j, tmp );
    zCVecMulNCDRC( tmp, 0.5 ); /* amplitude correction */
    zComplexMulDRC( zCVecElemNC(tmp,0), 2 );
    _zFFTInner( zCVecBuf(tmp), buf, zVecSizeNC(tmp), -zPIx2/zVecSizeNC(tmp) );
    zCMatPutColNC( cres, j, tmp );
  }
  /* rowwise FFT */
  zCVecSetSize( tmp, zMatColSizeNC(data) );
  for( i=0; i<zMatRowSizeNC(result); i++ ){
    zCMatGetRowNC( cres, i, tmp );
    zCVecMulNCDRC( tmp, 0.5 ); /* amplitude correction */
    zComplexMulDRC( zCVecElemNC(tmp,0), 2 );
    _zFFTInner( zCVecBuf(tmp), buf, zVecSizeNC(tmp), -zPIx2/zVecSizeNC(tmp) );
    for( j=0; j<zMatColSizeNC(result); j++ )
      zMatSetElemNC( result, i, j, zCVecElemNC(tmp,j)->re );
  }
  ret = true;

 TERMINATE:
  zCVecFree( tmp );
  zCMatFree( cres );
  free( buf );
  return ret;
}
