/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_fft - Fast Fourier Transform.
 */

#include <zm/zm_fft.h>

static void _zFFTInner(zComplex data[], zComplex buf[], int n, double theta);

/* (static)
 * _zFFTInner
 * - internal sub-FFT procedure.
 */
void _zFFTInner(zComplex data[], zComplex buf[], int n, double theta)
{
  register int i, j, k;
  int radix, n_radix;
  zComplex x, w, tmp;

  if( n <= 1 ) return;
  /* factorization */
  for( radix=2; radix*radix<=n; radix++ )
    if( n % radix == 0 ) break;
  if( n % radix != 0 ) radix = n;
  n_radix = n / radix;
  /* butterflies */
  for( i=0; i<n_radix; i++ )
    for( j=0; j<radix; j++ ){
      x = data[i];
      for( k=n_radix; k<n; k+=n_radix ){
        zComplexCreate( &w, cos(theta*j*k), sin(theta*j*k) );
	zComplexCMul( &w, &data[k+i], &tmp );
        zComplexAdd( &x, &tmp, &x );
      }
      zComplexCreate( &w, cos(theta*j*i), sin(theta*j*i) );
      zComplexCMul( &x, &w, &buf[n_radix*j+i] );
    }
  /* recursive sub-FFT */
  for( i=0; i<n; i+=n_radix )
    _zFFTInner( &buf[i], data, n_radix, theta*radix );
  /* rearrangement */
  for( i=0; i<n_radix; i++ )
    for( j=0; j<radix; j++ )
      data[radix*i+j] = buf[n_radix*j+i];
}

/* zFFT
 * - Fast Fourier Transform
 */
bool zFFT(double data[], int n, zComplex res[])
{
  register int i;
  zComplex *buf;

  if( !( buf = zAlloc( zComplex, n ) ) ){
    ZALLOCERROR();
    return false;
  }
  for( i=0; i<n; i++ )
    zComplexCreate( &res[i], data[i], 0 );
  /* internal FFT */
  _zFFTInner( res, buf, n, 2*zPI/n );
  /* amplitude correction */
  zComplexMul( &res[0], 1.0/n, &res[0] );
  for( i=1; i<n; i++ )
    zComplexMul( &res[i], 2.0/n, &res[i] );
  zFree( buf );
  return true;
}

/* zFFTInv
 * - Inverse Fast Fourier Transform
 */
bool zFFTInv(zComplex data[], int n, zComplex res[])
{
  register int i;
  zComplex *buf;

  if( !( buf = zAlloc( zComplex, n ) ) ){
    ZALLOCERROR();
    return false;
  }
  /* amplitude correction */
  res[0] = data[0];
  for( i=1; i<n; i++ )
    zComplexMul( &data[i], 0.5, &res[i] );
  /* internal inverse FFT */
  _zFFTInner( res, buf, n,-2*zPI/n );
  zFree( buf );
  return true;
}
