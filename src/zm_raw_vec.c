/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_raw_vec - raw vector and matrix : vector.
 */

#include <zm/zm_raw.h>

static void _zRawVecCats(double *v, int size, int n, va_list args);

/* touchup a raw vector. */
void zRawVecTouchup(double *v, int size)
{
  for( ; size>0; v++, size-- )
    if( zIsTiny( *v ) ) *v = 0;
}

/* create a uniform raw vector. */
void zRawVecSetAll(double *v, int size, double val)
{
  while( size-- > 0 ) *v++ = val;
}

/* create a linear space in a raw vector. */
void zRawVecLinSpace(double *v, int size, double from, double to)
{
  register int i = 0, n;
  double range;

  range = to - from;
  n = size - 1;
  while( size-- > 0 )
    *v++ = range * i++ / n + from;
}

/* create a random raw vector with a uniform range. */
void zRawVecRandUniform(double *v, int size, double min, double max)
{
  while( size-- > 0 ) *v++ = zRandF( min, max );
}

/* create a random raw vector with an arrayed range. */
__EXPORT void zRawVecRand(double *v, double *min, double *max, int size)
{
  while( size-- > 0 ) *v++ = zRandF( *min++, *max++ );
}

/* shift a raw vector by a scalar constant. */
void zRawVecShift(double *v, int size, double shift)
{
  while( size-- > 0 ) *v++ += shift;
}

/* swap components of a raw vector. */
double *zRawVecSwap(double *v, int i1, int i2)
{
  zSwap( double, v[i1], v[i2] );
  return v;
}

/* check if a raw vector is tiny. */
bool zRawVecIsTol(double *v, int size, double tol)
{
  while( size-- > 0 )
    if( !zIsTol( *v++, tol ) ) return false;
  return true;
}

/* check if a raw vector includes NaN. */
bool zRawVecIsNan(double *v, int size)
{
  for( ; size-->0; v++ )
    if( zIsNan( *v ) || zIsInf( *v ) ) return true;
  return false;
}

/* add two raw vectors. */
void zRawVecAdd(double *v1, double *v2, double *v, int size)
{
  while( size-- > 0 ) *v++ = *v1++ + *v2++;
}

/* subtract a raw vector from another. */
void zRawVecSub(double *v1, double *v2, double *v, int size)
{
  while( size-- > 0 ) *v++ = *v1++ - *v2++;
}

/* reverse a raw vector. */
void zRawVecRev(double *v1, double *v, int size)
{
  while( size-- > 0 ) *v++ = -*v1++;
}

/* multiply a raw vector by a scalar value. */
void zRawVecMul(double *v1, double k, double *v, int size)
{
  while( size-- > 0 ) *v++ = *v1++ * k;
}

/* divide a raw vector by a scalar value. */
void zRawVecDiv(double *v1, double k, double *v, int size)
{
  zRawVecMul( v1, 1.0/k, v, size );
}

/* amplify a raw vector by another. */
void zRawVecAmp(double *v1, double *amp, double *v, int size)
{
  while( size-- > 0 ) *v++ = *v1++ * *amp++;
}

/* demagnify a raw vector. */
void zRawVecDem(double *v1, double *dem, double *v, int size)
{
  while( size-- > 0 ) *v++ = *v1++ / *dem++;
}

/* concatenate a vector by another multiplied vector. */
void zRawVecCat(double *v1, double k, double *v2, double *v, int size)
{
  while( size-- > 0 ) *v++ = *v1++ + *v2++ * k;
}

/* add two raw vectors directly. */
void zRawVecAddDRC(double *v1, double *v2, int size)
{
  while( size-- > 0 ) *v1++ += *v2++;
}

/* subtract a raw vector directly from another. */
void zRawVecSubDRC(double *v1, double *v2, int size)
{
  while( size-- > 0 ) *v1++ -= *v2++;
}

/* reverse a raw vector directly. */
void zRawVecRevDRC(double *v, int size)
{
  for( ; size-->0; v++ ) *v = -*v;
}

/* multiply a raw vector directly by a scalar value. */
void zRawVecMulDRC(double *v, double k, int size)
{
  while( size-- > 0 ) *v++ *= k;
}

/* divide a raw vector directly by a scalar value. */
void zRawVecDivDRC(double *v, double k, int size)
{
  zRawVecMulDRC( v, 1.0/k, size );
}

/* amplify a raw vector directly by another. */
void zRawVecAmpDRC(double *v, double *amp, int size)
{
  while( size-- > 0 ) *v++ *= *amp++;
}

/* demagnify a raw vector directly by another. */
void zRawVecDemDRC(double *v, double *dem, int size)
{
  while( size-- > 0 ) *v++ /= *dem++;
}

/* concatenate a raw vector directly by another vector
 * multiplied by a scalar. */
void zRawVecCatDRC(double *v1, double k, double *v2, int size)
{
  while( size-- > 0 ) *v1++ += *v2++ * k;
}

/* (static)
 * concatenate multiple sets of a scalar and a raw vector,
 * which is internally called to manage variable arguments. */
void _zRawVecCats(double *v, int size, int n, va_list args)
{
  register int i;
  double k, *vec;

  for( i=0; i<n; i++ ){
    k = (double)va_arg( args, double );
    vec = (double *)va_arg( args, double* );
    zRawVecCatDRC( v, k, vec, size );
  }
}

/* concatenate multiple sets of a scalar and a raw vector. */
void zRawVecCats(double *v, int size, int n, ...)
{
  va_list args;

  va_start( args, n );
  _zRawVecCats( v, size, n, args );
  va_end( args );
}

/* linear sum of the multiple sets of a scalar and a raw vector. */
void zRawVecLS(double *v, int size, int n, ...)
{
  va_list args;

  zRawVecZero( v, size );
  va_start( args, n );
  _zRawVecCats( v, size, n, args );
  va_end( args );
}

/* interior division of two raw vectors. */
void zRawVecInterDiv(double *v1, double *v2, double ratio, double *v, int size)
{
  while( size-- > 0 ){
    *v++ = *v1 + ratio * ( *v2 - *v1 );
    v1++;
    v2++;
  }
}

/* scale a raw vector with two boundary vectors. */
void zRawVecScale(double *x, double *min, double *max, double *v, int size)
{
  while( size-- > 0 ){
    *v++ = *min + *x++ * ( *max - *min );
    min++;
    max++;
  }
}

/* inner product of two raw vectors. */
double zRawVecInnerProd(double *v1, double *v2, int size)
{
  double s=0, s_prev=0, v, q=0, r;

  while( size-- > 0 ){
    s = s_prev + ( v = *v1++ * *v2++ );
    r = s - s_prev;
    q += v - r;
    s_prev = s;
  }
  return s + q;
}

/* squared norm of a raw vector. */
double zRawVecSqrNorm(double *v, int size)
{
  return zRawVecInnerProd( v, v, size );
}

/* weighted squared norm of a raw vector. */
double zRawVecWSqrNorm(double *v, double *w, int size)
{
  double s=0, s_prev=0, vs, q=0, r;

  for( ; size-->0; v++, w++ ){
    s = s_prev + ( vs = *w * *v * *v );
    r = s - s_prev;
    q += vs - r;
    s_prev = s;
  }
  return s + q;
}

/* normalize a raw vector. */
double *zRawVecNormalize(double *src, int size, double *dest)
{
  double l;

  l = zRawVecNorm( src, size );
  if( zIsTiny( l ) ){
    ZRUNERROR( ZM_ERR_ZERODIV );
    return NULL;
  }
  zRawVecMul( src, 1.0/l, dest, size );
  return dest;
}

/* squared distance between two raw vectors. */
double zRawVecSqrDist(double *v1, double *v2, int size)
{
  register int i;
  double d;

  for( d=0, i=0; i<size; i++ )
    d += zSqr( v1[i] - v2[i] );
  return d;
}

/* print a raw vector out to a file. */
void zRawVecFPrint(FILE *fp, double *v, int size)
{
  while( size-- > 0 )
    fprintf( fp, "%.10g ", *v++ );
  fprintf( fp, "\n" );
}
