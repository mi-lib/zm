/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_raw_vec - raw vector and matrix : vector.
 */

#include <zm/zm_raw.h>

static void _zRawVecCats(double *v, int size, int n, va_list args);

/* zRawVecTouchup
 * - touchup a raw vector.
 */
void zRawVecTouchup(double *v, int size)
{
  for( ; size>0; v++, size-- )
    if( zIsTiny( *v ) ) *v = 0;
}

/* zRawVecSetAll
 * - create a uniform raw vector.
 */
void zRawVecSetAll(double *v, int size, double val)
{
  while( size-- > 0 ) *v++ = val;
}

/* zRawVecLinSpace
 * - create a linear space in a raw vector.
 */
void zRawVecLinSpace(double *v, int size, double from, double to)
{
  register int i = 0, n;
  double range;

  range = to - from;
  n = size - 1;
  while( size-- > 0 )
    *v++ = range * i++ / n + from;
}

/* zRawVecRandUniform
 * - create a random raw vector with a uniform range.
 */
void zRawVecRandUniform(double *v, int size, double min, double max)
{
  while( size-- > 0 ) *v++ = zRandF( min, max );
}

/* zRawVecRand
 * - create a random raw vector with an arrayed range.
 */
__EXPORT void zRawVecRand(double *v, double *min, double *max, int size)
{
  while( size-- > 0 ) *v++ = zRandF( *min++, *max++ );
}

/* zRawVecShift
 * - shift a raw vector by a scalar constant.
 */
void zRawVecShift(double *v, int size, double shift)
{
  while( size-- > 0 ) *v++ += shift;
}

/* zRawVecSwap
 * - swap components of a raw vector.
 */
double *zRawVecSwap(double *v, int i1, int i2)
{
  zSwap( double, v[i1], v[i2] );
  return v;
}

/* zRawVecIsTol
 * - test if vector is tiny.
 */
bool zRawVecIsTol(double *v, int size, double tol)
{
  while( size-- > 0 )
    if( !zIsTol( *v++, tol ) ) return false;
  return true;
}

/* zRawVecIsNan
 * - test if vector includes NaN.
 */
bool zRawVecIsNan(double *v, int size)
{
  for( ; size-->0; v++ )
    if( zIsNan( *v ) || zIsInf( *v ) ) return true;
  return false;
}

/* zRawVecAdd
 * - add two raw vectors.
 */
void zRawVecAdd(double *v1, double *v2, double *v, int size)
{
  while( size-- > 0 ) *v++ = *v1++ + *v2++;
}

/* zRawVecSub
 * - subtract raw vectors.
 */
void zRawVecSub(double *v1, double *v2, double *v, int size)
{
  while( size-- > 0 ) *v++ = *v1++ - *v2++;
}

/* zRawVecRev
 * - reverse a raw vector.
 */
void zRawVecRev(double *v1, double *v, int size)
{
  while( size-- > 0 ) *v++ = -*v1++;
}

/* zRawVecMul
 * - multiply a raw vector by a scalar value.
 */
void zRawVecMul(double *v1, double k, double *v, int size)
{
  while( size-- > 0 ) *v++ = *v1++ * k;
}

/* zRawVecDiv
 * - divide a raw vector by a scalar value.
 */
void zRawVecDiv(double *v1, double k, double *v, int size)
{
  zRawVecMul( v1, 1.0/k, v, size );
}

/* zRawVecAmp
 * - amplify a raw vector.
 */
void zRawVecAmp(double *v1, double *amp, double *v, int size)
{
  while( size-- > 0 ) *v++ = *v1++ * *amp++;
}

/* zRawVecDem
 * - demagnify a raw vector.
 */
void zRawVecDem(double *v1, double *dem, double *v, int size)
{
  while( size-- > 0 ) *v++ = *v1++ / *dem++;
}

/* zRawVecCat
 * - concatenate a vector by another multiplied vector.
 */
void zRawVecCat(double *v1, double k, double *v2, double *v, int size)
{
  while( size-- > 0 ) *v++ = *v1++ + *v2++ * k;
}

/* zRawVecAddDRC
 * - directly add two raw vectors.
 */
void zRawVecAddDRC(double *v1, double *v2, int size)
{
  while( size-- > 0 ) *v1++ += *v2++;
}

/* zRawVecSubDRC
 * - directly subtract raw vectors.
 */
void zRawVecSubDRC(double *v1, double *v2, int size)
{
  while( size-- > 0 ) *v1++ -= *v2++;
}

/* zRawVecRevDRC
 * - directly reverse a raw vector.
 */
void zRawVecRevDRC(double *v, int size)
{
  for( ; size-->0; v++ ) *v = -*v;
}

/* zRawVecMulDRC
 * - directly multiply a raw vector by a scalar value.
 */
void zRawVecMulDRC(double *v, double k, int size)
{
  while( size-- > 0 ) *v++ *= k;
}

/* zRawVecDivDRC
 * - directly divide a raw vector by a scalar value.
 */
void zRawVecDivDRC(double *v, double k, int size)
{
  zRawVecMulDRC( v, 1.0/k, size );
}

/* zRawVecAmpDRC
 * - directly amplify a raw vector.
 */
void zRawVecAmpDRC(double *v, double *amp, int size)
{
  while( size-- > 0 ) *v++ *= *amp++;
}

/* zRawVecDemDRC
 * - directly demagnify a raw vector.
 */
void zRawVecDemDRC(double *v, double *dem, int size)
{
  while( size-- > 0 ) *v++ /= *dem++;
}

/* zRawVecCatDRC
 * - directly concatenates a vector by another vector
 *   multiplied by a scalar.
 */
void zRawVecCatDRC(double *v1, double k, double *v2, int size)
{
  while( size-- > 0 ) *v1++ += *v2++ * k;
}

/* (static)
 * _zRawVecCats
 * - concatenate multiple sets of a scalar and a vector.
 *   (internally called to manage variable arguments.)
 */
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

/* zRawVecCats
 * - concatenate multiple sets of a scalar and a vector.
 */
void zRawVecCats(double *v, int size, int n, ...)
{
  va_list args;

  va_start( args, n );
  _zRawVecCats( v, size, n, args );
  va_end( args );
}

/* zRawVecLS
 * - linear sum of multiple sets of a scalar and a vector.
 */
void zRawVecLS(double *v, int size, int n, ...)
{
  va_list args;

  zRawVecClear( v, size );
  va_start( args, n );
  _zRawVecCats( v, size, n, args );
  va_end( args );
}

/* zRawVecInnerProd
 * - inner product of two raw vectors.
 */
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

/* zRawVecSqrNorm
 * - squared norm of a raw vector.
 */
double zRawVecSqrNorm(double *v, int size)
{
  return zRawVecInnerProd( v, v, size );
}

/* zRawVecWSqrNorm
 * - weighted squared norm of a raw vector.
 */
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

/* zRawVecNormalize
 * - normalization of a raw vector.
 */
double *zRawVecNormalize(double *src, int size, double *dest)
{
  double l;

  l = zRawVecNorm( src, size );
  if( zIsTiny( l ) ){
    ZRUNERROR( ZM_ERR_ZERODIV );
    return NULL;
  }
  zRawVecDiv( src, l, dest, size );
  return dest;
}

/* zRawVecSqrDist
 * - squared distance between two raw vectors.
 */
double zRawVecSqrDist(double *v1, double *v2, int size)
{
  register int i;
  double d, e;

  for( d=0, i=0; i<size; i++ ){
    e = v1[i] - v2[i];
    d += zSqr(e);
  }
  return d;
}

/* zRawVecFWrite
 * - output raw vector.
 */
void zRawVecFWrite(FILE *fp, double *v, int size)
{
  while( size-- > 0 )
    fprintf( fp, "%.10g ", *v++ );
  fprintf( fp, "\n" );
}
