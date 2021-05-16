/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_vec_vec - vector and matrix class : vector class.
 */

#include <zm/zm_vec.h>

/* ********************************************************** */
/* CLASS: zVec
 * double precision floating point value vector class
 * ********************************************************** */

/* set vector components from variable argument list. */
zVec zVecSetElemVList(zVec v, va_list args)
{
  register int i;

  for( i=0; i<zVecSizeNC(v); i++ )
    zVecSetElemNC( v, i, (double)va_arg( args, double ) );
  return v;
}

/* set vector components from argument list. */
zVec zVecSetElemList(zVec v, ... )
{
  va_list args;

  va_start( args, v );
  zVecSetElemVList( v, args );
  va_end( args );
  return v;
}

/* allocate memory for a vector. */
zVec zVecAlloc(int size)
{
  zVec v;

  if( !( v = zAlloc( zVecStruct, 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  if( !( zVecBufNC(v) = zAlloc( double, size ) ) ){
    ZALLOCERROR();
    free( v );
    return NULL;
  }
  zVecSetSize( v, size );
  return v;
}

/* create a vector from argument list. */
zVec zVecCreateList(int size, ... )
{
  zVec v;
  va_list args;

  if( !( v = zVecAlloc( size ) ) ) return NULL;
  va_start( args, size );
  zVecSetElemVList( v, args );
  va_end( args );
  return v;
}

/* free memory for a vector. */
void zVecFree(zVec v)
{
  if( !v ) return;
  zFree( zVecBufNC( v ) );
  free( v );
}

/* free memory for multiple vectors at once. */
void zVecFreeAO(int n, ...)
{
  va_list arg;
  zVec v;
  register int i;

  va_start( arg, n );
  for( i=0; i<n; i++ ){
    v = va_arg( arg, zVec );
    zVecFree( v );
  }
  va_end( arg );
}

/* zero a vector. */
zVec zVecZero(zVec v)
{
  zRawVecZero( zVecBufNC(v), zVecSizeNC(v) );
  return v;
}

/* touchup a vector. */
zVec zVecTouchup(zVec v)
{
  zRawVecTouchup( zVecBufNC(v), zVecSizeNC(v) );
  return v;
}

/* copy a vector without checking size consistency. */
zVec zVecCopyNC(zVec src, zVec dest)
{
  zRawVecCopy( zVecBufNC(src), zVecBufNC(dest), zVecSizeNC(dest) );
  return dest;
}

/* copy a vector. */
zVec zVecCopy(zVec src, zVec dest)
{
  return zVecSizeIsEqual( src, dest ) ?
    zVecCopyNC( src, dest ) : NULL;
}

/* copy a vector from an array of double-precision floating-point values. */
zVec zVecCopyArray(double array[], int s, zVec v)
{
  if( zVecSizeNC(v) != s ) return NULL;
  zRawVecCopy( array, zVecBufNC(v), s );
  return v;
}

/* clone a vector. */
zVec zVecClone(zVec src)
{
  zVec dest;

  if( ( dest = zVecAlloc( zVecSizeNC(src) ) ) )
    zVecCopyNC( src, dest );
  return dest;
}

/* create a vector from an array of double-precision floating-point values. */
zVec zVecCloneArray(double array[], int s)
{
  zVec v;

  if( ( v = zVecAlloc( s ) ) )
    zVecCopyArray( array, s, v );
  return v;
}

/* get a partial vector without checking the size validity. */
zVec zVecGetNC(zVec src, int pos, zVec dest)
{
  zRawVecGet( zVecBufNC(src), pos, zVecBufNC(dest), zVecSizeNC(dest) );
  return dest;
}

/* get a parttial vector. */
zVec zVecGet(zVec src, int pos, zVec dest)
{
  if( pos + zVecSizeNC(dest) > zVecSizeNC(src) ){
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );
    return NULL;
  }
  return zVecGetNC( src, pos, dest );
}

/* put a partial vector without checking the size validity. */
zVec zVecPutNC(zVec dest, int pos, zVec src)
{
  zRawVecPut( zVecBufNC(dest), pos, zVecBufNC(src), zVecSizeNC(src) );
  return dest;
}

/* put a partial vector. */
zVec zVecPut(zVec dest, int pos, zVec src)
{
  if( pos + zVecSizeNC(src) > zVecSizeNC(dest) ){
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );
    return NULL;
  }
  return zVecPutNC( dest, pos, src );
}

/* create a uniform vector. */
zVec zVecSetAll(zVec v, double val)
{
  zRawVecSetAll( zVecBufNC(v), zVecSizeNC(v), val );
  return v;
}

/* create a linear space vector. */
zVec zVecLinSpace(zVec v, double from, double to)
{
  zRawVecLinSpace( zVecBufNC(v), zVecSizeNC(v), from, to );
  return v;
}

/* create a random vector with a uniform range. */
zVec zVecRandUniform(zVec v, double min, double max)
{
  zRawVecRandUniform( zVecBufNC(v), zVecSizeNC(v), min, max );
  return v;
}

/* create a random vector with range vectors. */
zVec zVecRand(zVec v, zVec min, zVec max)
{
  zRawVecRand( zVecBufNC(v), zVecBufNC(min), zVecBufNC(max), zVecSizeNC(v) );
  return v;
}

/* shift vector by a constant scalar value. */
zVec zVecShift(zVec v, double shift)
{
  zRawVecShift( zVecBufNC(v), zVecSizeNC(v), shift );
  return v;
}

/* swap vector components without checking size. */
zVec zVecSwapNC(zVec v, int i1, int i2)
{
  zRawVecSwap( zVecBufNC(v), i1, i2 );
  return v;
}

/* swap vector components. */
zVec zVecSwap(zVec v, int i1, int i2)
{
  if( i1 > zVecSizeNC(v) || i2 > zVecSizeNC(v) ){
    ZRUNWARN( ZM_ERR_INV_INDEX );
    return NULL;
  }
  return zVecSwapNC( v, i1, i2 );
}

/* comparison function for zVecSort. */
static int _zVecSortCmp(void *p1, void *p2, void *priv)
{
  double d;

  d = ((double *)priv)[*(int*)p1] - ((double *)priv)[*(int*)p2];
  if( d > 0 ) return 1;
  if( d < 0 ) return -1;
  return 0;
}

/* rearrange index so as to sort vector. */
void zVecSort(zVec v, zIndex idx)
{
  if( zVecSizeNC(v) != zArraySize(idx) ){
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );
    return;
  }
  zQuickSort( zVecBufNC(idx), zArraySize(idx), sizeof(int), _zVecSortCmp, zVecBufNC(v) );
}

/* check if two vectors are equal. */
bool zVecIsEqual(zVec v1, zVec v2, double tol)
{
  register int i;

  if( !zVecSizeIsEqual( v1, v2 ) ) return false;
  for( i=0; i<zVecSizeNC(v1); i++ )
    if( !zIsEqual( zVecElemNC(v1,i), zVecElemNC(v2,i), tol ) ) return false;
  return true;
}

/* check if a vector is tiny. */
bool zVecIsTol(zVec v, double tol)
{
  return zRawVecIsTol( zVecBufNC(v), zVecSizeNC(v), tol );
}

/* check if a vector includes NaN. */
bool zVecIsNan(zVec v)
{
  return zRawVecIsNan( zVecBufNC(v), zVecSizeNC(v) );
}

/* add two vectors without checking size consistency. */
zVec zVecAddNC(zVec v1, zVec v2, zVec v)
{
  zRawVecAdd( zVecBufNC(v1), zVecBufNC(v2), zVecBufNC(v), zVecSizeNC(v) );
  return v;
}

/* subtract a vector from another without checking size consistency. */
zVec zVecSubNC(zVec v1, zVec v2, zVec v)
{
  zRawVecSub( zVecBufNC(v1), zVecBufNC(v2), zVecBufNC(v), zVecSizeNC(v) );
  return v;
}

/* reverse a vector without checking size consistency. */
zVec zVecRevNC(zVec v1, zVec v)
{
  zRawVecRev( zVecBufNC(v1), zVecBufNC(v), zVecSizeNC(v) );
  return v;
}

/* multiply a vector by a scalar value without checking size consistency. */
zVec zVecMulNC(zVec v1, double k, zVec v)
{
  zRawVecMul( zVecBufNC(v1), k, zVecBufNC(v), zVecSizeNC(v) );
  return v;
}

/* divide a vector by a scalar value without checking size consistency. */
zVec zVecDivNC(zVec v1, double k, zVec v)
{
  zRawVecDiv( zVecBufNC(v1), k, zVecBufNC(v), zVecSizeNC(v) );
  return v;
}

/* amplify each component of a vector by the corresponding component of
 * another vector without checking size consistency. */
zVec zVecAmpNC(zVec v1, zVec amp, zVec v)
{
  zRawVecAmp( zVecBufNC(v1), zVecBufNC(amp), zVecBufNC(v), zVecSizeNC(v) );
  return v;
}

/* demagnify each component of a vector by the corresponding component of
 * another vector without checking size consistency. */
zVec zVecDemNC(zVec v1, zVec dem, zVec v)
{
  zRawVecDem( zVecBufNC(v1), zVecBufNC(dem), zVecBufNC(v), zVecSizeNC(v) );
  return v;
}

/* concatenate a vector with another vector multiplied by a scalar value without checking size consistency. */
zVec zVecCatNC(zVec v1, double k, zVec v2, zVec v)
{
  zRawVecCat( zVecBufNC(v1), k, zVecBufNC(v2), zVecBufNC(v), zVecSizeNC(v) );
  return v;
}

#define __z_vec_size_check_2(v1,v2) \
  if( !zVecSizeIsEqual(v1,v2) ){\
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );\
    return NULL;\
  }
#define __z_vec_size_check_3(v1,v2,v) \
  if( !zVecSizeIsEqual(v1,v2) || !zVecSizeIsEqual(v1,v) ){\
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );\
    return NULL;\
  }

/* add two vectors. */
zVec zVecAdd(zVec v1, zVec v2, zVec v)
{
  __z_vec_size_check_3( v1, v2, v );
  return zVecAddNC( v1, v2, v );
}

/* subtract a vector from another. */
zVec zVecSub(zVec v1, zVec v2, zVec v)
{
  __z_vec_size_check_3( v1, v2, v );
  return zVecSubNC( v1, v2, v );
}

/* reverse a vector. */
zVec zVecRev(zVec v1, zVec v)
{
  __z_vec_size_check_2( v1, v );
  return zVecRevNC( v1, v );
}

/* multiply a vector by a scalar value. */
zVec zVecMul(zVec v1, double k, zVec v)
{
  __z_vec_size_check_2( v1, v );
  return zVecMulNC( v1, k, v );
}

/* divide a vector by a scalar value. */
zVec zVecDiv(zVec v1, double k, zVec v)
{
  __z_vec_size_check_2( v1, v );
  if( zIsTiny( k ) ){
    ZRUNERROR( ZM_ERR_ZERODIV );
    return NULL;
  }
  return zVecDivNC( v1, k, v );
}

/* amplify each component of a vector by the corresponding component of another vector. */
zVec zVecAmp(zVec v1, zVec amp, zVec v)
{
  __z_vec_size_check_3( v1, amp, v );
  return zVecAmpNC( v1, amp, v );
}

/* demagnify each component of a vector by the corresponding component of another vector. */
zVec zVecDem(zVec v1, zVec dem, zVec v)
{
  __z_vec_size_check_3( v1, dem, v );
  return zVecDemNC( v1, dem, v );
}

/* concatenate a vector with another vector multiplied by a scalar value. */
zVec zVecCat(zVec v1, double k, zVec v2, zVec v)
{
  __z_vec_size_check_3( v1, v2, v );
  return zVecCatNC( v1, k, v2, v );
}

/* concatenate multiple vectors multiplied by scalar values.
 * this function is internally called to mange variable arguments.
 */
static void _zVecCats(zVec v, int n, va_list args)
{
  register int i;
  double k;
  zVec vec;

  for( i=0; i<n; i++ ){
    k = (double)va_arg( args, double );
    vec = (zVec)va_arg( args, zVec );
    if( !zVecSizeIsEqual( v, vec ) )
      ZRUNWARN( ZM_ERR_SIZMIS_VEC );
    zVecCatDRC( v, k, vec );
  }
}

/* concatenate multiple vectors multiplied by scalar values. */
zVec zVecCats(zVec v, int n, ...)
{
  va_list args;

  va_start( args, n );
  _zVecCats( v, n, args );
  va_end( args );
  return v;
}

/* linear sum of multiple vectors multiplied by scalar values. */
zVec zVecLS(zVec v, int n, ...)
{
  va_list args;

  zVecZero( v );
  va_start( args, n );
  _zVecCats( v, n, args );
  va_end( args );
  return v;
}

/* interior division of two vectors. */
zVec zVecInterDiv(zVec v1, zVec v2, double ratio, zVec v)
{
  zRawVecInterDiv( zVecBufNC(v1), zVecBufNC(v2), ratio, zVecBufNC(v), zVecSizeNC(v) );
  return v;
}

/* midpoint of two vectors. */
zVec zVecMid(zVec v1, zVec v2, zVec v)
{
  zRawVecMid( zVecBufNC(v1), zVecBufNC(v2), zVecBufNC(v), zVecSizeNC(v) );
  return v;
}

/* scale a raw vector with two boundary vectors. */
zVec zVecScale(zVec x, zVec min, zVec max, zVec v)
{
  zRawVecScale( zVecBuf(x), zVecBuf(min), zVecBuf(max), zVecBuf(v), zVecSizeNC(x) );
  return v;
}

/* inner product of two vectors without checking size consistency. */
double zVecInnerProdNC(zVec v1, zVec v2)
{
  return zRawVecInnerProd( zVecBufNC(v1), zVecBufNC(v2), zVecSizeNC(v1) );
}

/* inner product of two vectors. */
double zVecInnerProd(zVec v1, zVec v2)
{
  if( !zVecSizeIsEqual(v1,v2) ){
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );
    return 0;
  }
  return zVecInnerProdNC( v1, v2 );
}

/* squared norm of a vector. */
double zVecSqrNorm(zVec v)
{
  return zVecInnerProdNC( v, v );
}

/* weighted squared norm of a vector without checking size consistency. */
double zVecWSqrNormNC(zVec v, zVec w)
{
  return zRawVecWSqrNorm( zVecBufNC(v), zVecBufNC(w), zVecSizeNC(v) );
}

/* weighted squared norm of a vector. */
double zVecWSqrNorm(zVec v, zVec w)
{
  if( !zVecSizeIsEqual(v,w) ){
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );
    return 0;
  }
  return zVecWSqrNormNC( v, w );
}

/* infinity norm of a vector. */
double zVecInfNorm(zVec v)
{
  return zDataAbsMax( zVecBufNC(v), zVecSizeNC(v), NULL );
}

/* normalize a vector. */
zVec zVecNormalize(zVec src, zVec dest)
{
  return zRawVecNormalize( zVecBufNC(src), zVecSizeNC(src), zVecBufNC(dest) ) ? dest : NULL;
}

/* ********************************************************** */
/* I/O
 * ********************************************************** */

/* read a vector from a ZTK format processor. */
zVec zVecFromZTK(ZTK *ztk)
{
  register int i, size;
  zVec v;

  size = ZTKInt(ztk);
  if( !( v = zVecAlloc( size ) ) ) return NULL;
  for( i=0; i<size; i++ )
    zVecSetElemNC( v, i, ZTKDouble(ztk) );
  return v;
}

/* scan a vector from a file. */
zVec zVecFScan(FILE *fp)
{
  register int i;
  int size;
  zVec v;

  if( !zFInt( fp, &size ) ){
    ZRUNERROR( ZM_ERR_SIZUNFOUND_VEC );
    return NULL;
  }
  if( !( v = zVecAlloc( size ) ) ) return NULL;
  for( i=0; i<size; i++ ){
    if( zFDouble( fp, &zVecElemNC(v,i) ) ){
      ZRUNERROR( ZM_WARN_VEC_SIZMIS, i, size );
      break;
    }
  }
  return v;
}

/* print a vector out to a file. */
void zVecFPrint(FILE *fp, zVec v)
{
  register int i;

  if( !v )
    fprintf( fp, "(null vector)\n" );
  else{
    fprintf( fp, "%d (", zVecSizeNC(v) );
    for( i=0; i<zVecSizeNC(v); i++ )
      fprintf( fp, " %.10g", zVecElemNC(v,i) );
    fprintf( fp, " )\n" );
  }
}

/* print a vector out to a file. */
void zVecDataFPrint(FILE *fp, zVec v)
{
  register int i;

  if( !v ) return;
  for( i=0; i<zVecSizeNC(v); i++ )
    fprintf( fp, "%.10g ", zVecElemNC(v,i) );
  fprintf( fp, "\n" );
}
