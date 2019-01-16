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

static int _zVecSortCmp(void *p1, void *p2, void *priv);
static void _zVecCats(zVec v, int n, va_list args);

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
  if( !( zVecBuf(v) = zAlloc( double, size ) ) ){
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
  zFree( zVecBuf( v ) );
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

/* cleanup a vector. */
zVec zVecClear(zVec v)
{
  zRawVecClear( zVecBuf(v), zVecSizeNC(v) );
  return v;
}

/* touchup a vector. */
zVec zVecTouchup(zVec v)
{
  zRawVecTouchup( zVecBuf(v), zVecSizeNC(v) );
  return v;
}

/* copy a vector without checking size consistency. */
zVec zVecCopyNC(zVec src, zVec dest)
{
  zRawVecCopy( zVecBuf(src), zVecBuf(dest), zVecSizeNC(dest) );
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
  zRawVecCopy( array, zVecBuf(v), s );
  return v;
}

/* clone a vector. */
zVec zVecClone(zVec src)
{
  zVec dest;

  if( ( dest = zVecAlloc( zVecSize(src) ) ) )
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
  zRawVecGet( zVecBuf(src), pos, zVecBuf(dest), zVecSize(dest) );
  return dest;
}

/* get a parttial vector. */
zVec zVecGet(zVec src, int pos, zVec dest)
{
  if( pos+zVecSize(dest) > zVecSize(src) ){
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );
    return NULL;
  }
  return zVecGetNC( src, pos, dest );
}

/* put a partial vector without checking the size validity. */
zVec zVecPutNC(zVec dest, int pos, zVec src)
{
  zRawVecPut( zVecBuf(dest), pos, zVecBuf(src), zVecSize(src) );
  return dest;
}

/* put a partial vector. */
zVec zVecPut(zVec dest, int pos, zVec src)
{
  if( pos+zVecSize(src) > zVecSize(dest) ){
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );
    return NULL;
  }
  return zVecPutNC( dest, pos, src );
}

/* create a uniform vector. */
zVec zVecSetAll(zVec v, double val)
{
  zRawVecSetAll( zVecBuf(v), zVecSizeNC(v), val );
  return v;
}

/* create a linear space vector. */
zVec zVecLinSpace(zVec v, double from, double to)
{
  zRawVecLinSpace( zVecBuf(v), zVecSizeNC(v), from, to );
  return v;
}

/* create a random vector with a uniform range. */
zVec zVecRandUniform(zVec v, double min, double max)
{
  zRawVecRandUniform( zVecBuf(v), zVecSizeNC(v), min, max );
  return v;
}

/* create a random vector with range vectors. */
zVec zVecRand(zVec v, zVec min, zVec max)
{
  zRawVecRand( zVecBuf(v), zVecBuf(min), zVecBuf(max), zVecSizeNC(v) );
  return v;
}

/* shift vector by a constant scalar value. */
zVec zVecShift(zVec v, double shift)
{
  zRawVecShift( zVecBuf(v), zVecSizeNC(v), shift );
  return v;
}

/* swap vector components without checking size. */
zVec zVecSwapNC(zVec v, int i1, int i2)
{
  zRawVecSwap( zVecBuf(v), i1, i2 );
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
int _zVecSortCmp(void *p1, void *p2, void *priv)
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
  zQuickSort( zArrayBuf(idx), zArraySize(idx), sizeof(int), _zVecSortCmp, zVecBuf(v) );
}

/* check if two vectors are equal. */
bool zVecIsEqual(zVec v1, zVec v2)
{
  register int i;

  if( !zVecSizeIsEqual( v1, v2 ) ) return false;
  for( i=0; i<zVecSizeNC(v1); i++ )
    if( !zIsTiny( zVecElem(v1,i) - zVecElem(v2,i) ) ) return false;
  return true;
}

/* check if a vector is tiny. */
bool zVecIsTol(zVec v, double tol)
{
  return zRawVecIsTol( zVecBuf(v), zVecSizeNC(v), tol );
}

/* check if a vector includes NaN. */
bool zVecIsNan(zVec v)
{
  return zRawVecIsNan( zVecBuf(v), zVecSizeNC(v) );
}

/* add two vectors without checking size consistency. */
zVec zVecAddNC(zVec v1, zVec v2, zVec v)
{
  zRawVecAdd( zVecBuf(v1), zVecBuf(v2), zVecBuf(v), zVecSizeNC(v) );
  return v;
}

/* subtract a vector from another without checking size consistency. */
zVec zVecSubNC(zVec v1, zVec v2, zVec v)
{
  zRawVecSub( zVecBuf(v1), zVecBuf(v2), zVecBuf(v), zVecSizeNC(v) );
  return v;
}

/* reverse a vector without checking size consistency. */
zVec zVecRevNC(zVec v1, zVec v)
{
  zRawVecRev( zVecBuf(v1), zVecBuf(v), zVecSizeNC(v) );
  return v;
}

/* multiply a vector by a scalar value without checking size consistency. */
zVec zVecMulNC(zVec v1, double k, zVec v)
{
  zRawVecMul( zVecBuf(v1), k, zVecBuf(v), zVecSizeNC(v) );
  return v;
}

/* divide a vector by a scalar value without checking size consistency. */
zVec zVecDivNC(zVec v1, double k, zVec v)
{
  zRawVecDiv( zVecBuf(v1), k, zVecBuf(v), zVecSizeNC(v) );
  return v;
}

/* amplify each component of a vector by the corresponding component of
 * another vector without checking size consistency. */
zVec zVecAmpNC(zVec v1, zVec amp, zVec v)
{
  zRawVecAmp( zVecBuf(v1), zVecBuf(amp), zVecBuf(v), zVecSizeNC(v) );
  return v;
}

/* demagnify each component of a vector by the corresponding component of
 * another vector without checking size consistency. */
zVec zVecDemNC(zVec v1, zVec dem, zVec v)
{
  zRawVecDem( zVecBuf(v1), zVecBuf(dem), zVecBuf(v), zVecSizeNC(v) );
  return v;
}

/* concatenate a vector with another vector multiplied by a scalar value without checking size consistency. */
zVec zVecCatNC(zVec v1, double k, zVec v2, zVec v)
{
  zRawVecCat( zVecBuf(v1), k, zVecBuf(v2), zVecBuf(v), zVecSizeNC(v) );
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
void _zVecCats(zVec v, int n, va_list args)
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

  zVecClear( v );
  va_start( args, n );
  _zVecCats( v, n, args );
  va_end( args );
  return v;
}

/* inner product of two vectors without checking size consistency. */
double zVecInnerProdNC(zVec v1, zVec v2)
{
  return zRawVecInnerProd( zVecBuf(v1), zVecBuf(v2), zVecSizeNC(v1) );
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
  return zRawVecWSqrNorm( zVecBuf(v), zVecBuf(w), zVecSizeNC(v) );
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
  return zDataAbsMax( zVecBuf(v), zVecSizeNC(v), NULL );
}

/* normalize a vector. */
zVec zVecNormalize(zVec src, zVec dest)
{
  return zRawVecNormalize( zVecBuf(src), zVecSizeNC(src), zVecBuf(dest) ) ? dest : NULL;
}

/* read information of a vector from file. */
zVec zVecFRead(FILE *fp)
{
  register int i, size;
  zVec v;

  size = zFInt( fp );
  if( !( v = zVecAlloc( size ) ) ) return NULL;
  for( i=0; i<size; i++ )
    zVecSetElemNC( v, i, zFDouble( fp ) );
  return v;
}

/* output information of a vector to file. */
void zVecFWrite(FILE *fp, zVec v)
{
  register int i;

  if( !v )
    fprintf( fp, "(null vector)\n" );
  else{
    fprintf( fp, "%d (", zVecSizeNC(v) );
    for( i=0; i<zVecSizeNC(v); i++ )
      fprintf( fp, " %.10g", zVecElem(v,i) );
    fprintf( fp, " )\n" );
  }
}

/* output information of a vector to file. */
void zVecDataFWrite(FILE *fp, zVec v)
{
  register int i;

  if( !v ) return;
  for( i=0; i<zVecSizeNC(v); i++ )
    fprintf( fp, "%.10g ", zVecElem(v,i) );
  fprintf( fp, "\n" );
}
