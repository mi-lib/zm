/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_cvec - complex vector class.
 */

#include <zm/zm_cvec.h>

/* allocate a complex vector. */
zCVec zCVecAlloc(int size)
{
  zCVec v;

  if( !( v = zAlloc( zCVecStruct, 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  if( !( zArrayBuf(v) = zAlloc( zComplex, size ) ) ){
    ZALLOCERROR();
    free( v );
    return NULL;
  }
  zCVecSetSizeNC( v, size );
  return v;
}

/* free a complex vector. */
void zCVecFree(zCVec v)
{
  if( !v ) return;
  zFree( zArrayBuf(v) );
  free( v );
}

/* zero a complex vector. */
zCVec zCVecZero(zCVec v)
{
  int i;

  for( i=0; i<zCVecSizeNC(v); i++ )
    zComplexZero( zCVecElemNC( v, i ) );
  return v;
}

/* touchup a complex vector. */
zCVec zCVecTouchup(zCVec v, double tol)
{
  int i;

  for( i=0; i<zCVecSizeNC(v); i++ )
    zComplexTouchup( zCVecElemNC(v,i), tol );
  return v;
}

/* copy a complex vector without checking size consistency. */
zCVec zCVecCopyNC(const zCVec src, zCVec dest)
{
  int i;

  for( i=0; i<zCVecSizeNC(dest); i++ )
    zComplexCopy( zCVecElemNC(src,i), zCVecElemNC(dest,i) );
  return dest;
}

/* copy a complex vector. */
zCVec zCVecCopy(const zCVec src, zCVec dest)
{
  return zCVecSizeEqual( src, dest ) ? zCVecCopyNC( src, dest ) : NULL;
}

/* clone a complex vector. */
zCVec zCVecClone(const zCVec src)
{
  zCVec dest;

  if( !src ) return NULL;
  if( ( dest = zCVecAlloc( zCVecSizeNC(src) ) ) )
    zCVecCopyNC( src, dest );
  return dest;
}

/* convert a real vector to a complex vector. */
zCVec zVecToCVec(const zVec v, zCVec cv)
{
  int i;

  if( zVecSize(v) != zCVecSize(cv) ){
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );
    return NULL;
  }
  for( i=0; i<zCVecSizeNC(cv); i++ )
    zComplexCreate( zCVecElemNC(cv,i), zVecElemNC(v,i), 0 );
  return cv;
}

/* abstract the real-part vector from a complex vector. */
zVec zCVecToReVec(const zCVec cv, zVec rv)
{
  int i;

  if( zVecSize(rv) != zCVecSize(cv) ){
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );
    return NULL;
  }
  for( i=0; i<zCVecSizeNC(cv); i++ )
    zVecSetElemNC( rv, i, zCVecElemNC(cv,i)->re );
  return rv;
}

/* abstract the imaginary-part vector from a complex vector. */
zVec zCVecToImVec(const zCVec cv, zVec iv)
{
  int i;

  if( zVecSize(iv) != zCVecSize(cv) ){
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );
    return NULL;
  }
  for( i=0; i<zCVecSizeNC(cv); i++ )
    zVecSetElemNC( iv, i, zCVecElemNC(cv,i)->im );
  return iv;
}

/* generate a uniformly random complex vector. */
zCVec zCVecRandUniform(zCVec v, double rmin, double imin, double rmax, double imax)
{
  int i;

  for( i=0; i<zCVecSizeNC(v); i++ )
    zComplexCreate( zCVecElemNC(v,i), zRandF(rmin,rmax), zRandF(imin,imax) );
  return v;
}

/* check if two complex vectors are equal. */
bool zCVecEqual(const zCVec v1, const zCVec v2, double tol)
{
  int i;

  if( !zCVecSizeEqual( v1, v2 ) ) return false;
  for( i=0; i<zCVecSizeNC(v1); i++ )
    if( !zComplexEqual( zCVecElemNC(v1,i), zCVecElemNC(v2,i), tol ) ) return false;
  return true;
}

/* check if a complex vector is tiny. */
bool zCVecIsTol(const zCVec v, double tol)
{
  int i;

  for( i=0; i<zCVecSizeNC(v); i++ )
    if( !zComplexIsTol( zCVecElemNC(v,i), tol ) ) return false;
  return true;
}

/* split a complex vector into real and imaginary vectors. */
bool zCVecToReImVec(const zCVec cvec, zVec *rvec, zCVec *ivec, double tol)
{
  zIndex ridx, iidx;
  int i, rsize, isize;
  bool ret = true;

  ridx = zIndexCreate( zCVecSizeNC(cvec) );
  iidx = zIndexCreate( zCVecSizeNC(cvec) );
  if( !ridx || !iidx ){
    ret = false;
    goto TERMINATE;
  }
  for( i=rsize=isize=0; i<zCVecSizeNC(cvec); i++ ){
    if( zComplexIsReal( zCVecElemNC(cvec,i), tol ) )
      zIndexSetElemNC( ridx, rsize++, i );
    else{
      zIndexSetElemNC( iidx, isize++, i );
    }
  }
  *rvec = rsize > 0 ? zVecAlloc( ( zIndexSizeNC(ridx) = rsize ) ) : NULL;
  *ivec = isize > 0 ? zCVecAlloc( ( zIndexSizeNC(iidx) = isize ) ) : NULL;
  if( !*rvec && !*ivec ){
    zVecFree( *rvec );
    zCVecFree( *ivec );
    ret = false;
    goto TERMINATE;
  }
  for( i=0; i<rsize; i++ )
    zVecSetElemNC( *rvec, i, zCVecElemNC(cvec,zIndexElemNC(ridx,i))->re );
  for( i=0; i<isize; i++ )
    zCVecSetElemNC( *ivec, i, zCVecElemNC(cvec,zIndexElemNC(iidx,i)) );
  if( *ivec && !zCVecConjPair( *ivec, tol ) ) ret = false;

 TERMINATE:
  zIndexFree( ridx );
  zIndexFree( iidx );
  return ret;
}

/* reorder a complex vector as co-conjugate numbers are paired as adjacencies. */
zCVec zCVecConjPair(zCVec v, double tol)
{
  int i, j;

  if( zIsOdd( zCVecSizeNC(v) ) ){
    ZRUNERROR( ZM_ERR_CVEC_CONJPAIR_UNABLE );
    return NULL;
  }
  for( i=0; i<zCVecSizeNC(v); i+=2 ){
    for( j=i+1; j<zCVecSizeNC(v); j++ ){
      if( zComplexCoconj( zCVecElemNC(v,i), zCVecElemNC(v,j), tol ) ){
        if( j > i + 1 ) zSwap( zComplex, *zCVecElemNC(v,i+1), *zCVecElemNC(v,j) );
        break;
      }
    }
    if( j == zCVecSizeNC(v) ){
      ZRUNERROR( ZM_ERR_CVEC_CONJPAIR_UNABLE );
      return NULL;
    }
  }
  return v;
}

/* basic arithmetics for the complex vector. */

/* add two complex vectors without checking size consistency. */
zCVec zCVecAddNC(const zCVec v1, const zCVec v2, zCVec v)
{
  int i;

  for( i=0; i<zCVecSizeNC(v); i++ )
    zComplexAdd( zCVecElemNC(v1,i), zCVecElemNC(v2,i), zCVecElemNC(v,i) );
  return v;
}

/* subtract a complex vector from another without checking size consistency. */
zCVec zCVecSubNC(const zCVec v1, const zCVec v2, zCVec v)
{
  int i;

  for( i=0; i<zCVecSizeNC(v); i++ )
    zComplexSub( zCVecElemNC(v1,i), zCVecElemNC(v2,i), zCVecElemNC(v,i) );
  return v;
}

/* reverse a complex vector without checking size consistency. */
zCVec zCVecRevNC(const zCVec v1, zCVec v)
{
  int i;

  for( i=0; i<zCVecSizeNC(v); i++ )
    zComplexRev( zCVecElemNC(v1,i), zCVecElemNC(v,i) );
  return v;
}

/* multiply a complex vector by a complex scalar value without checking size consistency. */
zCVec zCVecMulNC(const zCVec v1, double k, zCVec v)
{
  int i;

  for( i=0; i<zCVecSizeNC(v); i++ )
    zComplexMul( zCVecElemNC(v1,i), k, zCVecElemNC(v,i) );
  return v;
}

/* divide a complex vector by a complex scalar value without checking size consistency. */
zCVec zCVecDivNC(const zCVec v1, double k, zCVec v)
{
  int i;

  k = 1 / k;
  for( i=0; i<zCVecSizeNC(v); i++ )
    zComplexMul( zCVecElemNC(v1,i), k, zCVecElemNC(v,i) );
  return v;
}

/* multiply a complex vector by a complex scalar value without checking size consistency. */
zCVec zCVecCMulNC(const zCVec v1, const zComplex *z, zCVec v)
{
  int i;

  for( i=0; i<zCVecSizeNC(v); i++ )
    zComplexCMul( zCVecElemNC(v1,i), z, zCVecElemNC(v,i) );
  return v;
}

/* divide a complex vector by a complex scalar value without checking size consistency. */
zCVec zCVecCDivNC(const zCVec v1, const zComplex *z, zCVec v)
{
  int i;
  double r;
  zComplex dz;

  r = zComplexSqrAbs( z );
  zComplexConj( z, &dz );
  zComplexDiv( &dz, r, &dz );
  for( i=0; i<zCVecSizeNC(v); i++ )
    zComplexCMul( zCVecElemNC(v1,i), &dz, zCVecElemNC(v,i) );
  return v;
}

/* concatenate a complex vector with another multiplied by a complex scalar value without checking
 * size consistency. */
zCVec zCVecCatNC(const zCVec v1, const zComplex *z, const zCVec v2, zCVec v)
{
  int i;
  zComplex dz;

  for( i=0; i<zCVecSizeNC(v); i++ ){
    zComplexCMul( zCVecElemNC(v2,i), z, &dz );
    zComplexAdd( zCVecElemNC(v1,i), &dz, zCVecElemNC(v,i) );
  }
  return v;
}

#define __z_cvec_size_check_2(v1,v2) \
  if( !zCVecSizeEqual(v1,v2) ){\
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );\
    return NULL;\
  }
#define __z_cvec_size_check_3(v1,v2,v) \
  if( !zCVecSizeEqual(v1,v2) || !zCVecSizeEqual(v1,v) ){\
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );\
    return NULL;\
  }

/* add two complex vectors. */
zCVec zCVecAdd(const zCVec v1, const zCVec v2, zCVec v)
{
  __z_cvec_size_check_3( v1, v2, v );
  return zCVecAddNC( v1, v2, v );
}

/* subtract a complex vector from another. */
zCVec zCVecSub(const zCVec v1, const zCVec v2, zCVec v)
{
  __z_cvec_size_check_3( v1, v2, v );
  return zCVecSubNC( v1, v2, v );
}

/* reverse a complex vector. */
zCVec zCVecRev(const zCVec v1, zCVec v)
{
  __z_cvec_size_check_2( v1, v );
  return zCVecRevNC( v1, v );
}

/* multiply a complex vector by a scalar value. */
zCVec zCVecMul(const zCVec v1, double k, zCVec v)
{
  __z_cvec_size_check_2( v1, v );
  return zCVecMulNC( v1, k, v );
}

/* divide a complex vector by a scalar value. */
zCVec zCVecDiv(const zCVec v1, double k, zCVec v)
{
  __z_cvec_size_check_2( v1, v );
  if( zIsTiny( k ) ){
    ZRUNERROR( ZM_ERR_ZERODIV );
    return NULL;
  }
  return zCVecDivNC( v1, k, v );
}

/* multiply a complex vector by a complex number. */
zCVec zCVecCMul(const zCVec v1, const zComplex *z, zCVec v)
{
  __z_cvec_size_check_2( v1, v );
  return zCVecCMulNC( v1, z, v );
}

/* divide a complex vector by a complex number. */
zCVec zCVecCDiv(const zCVec v1, const zComplex *z, zCVec v)
{
  __z_cvec_size_check_2( v1, v );
  if( zComplexIsTiny( z ) ){
    ZRUNERROR( ZM_ERR_ZERODIV );
    return NULL;
  }
  return zCVecCDivNC( v1, z, v );
}

/* concatenate a complex vector with another multiplied by a complex scalar value. */
zCVec zCVecCat(const zCVec v1, const zComplex *z, const zCVec v2, zCVec v)
{
  __z_cvec_size_check_3( v1, v2, v );
  return zCVecCatNC( v1, z, v2, v );
}

/* inner product of two complex vectors without checking size consistency. */
zComplex *zCVecInnerProdNC(const zCVec v1, const zCVec v2, zComplex *z)
{
  int i;
  zComplex c;

  zComplexZero( z );
  for( i=0; i<zCVecSizeNC(v1); i++ ){
    zComplexCMulConj( zCVecElemNC(v1,i), zCVecElemNC(v2,i), &c );
    zComplexAdd( z, &c, z );
  }
  return z;
}

/* inner product of two complex vector. */
zComplex *zCVecInnerProd(const zCVec v1, const zCVec v2, zComplex *z)
{
  if( !zCVecSizeEqual(v1,v2) ){
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );
    return NULL;
  }
  return zCVecInnerProdNC( v1, v2, z );
}

/* squared norm of a complex vector. */
double zCVecSqrNorm(const zCVec v)
{
  zComplex z;

  return zComplexAbs( zCVecInnerProdNC( v, v, &z ) );
}

/* normalize a complex vector. */
zCVec zCVecNormalize(const zCVec src, zCVec dest)
{
  int i;
  double r;

  if( zIsTiny( ( r = zCVecNorm( src ) ) ) ){
    ZRUNWARN( ZM_ERR_VEC_ZERONORM );
    return NULL;
  }
  for( i=0; i<zCVecSizeNC(dest); i++ )
    zComplexDiv( zCVecElemNC(src,i), r, zCVecElemNC(dest,i) );
  return dest;
}

/* print a complex vector to a file. */
void zCVecFPrint(FILE *fp, const zCVec v)
{
  int i;

  if( !v )
    fprintf( fp, "(null vector)\n" );
  else{
    fprintf( fp, "%d (\n", zCVecSizeNC(v) );
    for( i=0; i<zCVecSizeNC(v); i++ ){
      fprintf( fp, "  " );
      zComplexFPrint( fp, zCVecElemNC(v,i) );
      fprintf( fp, "\n" );
    }
    fprintf( fp, ")\n" );
  }
}
