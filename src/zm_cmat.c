/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_cmat - complex matrix class.
 */

#include <zm/zm_cmat.h>

/* allocate a complex matrix. */
zCMat zCMatAlloc(int row, int col)
{
  zCMat m;

  if( !( m = zAlloc( zCMatStruct, 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  if( !( zCMatBufNC(m) = zAlloc( zComplex, row*col ) ) ){
    ZALLOCERROR();
    zFree( m );
    return NULL;
  }
  zCMatSetSizeNC( m, row, col );
  return m;
}

/* free a complex matrix. */
void zCMatFree(zCMat m)
{
  if( m ){
    zFree( zCMatBufNC(m) );
    free( m );
  }
}

/* zero a complex matrix. */
zCMat zCMatZero(zCMat m)
{
  int i, n;

  n = zCMatRowSizeNC(m) * zCMatColSizeNC(m);
  for( i=0; i<n; i++ )
    zComplexZero( &zCMatBufNC(m)[i] );
  return m;
}

/* touchup a complex matrix. */
zCMat zCMatTouchup(zCMat m)
{
  int i, j;

  for( i=0; i<zCMatRowSizeNC(m); i++ )
    for( j=0; j<zCMatColSizeNC(m); j++ )
      zComplexTouchup( zCMatElemNC(m,i,j) );
  return m;
}

/* create a random complex matrix with a uniform range. */
zCMat zCMatRandUniform(zCMat m, double rmin, double rmax, double imin, double imax)
{
  int i, j;

  for( i=0; i<zCMatRowSizeNC(m); i++ )
    for( j=0; j<zCMatColSizeNC(m); j++ )
      zComplexCreate( zCMatElemNC(m,i,j), zRandF(rmin,rmax), zRandF(imin,imax) );
  return m;
}

/* create a complex random matrix with range matrices. */
zCMat zCMatRand(zCMat m, zCMat min, zCMat max)
{
  int i, j;

  for( i=0; i<zCMatRowSizeNC(m); i++ )
    for( j=0; j<zCMatColSizeNC(m); j++ )
      zComplexCreate( zCMatElemNC(m,i,j),
        zRandF(zCMatElemNC(min,i,j)->re,zCMatElemNC(max,i,j)->re),
        zRandF(zCMatElemNC(min,i,j)->im,zCMatElemNC(max,i,j)->im) );
  return m;
}

/* copy a complex matrix to another without checking size consistency. */
zCMat zCMatCopyNC(const zCMat src, zCMat dest)
{
  int i, n;

  n = zCMatRowSizeNC(src) * zCMatColSizeNC(src);
  for( i=0; i<n; i++ )
    zComplexCopy( &zCMatBufNC(src)[i], &zCMatBufNC(dest)[i] );
  return dest;
}

/* copy a complex matrix. */
zCMat zCMatCopy(const zCMat src, zCMat dest)
{
  if( !zCMatSizeEqual( src, dest ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  return zCMatCopyNC( src, dest );
}

/* clone a complex matrix. */
zCMat zCMatClone(const zCMat src)
{
  zCMat dest;

  if( !src ) return NULL;
  if( ( dest = zCMatAlloc( zCMatRowSizeNC(src), zCMatColSizeNC(src) ) ) )
    zCMatCopyNC( src, dest );
  return dest;
}

/* convert a real matrix to a complex matrix. */
zCMat zMatToCMat(const zMat m, zCMat cm)
{
  int i, n;

  n = zMatRowSizeNC(m) * zMatColSizeNC(m);
  for( i=0; i<n; i++ )
    zComplexCreate( &zCMatBufNC(cm)[i], zMatBufNC(m)[i], 0 );
  return cm;
}

/* abstract row vector of a complex matrix without checking size consistency. */
zCVec zCMatGetRowNC(const zCMat m, int row, zCVec v)
{
  memcpy( zCVecBufNC(v), zCMatBufNC(m)+row*zCMatColSize(m), sizeof(zComplex)*zCMatColSizeNC(m) );
  return v;
}

/* abstract column vector of a complex matrix without checking size consistency. */
zCVec zCMatGetColNC(const zCMat m, int col, zCVec v)
{
  int i;

  for( i=0; i<zCMatRowSizeNC(m); i++ )
    zComplexCopy( zCMatElemNC(m,i,col), zCVecElemNC(v,i) );
  return v;
}

/* abstract row vector of a complex matrix. */
zCVec zCMatGetRow(const zCMat m, int row, zCVec v)
{
  if( !zCMatColCVecSizeEqual( m, v ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  if( row >= zCMatRowSizeNC(m) ){
    ZRUNERROR( ZM_ERR_INVALID_ROW );
    return NULL;
  }
  return zCMatGetRowNC( m, row, v );
}

/* abstract column vector of a complex matrix. */
zCVec zCMatGetCol(const zCMat m, int col, zCVec v)
{
  if( !zCMatRowCVecSizeEqual( m, v ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  if( col >= zCMatColSizeNC(m) ){
    ZRUNERROR( ZM_ERR_INVALID_COL );
    return NULL;
  }
  return zCMatGetColNC( m, col, v );
}

/* put a row vector to a complex matrix without checking size consistency. */
zCMat zCMatPutRowNC(zCMat m, int row, const zCVec v)
{
  memcpy( zCMatBufNC(m)+row*zCMatColSize(m), zCVecBufNC(v), sizeof(zComplex)*zCVecSizeNC(v) );
  return m;
}

/* put a column vector to a complex matrix without checking size consistency. */
zCMat zCMatPutColNC(zCMat m, int col, const zCVec v)
{
  int i;

  for( i=0; i<zCVecSizeNC(v); i++ )
    zComplexCopy( zCVecElemNC(v,i), zCMatElemNC(m,i,col) );
  return m;
}

/* put a row vector to a complex matrix. */
zCMat zCMatPutRow(zCMat m, int row, const zCVec v)
{
  if( !zCMatColCVecSizeEqual( m, v ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  if( row >= zCMatRowSizeNC(m) ){
    ZRUNERROR( ZM_ERR_INVALID_ROW );
    return NULL;
  }
  return zCMatPutRowNC( m, row, v );
}

/* put a column vector to a complex matrix. */
zCMat zCMatPutCol(zCMat m, int col, const zCVec v)
{
  if( !zCMatRowCVecSizeEqual( m, v ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  if( col >= zCMatColSizeNC(m) ){
    ZRUNERROR( ZM_ERR_INVALID_COL );
    return NULL;
  }
  return zCMatPutColNC( m, col, v );
}

/* check if a complex matrix is tiny. */
bool zCMatIsTol(const zCMat m, double tol)
{
  int i, n;

  n = zCMatRowSizeNC(m) * zCMatColSizeNC(m);
  for( i=0; i<n; i++ )
    if( !zComplexIsTol( &zCMatBufNC(m)[i], tol ) ) return false;
  return true;
}

/* add two complex matrices without checking size consistency. */
zCMat zCMatAddNC(const zCMat m1, const zCMat m2, zCMat m)
{
  int i, n;

  n = zCMatRowSizeNC(m) * zCMatColSizeNC(m);
  for( i=0; i<n; i++ )
    zComplexAdd( &zCMatBufNC(m1)[i], &zCMatBufNC(m2)[i], &zCMatBufNC(m)[i] );
  return m;
}

/* subtract a complex matrix from another without checking size consistency. */
zCMat zCMatSubNC(const zCMat m1, const zCMat m2, zCMat m)
{
  int i, n;

  n = zCMatRowSizeNC(m) * zCMatColSizeNC(m);
  for( i=0; i<n; i++ )
    zComplexSub( &zCMatBufNC(m1)[i], &zCMatBufNC(m2)[i], &zCMatBufNC(m)[i] );
  return m;
}

/* reverse a complex matrix without checking size consistency. */
zCMat zCMatRevNC(const zCMat m1, zCMat m)
{
  int i, n;

  n = zCMatRowSizeNC(m) * zCMatColSizeNC(m);
  for( i=0; i<n; i++ )
    zComplexRev( &zCMatBufNC(m1)[i], &zCMatBufNC(m)[i] );
  return m;
}

/* multiply a complex matrix by a real scalar value without checking size consistency. */
zCMat zCMatMulNC(const zCMat m1, double k, zCMat m)
{
  int i, n;

  n = zCMatRowSizeNC(m) * zCMatColSizeNC(m);
  for( i=0; i<n; i++ )
    zComplexMul( &zCMatBufNC(m1)[i], k, &zCMatBufNC(m)[i] );
  return m;
}

/* divide a complex matrix by a real scalar value without checking size consistency. */
zCMat zCMatDivNC(const zCMat m1, double k, zCMat m)
{
  int i, n;

  k = 1 / k;
  n = zCMatRowSizeNC(m) * zCMatColSizeNC(m);
  for( i=0; i<n; i++ )
    zComplexMul( &zCMatBufNC(m1)[i], k, &zCMatBufNC(m)[i] );
  return m;
}

/* multiply a complex matrix by a complex scalar without checking size consistency. */
zCMat zCMatCMulNC(const zCMat m1, const zComplex *z, zCMat m)
{
  int i, n;

  n = zCMatRowSizeNC(m) * zCMatColSizeNC(m);
  for( i=0; i<n; i++ )
    zComplexCMul( &zCMatBufNC(m1)[i], z, &zCMatBufNC(m)[i] );
  return m;
}

/* divide a complex matrix by a complex scalar without checking size consistency. */
zCMat zCMatCDivNC(const zCMat m1, const zComplex *z, zCMat m)
{
  int i, n;
  double r;
  zComplex dz;

  r = zComplexSqrAbs( z );
  zComplexConj( z, &dz );
  zComplexDiv( &dz, r, &dz );
  n = zCMatRowSizeNC(m) * zCMatColSizeNC(m);
  for( i=0; i<n; i++ )
    zComplexCMul( &zCMatBufNC(m1)[i], &dz, &zCMatBufNC(m)[i] );
  return m;
}

/* add two complex matrices. */
zCMat zCMatAdd(const zCMat m1, const zCMat m2, zCMat m)
{
  if( !zCMatSizeEqual(m1,m2) || !zCMatSizeEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  return zCMatAddNC( m1, m2, m );
}

/* substract a complex matrix from another. */
zCMat zCMatSub(const zCMat m1, const zCMat m2, zCMat m)
{
  if( !zCMatSizeEqual(m1,m2) || !zCMatSizeEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  return zCMatSubNC( m1, m2, m );
}

/* reverse a complex matrix. */
zCMat zCMatRev(const zCMat m1, zCMat m)
{
  if( !zCMatSizeEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  return zCMatRevNC( m1, m );
}

/* multiply a complex matrix by a real scalar value. */
zCMat zCMatMul(const zCMat m1, double k, zCMat m)
{
  if( !zCMatSizeEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  return zCMatMulNC( m1, k, m );
}

/* divide a complex matrix by a real scalar value. */
zCMat zCMatDiv(const zCMat m1, double k, zCMat m)
{
  if( !zCMatSizeEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  if( zIsTiny( k ) ){
    ZRUNERROR( ZM_ERR_ZERODIV );
    return NULL;
  }
  return zCMatDivNC( m1, k, m );
}

/* multiply a complex matrix by a complex number. */
zCMat zCMatCMul(const zCMat m1, const zComplex *z, zCMat m)
{
  if( !zCMatSizeEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  return zCMatCMulNC( m1, z, m );
}

/* divide a complex matrix by a complex number. */
zCMat zCMatCDiv(const zCMat m1, const zComplex *z, zCMat m)
{
  if( !zCMatSizeEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  if( zComplexIsTiny( z ) ){
    ZRUNERROR( ZM_ERR_ZERODIV );
    return NULL;
  }
  return zCMatCDivNC( m1, z, m );
}

/* multiply a complex matrix and a complex column vector without checking size consistency. */
zCVec zCMulMatVecNC(const zCMat m, const zCVec v1, zCVec v)
{
  int i, j;
  zComplex *e, z;

  e = zCMatBufNC(m);
  for( i=0; i<zCMatRowSizeNC(m); i++, e+=zCMatColSizeNC(m) ){
    zComplexZero( zCVecElemNC(v,i) );
    for( j=0; j<zCMatColSizeNC(m); j++ ){
      zComplexCMul( &e[j], zCVecElemNC(v1,j), &z );
      zComplexAdd( zCVecElemNC(v,i), &z, zCVecElemNC(v,i) );
    }
  }
  return v;
}

/* multiply a complex matrix and a complex column vector. */
zCVec zCMulMatVec(const zCMat m, const zCVec v1, zCVec v)
{
  if( !zCMatColCVecSizeEqual(m,v1) || !zCMatRowCVecSizeEqual(m,v) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  return zCMulMatVecNC( m, v1, v );
}

/* print a complex matrix out to a file. */
void zCMatFPrint(FILE *fp, const zCMat m)
{
  int i, j;

  if( !m )
    fprintf( fp, "(null matrix)\n" );
  else{
    fprintf( fp,"(%d, %d) {\n",
      zCMatRowSizeNC(m), zCMatColSizeNC(m) );
    for( i=0; i<zCMatRowSizeNC(m); i++ ){
      for( j=0; j<zCMatColSizeNC(m); j++ ){
        zComplexFPrint( fp, zCMatElemNC(m,i,j) );
        fprintf( fp, ", " );
      }
      fprintf( fp, "\n" );
    }
    fprintf( fp, "}\n" );
  }
}
