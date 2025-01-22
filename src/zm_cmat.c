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
  if( !zCMatSizeIsEqual( src, dest ) ){
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

/* multiply a complex matrix by a complex scalar without checking size consistency. */
zCMat zCMatMulNC(const zCMat m1, const zComplex *z, zCMat m)
{
  int i, n;

  n = zCMatRowSizeNC(m) * zCMatColSizeNC(m);
  for( i=0; i<n; i++ )
    zComplexCMul( &zCMatBufNC(m1)[i], z, &zCMatBufNC(m)[i] );
  return m;
}

/* divide a complex matrix by a complex scalar without checking size consistency. */
zCMat zCMatDivNC(const zCMat m1, const zComplex *z, zCMat m)
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
  if( !zCMatSizeIsEqual(m1,m2) || !zCMatSizeIsEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  return zCMatAddNC( m1, m2, m );
}

/* substract a complex matrix from another. */
zCMat zCMatSub(const zCMat m1, const zCMat m2, zCMat m)
{
  if( !zCMatSizeIsEqual(m1,m2) || !zCMatSizeIsEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  return zCMatSubNC( m1, m2, m );
}

/* reverse a complex matrix. */
zCMat zCMatRev(const zCMat m1, zCMat m)
{
  if( !zCMatSizeIsEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  return zCMatRevNC( m1, m );
}

/* multiply a complex matrix by a complex scalar. */
zCMat zCMatMul(const zCMat m1, const zComplex *z, zCMat m)
{
  if( !zCMatSizeIsEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  return zCMatMulNC( m1, z, m );
}

/* divide a complex matrix by a complex scalar. */
zCMat zCMatDiv(const zCMat m1, const zComplex *z, zCMat m)
{
  if( !zCMatSizeIsEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  if( zComplexIsTiny( z ) ){
    ZRUNERROR( ZM_ERR_ZERODIV );
    return NULL;
  }
  return zCMatDivNC( m1, z, m );
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
  if( zCMatColSizeNC(m) != zCVecSizeNC(v1) ||
      zCMatRowSizeNC(m) != zCVecSizeNC(v) ){
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
