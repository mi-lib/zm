/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_cmat - complex matrix class.
 */

#include <zm/zm_cmat.h>

/* ********************************************************** */
/* CLASS: zCMat
 * double precision floating point value matrix class
 * NOTES: each element of matrix(size=r*c) is at (0 - r-1,0 - c-1).
 * ********************************************************** */

/* zCMatAlloc
 * - allocate matrix.
 */
zCMat zCMatAlloc(int row, int col)
{
  zCMat m;

  if( !( m = zAlloc( zCMatStruct, 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  if( !( m->elem = zAlloc( zComplex, row*col ) ) ){
    ZALLOCERROR();
    zFree( m );
    return NULL;
  }
  zCMatSetSize( m, row, col );
  return m;
}

/* zCMatFree
 * - free matrix.
 */
void zCMatFree(zCMat m)
{
  if( m ){
    zFree( m->elem );
    free( m );
  }
}

/* zCMatClear
 * - cleanup matrix.
 */
zCMat zCMatClear(zCMat m)
{
  register int i, n;

  n = _zCMatRowSize(m) * _zCMatColSize(m);
  for( i=0; i<n; i++ )
    zComplexCopy( ZCOMPLEXZERO, &m->elem[i] );
  return m;
}

/* zCMatCopyNC
 * - copy matrix without checking size consistency.
 */
zCMat zCMatCopyNC(zCMat src, zCMat dest)
{
  register int i, n;

  n = _zCMatRowSize(src) * _zCMatColSize(src);
  for( i=0; i<n; i++ )
    zComplexCopy( &src->elem[i], &dest->elem[i] );
  return dest;
}

/* zCMatCopy
 * - copy matrix.
 */
zCMat zCMatCopy(zCMat src, zCMat dest)
{
  if( !zCMatSizeIsEqual( src, dest ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zCMatCopyNC( src, dest );
}

/* zCMatClone
 * - clone matrix.
 */
zCMat zCMatClone(zCMat src)
{
  zCMat dest;

  if( ( dest = zCMatAlloc( _zCMatRowSize(src), _zCMatColSize(src) ) ) )
    zCMatCopyNC( src, dest );
  return dest;
}

/* zMat2CMat
 * - convert real matrix to complex matrix.
 */
zCMat zMat2CMat(zMat m, zCMat cm)
{
  register int i, n;

  n = zMatRowSizeNC(m) * zMatColSizeNC(m);
  for( i=0; i<n; i++ )
    zComplexCreate( &cm->elem[i], zMatBuf(m)[i], 0 );
  return cm;
}

/* zCMatIsTol
 * - test for tiny matrix.
 */
bool zCMatIsTol(zCMat m, double tol)
{
  register int i, n;

  n = _zCMatRowSize(m) * _zCMatColSize(m);
  for( i=0; i<n; i++ )
    if( !zComplexIsTol( &m->elem[i], tol ) ) return false;
  return true;
}

/* zCMatAddNC
 * - add matrix without checking size consistency.
 */
zCMat zCMatAddNC(zCMat m1, zCMat m2, zCMat m)
{
  register int i, n;

  n = _zCMatRowSize(m) * _zCMatColSize(m);
  for( i=0; i<n; i++ )
    zComplexAdd( &m1->elem[i], &m2->elem[i], &m->elem[i] );
  return m;
}

/* zCMatSubNC
 * - subtract matrix without checking size consistency.
 */
zCMat zCMatSubNC(zCMat m1, zCMat m2, zCMat m)
{
  register int i, n;

  n = _zCMatRowSize(m) * _zCMatColSize(m);
  for( i=0; i<n; i++ )
    zComplexSub( &m1->elem[i], &m2->elem[i], &m->elem[i] );
  return m;
}

/* zCMatRevNC
 * - reverse matrix without checking size consistency.
 */
zCMat zCMatRevNC(zCMat m1, zCMat m)
{
  register int i, n;

  n = _zCMatRowSize(m) * _zCMatColSize(m);
  for( i=0; i<n; i++ )
    zComplexRev( &m1->elem[i], &m->elem[i] );
  return m;
}

/* zCMatMulNC
 * - multiply matrix without checking size consistency.
 */
zCMat zCMatMulNC(zCMat m1, zComplex *z, zCMat m)
{
  register int i, n;

  n = _zCMatRowSize(m) * _zCMatColSize(m);
  for( i=0; i<n; i++ )
    zComplexCMul( &m1->elem[i], z, &m->elem[i] );
  return m;
}

/* zCMatDivNC
 * - divide matrix without checking size consistency.
 */
zCMat zCMatDivNC(zCMat m1, zComplex *z, zCMat m)
{
  register int i, n;
  double r;
  zComplex dz;

  r = zComplexSqrAbs( z );
  zComplexConj( z, &dz );
  zComplexDiv( &dz, r, &dz );
  n = _zCMatRowSize(m) * _zCMatColSize(m);
  for( i=0; i<n; i++ )
    zComplexCMul( &m1->elem[i], &dz, &m->elem[i] );
  return m;
}

/* zCMatAdd
 * - add matrix.
 */
zCMat zCMatAdd(zCMat m1, zCMat m2, zCMat m)
{
  if( !zCMatSizeIsEqual(m1,m2) || !zCMatSizeIsEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zCMatAddNC( m1, m2, m );
}

/* zCMatSub
 * - substraction of matrix.
 */
zCMat zCMatSub(zCMat m1, zCMat m2, zCMat m)
{
  if( !zCMatSizeIsEqual(m1,m2) || !zCMatSizeIsEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zCMatSubNC( m1, m2, m );
}

/* zCMatRev
 * - reverse matrix.
 */
zCMat zCMatRev(zCMat m1, zCMat m)
{
  if( !zCMatSizeIsEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zCMatRevNC( m1, m );
}

/* zCMatMul
 * - multiply matrix.
 */
zCMat zCMatMul(zCMat m1, zComplex *z, zCMat m)
{
  if( !zCMatSizeIsEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zCMatMulNC( m1, z, m );
}

/* zCMatDiv
 * - divide matrix.
 */
zCMat zCMatDiv(zCMat m1, zComplex *z, zCMat m)
{
  if( !zCMatSizeIsEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  if( zComplexIsTiny( z ) ){
    ZRUNERROR( ZM_ERR_ZERODIV );
    return NULL;
  }
  return zCMatDivNC( m1, z, m );
}

/* zCMulMatVecNC
 * - multiply a matrix and a column vector without
 *   checking size consistency.
 */
zCVec zCMulMatVecNC(zCMat m, zCVec v1, zCVec v)
{
  register int i, j;
  zComplex *e, z;

  e = m->elem;
  for( i=0; i<_zCMatRowSize(m); i++, e+=_zCMatColSize(m) ){
    zComplexClear( zCVecElem(v,i) );
    for( j=0; j<_zCMatColSize(m); j++ ){
      zComplexCMul( &e[j], zCVecElem(v1,j), &z );
      zComplexAdd( zCVecElem(v,i), &z, zCVecElem(v,i) );
    }
  }
  return v;
}

/* zCMulMatVec
 * - multiply a matrix and a column vector.
 */
zCVec zCMulMatVec(zCMat m, zCVec v1, zCVec v)
{
  if( zCMatColSize(m) != zCVecSize(v1) ||
      zCMatRowSize(m) != zCVecSize(v) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  return zCMulMatVecNC( m, v1, v );
}

/* zCMatFWrite
 * - output of matrix to file.
 */
void zCMatFWrite(FILE *fp, zCMat m)
{
  register int i, j;

  if( !m )
    fprintf( fp, "(null matrix)\n" );
  else{
    fprintf( fp,"(%d, %d) {\n",
      _zCMatRowSize(m), _zCMatColSize(m) );
    for( i=0; i<_zCMatRowSize(m); i++ ){
      for( j=0; j<_zCMatColSize(m); j++ ){
        zComplexFWrite( fp, zCMatElem(m,i,j) );
        fprintf( fp, ", " );
      }
      fprintf( fp, "\n" );
    }
    fprintf( fp, "}\n" );
  }
}
