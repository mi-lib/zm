/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_vec_mat - vector and matrix class : matrix class.
 */

#include <zm/zm_mat.h>

/* ********************************************************** */
/* CLASS: zMat
 * double precision floating point value matrix class
 * NOTES: each element of matrix(size=r*c) is at (0 - r-1,0 - c-1).
 * ********************************************************** */

static zMat _zMatSetElemList(zMat m, va_list args);

/* (static)
 * _zMatSetElemList
 * - set matrix components from value list.
 *
 */
zMat _zMatSetElemList(zMat m, va_list args)
{
  register int i, j;

  for( i=0; i<zMatRowSizeNC(m); i++ )
    for( j=0; j<zMatColSizeNC(m); j++ )
      zMatSetElem( m, i, j, (double)va_arg( args, double ) );
  return m;
}

/* zMatSetElemList
 * - set matrix components from value list.
 */
zMat zMatSetElemList(zMat m, ... )
{
  va_list args;

  va_start( args, m );
  _zMatSetElemList( m, args );
  va_end( args );
  return m;
}

/* zMatAlloc
 * - allocate matrix.
 */
zMat zMatAlloc(int row, int col)
{
  zMat m;

  if( !( m = zAlloc( zMatStruct, 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  if( !( zMatBuf(m) = zAlloc( double, row*col ) ) ){
    ZALLOCERROR();
    zFree( m );
    return NULL;
  }
  zMatSetSize( m, row, col );
  return m;
}

/* zMatCreateList
 * create a matrix from value list.
 */
zMat zMatCreateList(int row, int col, ... )
{
  zMat m;
  va_list args;

  if( !( m = zMatAlloc( row, col ) ) ) return NULL;
  va_start( args, col );
  _zMatSetElemList( m, args );
  va_end( args );
  return m;
}

/* zMatFree
 * - free matrix.
 */
void zMatFree(zMat m)
{
  if( m ){
    zFree( zMatBuf(m) );
    free( m );
  }
}

/* zMatFreeAO
 * - free matrices at once.
 */
void zMatFreeAO(int n, ...)
{
  va_list arg;
  zMat m;
  register int i;

  va_start( arg, n );
  for( i=0; i<n; i++ ){
    m = va_arg( arg, zMat );
    zMatFree( m );
  }
  va_end( arg );
}

/* zMatClear
 * - cleanup matrix.
 */
zMat zMatClear(zMat m)
{
  zRawMatClear( zMatBuf(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
  return m;
}

/* zMatTouchup
 * - touchup matrix.
 */
zMat zMatTouchup(zMat m)
{
  zRawMatTouchup( zMatBuf(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
  return m;
}

/* zMatIdentNC
 * - create identity matrix without checking size consistency.
 */
zMat zMatIdentNC(zMat m)
{
  zRawMatIdent( zMatBuf(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
  return m;
}

/* zMatDiagNC
 * - create diagonal matrix without checking size consistency.
 */
zMat zMatDiagNC(zMat m, zVec d)
{
  zRawMatDiag( zMatBuf(m), zMatRowSizeNC(m), zMatColSizeNC(m), zVecBuf(d) );
  return m;
}

/* zMatIdent
 * - create identity matrix.
 */
zMat zMatIdent(zMat m)
{
  if( !zMatIsSqr( m ) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return NULL;
  }
  return zMatIdentNC( m );
}

/* zMatDiag
 * - create diagonal matrix.
 */
zMat zMatDiag(zMat m, zVec d)
{
  if( !zMatIsSqr( m ) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return NULL;
  }
  if( !zMatColVecSizeIsEqual( m, d ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  return zMatDiagNC( m, d );
}

/* zMatRandUniform
 * - create a random matrix with a uniform range.
 */
zMat zMatRandUniform(zMat m, double min, double max)
{
  zRawMatRandUniform( zMatBuf(m), zMatRowSizeNC(m), zMatColSizeNC(m), min, max );
  return m;
}

/* zMatRand
 * - create a random matrix with range matrices.
 */
zMat zMatRand(zMat m, zMat min, zMat max)
{
  zRawMatRand( zMatBuf(m), zMatBuf(min), zMatBuf(max), zMatRowSizeNC(m), zMatColSizeNC(m) );
  return m;
}

/* zMatCopyNC
 * - copy matrix without checking size consistency.
 */
zMat zMatCopyNC(zMat src, zMat dest)
{
  zRawMatCopy( zMatBuf(src), zMatBuf(dest),
    zMatRowSizeNC(dest), zMatColSizeNC(dest) );
  return dest;
}

/* zMatCopy
 * - copy matrix.
 */
zMat zMatCopy(zMat src, zMat dest)
{
  if( !zMatSizeIsEqual( src, dest ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMatCopyNC( src, dest );
}

/* zMatCopyArray
 * - copy matrix from 2-dim array of double precision
 *   floating point values.
 */
zMat zMatCopyArray(double array[], int r, int c, zMat m)
{
  if( zMatRowSize(m)!=r || zMatColSize(m)!=c ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  zRawMatCopy( array, zMatBuf(m), r, c );
  return m;
}

/* zMatClone
 * - clone matrix.
 */
zMat zMatClone(zMat src)
{
  zMat dest;

  if( ( dest = zMatAlloc( zMatRowSizeNC(src), zMatColSizeNC(src) ) ) )
    zMatCopyNC( src, dest );
  return dest;
}

/* zMatCloneArray
 * - create matrix from an array of double precision
 *   floating point values.
 */
zMat zMatCloneArray(double array[], int r, int c)
{
  zMat m;

  if( ( m = zMatAlloc( r, c ) ) )
    zMatCopyArray( array, r, c, m );
  return m;
}

/* zMatGetNC
 * - get submatrix without checking the size validity.
 */
zMat zMatGetNC(zMat src, int pr, int pc, zMat dest)
{
  zRawMatGet( zMatBuf(src), zMatRowSizeNC(src), zMatColSizeNC(src), pr, pc, zMatBuf(dest), zMatRowSizeNC(dest), zMatColSizeNC(dest) );
  return dest;
}

/* zMatGet
 * - get submatrix.
 */
zMat zMatGet(zMat src, int pr, int pc, zMat dest)
{
  if( pr+zMatRowSize(dest) > zMatRowSize(src) ||
      pc+zMatColSize(dest) > zMatColSize(src) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMatGetNC( src, pr, pc, dest );
}

/* zMatTGetNC
 * - get transpose of a submatrix without checking the size validity.
 */
zMat zMatTGetNC(zMat src, int pr, int pc, zMat dest)
{
  zRawMatTGet( zMatBuf(src), zMatRowSizeNC(src), zMatColSizeNC(src), pr, pc, zMatBuf(dest), zMatRowSizeNC(dest), zMatColSizeNC(dest) );
  return dest;
}

/* zMatTGet
 * - get transpose of a submatrix.
 */
zMat zMatTGet(zMat src, int pr, int pc, zMat dest)
{
  if( pr+zMatColSize(dest) > zMatRowSize(src) ||
      pc+zMatRowSize(dest) > zMatColSize(src) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMatTGetNC( src, pr, pc, dest );
}

/* zMatPutNC
 * - put submatrix without checking the size validity.
 */
zMat zMatPutNC(zMat dest, int pr, int pc, zMat src)
{
  zRawMatPut( zMatBuf(dest), zMatRowSizeNC(dest), zMatColSizeNC(dest), pr, pc, zMatBuf(src), zMatRowSizeNC(src), zMatColSizeNC(src) );
  return dest;
}

/* zMatPut
 * - put submatrix.
 */
zMat zMatPut(zMat dest, int pr, int pc, zMat src)
{
  if( pr+zMatRowSize(src) > zMatRowSize(dest) ||
      pc+zMatColSize(src) > zMatColSize(dest) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMatPutNC( dest, pr, pc, src );
}

/* zMatTPutNC
 * - put transpose of a submatrix without checking the size validity.
 */
zMat zMatTPutNC(zMat dest, int pr, int pc, zMat src)
{
  zRawMatTPut( zMatBuf(dest), zMatRowSizeNC(dest), zMatColSizeNC(dest), pr, pc, zMatBuf(src), zMatRowSizeNC(src), zMatColSizeNC(src) );
  return dest;
}

/* zMatTPut
 * - put transpose of a submatrix.
 */
zMat zMatTPut(zMat dest, int pr, int pc, zMat src)
{
  if( pr+zMatColSize(src) > zMatRowSize(dest) ||
      pc+zMatRowSize(src) > zMatColSize(dest) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMatTPutNC( dest, pr, pc, src );
}

/* zMatGetRowNC
 * - abstract row vector of matrix without checking
 *   size consistency.
 */
zVec zMatGetRowNC(zMat m, int row, zVec v)
{
  zRawMatGetRow( zMatBuf(m), zMatRowSizeNC(m), zMatColSizeNC(m),
    row, zVecBuf(v) );
  return v;
}

/* zMatGetColNC
 * - abstract column vector of matrix without checking
 *   size consistency.
 */
zVec zMatGetColNC(zMat m, int col, zVec v)
{
  zRawMatGetCol( zMatBuf(m), zMatRowSizeNC(m), zMatColSizeNC(m),
    col, zVecBuf(v) );
  return v;
}

/* zMatGetRow
 * - abstract row vector of matrix.
 */
zVec zMatGetRow(zMat m, int row, zVec v)
{
  if( !zMatColVecSizeIsEqual( m, v ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  if( row >= zMatRowSize(m) ){
    ZRUNERROR( ZM_ERR_INV_ROW );
    return NULL;
  }
  return zMatGetRowNC( m, row, v );
}

/* zMatGetCol
 * - abstract column vector of matrix.
 */
zVec zMatGetCol(zMat m, int col, zVec v)
{
  if( !zMatRowVecSizeIsEqual( m, v ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  if( col >= zMatColSize(m) ){
    ZRUNERROR( ZM_ERR_INV_COL );
    return NULL;
  }
  return zMatGetColNC( m, col, v );
}

/* zMatSetRowNC
 * - set of row vector to matrix without checking size consistency.
 */
zMat zMatSetRowNC(zMat m, int row, zVec v)
{
  zRawMatSetRow( zMatBuf(m), zMatRowSizeNC(m), zMatColSizeNC(m),
    row, zVecBuf(v) );
  return m;
}

/* zMatColNC
 * - set of column vector to matrix without checking size consistency.
 */
zMat zMatSetColNC(zMat m, int col, zVec v)
{
  zRawMatSetCol( zMatBuf(m), zMatRowSizeNC(m), zMatColSizeNC(m),
    col, zVecBuf(v) );
  return m;
}

/* zMatSetRow
 * - set of row vector to matrix.
 */
zMat zMatSetRow(zMat m, int row, zVec v)
{
  if( !zMatColVecSizeIsEqual( m, v ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  if( row >= zMatRowSize(m) ){
    ZRUNERROR( ZM_ERR_INV_ROW );
    return NULL;
  }
  return zMatSetRowNC( m, row, v );
}

/* zMatSetCol
 * - set of column vector to matrix.
 */
zMat zMatSetCol(zMat m, int col, zVec v)
{
  if( !zMatRowVecSizeIsEqual( m, v ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  if( col >= zMatColSize(m) ){
    ZRUNERROR( ZM_ERR_INV_COL );
    return NULL;
  }
  return zMatSetColNC( m, col, v );
}

/* zMatSwapRowNC
 * - swap of matrix row without checking size consistency.
 */
zMat zMatSwapRowNC(zMat m, int r1, int r2)
{
  zRawMatSwapRow( zMatBuf(m), zMatRowSizeNC(m), zMatColSizeNC(m), r1, r2 );
  return m;
}

/* zMatSwapRow
 * - swap of matrix row.
 */
zMat zMatSwapRow(zMat m, int r1, int r2)
{
  if( r1 >= zMatRowSizeNC(m) || r2 >= zMatRowSizeNC(m) ){
    ZRUNERROR( ZM_ERR_INV_ROW );
    return NULL;
  }
  return zMatSwapRowNC( m, r1, r2 );
}

/* zMatSwapColNC
 * - swap of matrix column without checking size consistency.
 */
zMat zMatSwapColNC(zMat m, int c1, int c2)
{
  zRawMatSwapCol( zMatBuf(m), zMatRowSizeNC(m), zMatColSizeNC(m), c1, c2 );
  return m;
}

/* zMatSwapColNC
 * - swap of matrix column.
 */
zMat zMatSwapCol(zMat m, int c1, int c2)
{
  if( c1 >= zMatColSizeNC(m) || c2 >= zMatColSizeNC(m) ){
    ZRUNERROR( ZM_ERR_INV_COL );
    return NULL;
  }
  return zMatSwapColNC( m, c1, c2 );
}

/* zMatShift
 * - shift diagonal values of a matrix.
 */
void zMatShift(zMat m, double shift)
{
  register int i;

  for( i=0; i<zMatRowSizeNC(m); i++ )
    zMatElem( m, i, i ) += shift;
}

/* zMatIsEqual
 * - see if two matrices are equal.
 */
bool zMatIsEqual(zMat m1, zMat m2)
{
  register int i, j;

  if( !zMatSizeIsEqual( m1, m2 ) ) return false;
  for( i=0; i<zMatRowSizeNC(m1); i++ )
    for( j=0; j<zMatColSizeNC(m1); j++ )
      if( !zIsTiny( zMatElem(m1,i,j) - zMatElem(m2,i,j) ) ) return false;
  return true;
}

/* zMatIsTol
 * - test for tiny matrix.
 */
bool zMatIsTol(zMat m, double tol)
{
  return zRawMatIsTol( zMatBuf(m), zMatRowSizeNC(m), zMatColSizeNC(m), tol );
}

/* zMatRowReg
 * - row regression of matrix.
 */
zMat zMatRowReg(zMat m, int rank)
{
  if( rank < zMatRowSizeNC(m) ) zMatSetRowSize( m, rank );
  return m;
}

/* zMatColReg
 * - column regression of matrix.
 */
zMat zMatColReg(zMat m, int rank)
{
  register int i;
  double *sp, *dp;

  if( rank < zMatColSizeNC(m) ){
    for( sp=dp=zMatBuf(m), i=0; i<zMatRowSizeNC(m); i++, sp+=zMatColSizeNC(m), dp+=rank )
      memmove( dp, sp, sizeof(double)*rank );
    zMatSetColSize( m, rank );
  }
  return m;
}

/* zMatAddNC
 * - add matrix without checking size consistency.
 */
zMat zMatAddNC(zMat m1, zMat m2, zMat m)
{
  zRawMatAdd( zMatBuf(m1), zMatBuf(m2),
    zMatBuf(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
  return m;
}

/* zMatSubNC
 * - subtract matrix without checking size consistency.
 */
zMat zMatSubNC(zMat m1, zMat m2, zMat m)
{
  zRawMatSub( zMatBuf(m1), zMatBuf(m2),
    zMatBuf(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
  return m;
}

/* zMatRevNC
 * - reverse matrix without checking size consistency.
 */
zMat zMatRevNC(zMat m1, zMat m)
{
  zRawMatRev( zMatBuf(m1), zMatBuf(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
  return m;
}

/* zMatMulNC
 * - multiply matrix without checking size consistency.
 */
zMat zMatMulNC(zMat m1, double k, zMat m)
{
  zRawMatMul( zMatBuf(m1), k, zMatBuf(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
  return m;
}

/* zMatDivNC
 * - divide matrix without checking size consistency.
 */
zMat zMatDivNC(zMat m1, double k, zMat m)
{
  zRawMatDiv( zMatBuf(m1), k, zMatBuf(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
  return m;
}

/* zMatCatNC
 * - concatenate matrix without checking size consistency.
 */
zMat zMatCatNC(zMat m1, double k, zMat m2, zMat m)
{
  zRawMatCat( zMatBuf(m1), k, zMatBuf(m2), zMatBuf(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
  return m;
}

/* zMatAdd
 * - add matrix.
 */
zMat zMatAdd(zMat m1, zMat m2, zMat m)
{
  if( !zMatSizeIsEqual(m1,m2) || !zMatSizeIsEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMatAddNC( m1, m2, m );
}

/* zMatSub
 * - substraction of matrix.
 */
zMat zMatSub(zMat m1, zMat m2, zMat m)
{
  if( !zMatSizeIsEqual(m1,m2) || !zMatSizeIsEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMatSubNC( m1, m2, m );
}

/* zMatRev
 * - reverse matrix.
 */
zMat zMatRev(zMat m1, zMat m)
{
  if( !zMatSizeIsEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMatRevNC( m1, m );
}

/* zMatMul
 * - multiply matrix.
 */
zMat zMatMul(zMat m1, double k, zMat m)
{
  if( !zMatSizeIsEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMatMulNC( m1, k, m );
}

/* zMatDiv
 * - divide matrix.
 */
zMat zMatDiv(zMat m1, double k, zMat m)
{
  if( !zMatSizeIsEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  if( zIsTiny( k ) ){
    ZRUNERROR( ZM_ERR_ZERODIV );
    return NULL;
  }
  return zMatDivNC( m1, k, m );
}

/* zMatCat
 * - concatenate matrix.
 */
zMat zMatCat(zMat m1, double k, zMat m2, zMat m)
{
  if( !zMatSizeIsEqual(m1,m2) || !zMatSizeIsEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMatCatNC( m1, k, m2, m );
}

/* zMatSqrNorm
 * - squared norm of matrix.
 */
double zMatSqrNorm(zMat m)
{
  return zRawMatSqrNorm( zMatBuf(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
}

/* zMatInfNorm
 * - infinity norm of matrix.
 */
double zMatInfNorm(zMat m)
{
  double *mp, rs, rsmax = 0;
  register int i, j;

  mp = zMatBuf(m);
  for( i=0; i<zMatRowSizeNC(m); i++ ){
    for( rs=0, j=0; j<zMatColSizeNC(m); j++ )
      rs += fabs( *mp++ );
    if( rs > rsmax ) rsmax = rs;
  }
  return rsmax;
}

/* zMatTNC
 * - transpose of matrix without checking size consistency.
 */
zMat zMatTNC(zMat m, zMat tm)
{
  zRawMatT( zMatBuf(m), zMatBuf(tm),
    zMatRowSizeNC(tm), zMatColSizeNC(tm) );
  return tm;
}

/* zMatT
 * - transpose of matrix.
 */
zMat zMatT(zMat m, zMat tm)
{
  if( !zMatColRowSizeIsEqual( tm, m ) ||
      !zMatColRowSizeIsEqual( m, tm ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMatTNC( m, tm );
}

/* zMatTDST
 * - transpose of matrix (destructive).
 */
zMat zMatTDST(zMat m)
{
  int row, col;

  row = zMatRowSizeNC(m);
  col = zMatColSizeNC(m);
  zRawMatTDST( zMatBuf(m), row, col );
  zMatSetSize( m, col, row );
  return m;
}

/* zMatTClone
 * - clone transpose of a matrix.
 */
zMat zMatTClone(zMat src)
{
  zMat dest;

  if( ( dest = zMatAlloc( zMatColSizeNC(src), zMatRowSizeNC(src) ) ) )
    zMatTNC( src, dest );
  return dest;
}

/* zVecDyadNC
 * - dyadic product of vectors without checking size
 *   consistency.
 */
zMat zVecDyadNC(zVec v1, zVec v2, zMat dyad)
{
  zRawVecDyad( zVecBuf(v1), zVecSizeNC(v1), zVecBuf(v2), zVecSizeNC(v2), zMatBuf(dyad) );
  return dyad;
}

/* zVecDyad
 * - dyadic product of vectors.
 */
zMat zVecDyad(zVec v1, zVec v2, zMat dyad)
{
  if( zVecSize(v1) != zMatRowSize(dyad) ||
      zVecSize(v2) != zMatColSize(dyad) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  return zVecDyadNC( v1, v2, dyad );
}

/* zMatAddDyadNC
 * - add dyadic matrix of vectors to matrix without checking
 *   size consistency.
 */
zMat zMatAddDyadNC(zMat m, zVec v1, zVec v2)
{
  zRawMatAddDyad( zMatBuf(m), zVecBuf(v1), zVecSizeNC(v1), zVecBuf(v2), zVecSizeNC(v2) );
  return m;
}

/* zMatAddDyad
 * - add dyadic matrix of vectors to matrix.
 */
zMat zMatAddDyad(zMat m, zVec v1, zVec v2)
{
  if( zVecSize(v1) != zMatRowSize(m) ||
      zVecSize(v2) != zMatColSize(m) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  return zMatAddDyadNC( m, v1, v2 );
}

/* zMatSubDyadNC
 * - subtract dyadic matrix of vectors to matrix without checking
 *   size consistency.
 */
zMat zMatSubDyadNC(zMat m, zVec v1, zVec v2)
{
  zRawMatSubDyad( zMatBuf(m), zVecBuf(v1), zVecSizeNC(v1), zVecBuf(v2), zVecSizeNC(v2) );
  return m;
}

/* zMatSubDyad
 * - subtract dyadic matrix of vectors to matrix.
 */
zMat zMatSubDyad(zMat m, zVec v1, zVec v2)
{
  if( zVecSize(v1) != zMatRowSize(m) ||
      zVecSize(v2) != zMatColSize(m) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  return zMatSubDyadNC( m, v1, v2 );
}

/* zMatCatDyadNC
 * - add multiplied dyadic matrix of vectors by scalor to matrix
 *   without checking size consistency.
 */
zMat zMatCatDyadNC(zMat m, double k, zVec v1, zVec v2)
{
  zRawMatCatDyad( zMatBuf(m), k, zVecBuf(v1), zVecSizeNC(v1), zVecBuf(v2), zVecSizeNC(v2) );
  return m;
}

/* zMatCatDyad
 * - add multiplied dyadic matrix of vectors by scalor to matrix.
 */
zMat zMatCatDyad(zMat m, double k, zVec v1, zVec v2)
{
  if( zVecSize(v1) != zMatRowSize(m) ||
      zVecSize(v2) != zMatColSize(m) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  return zMatCatDyadNC( m, k, v1, v2 );
}

/* zMatTrNC
 * - calculation of trace value of matrix without checking
 *   size consistency.
 */
double zMatTrNC(zMat m)
{
  return zRawMatTr( zMatBuf(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
}

/* zMatTr
 * - calculation of trace value of matrix.
 */
double zMatTr(zMat m)
{
  if( !zMatIsSqr( m ) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return 0;
  }
  return zMatTrNC( m );
}

/* zMulMatVecNC
 * - multiply a matrix and a column vector without
 *   checking size consistency.
 */
zVec zMulMatVecNC(zMat m, zVec v1, zVec v)
{
  zRawMulMatVec( zMatBuf(m), zVecBuf(v1),
    zMatRowSizeNC(m), zMatColSizeNC(m), zVecBuf(v) );
  return v;
}

/* zMulVecMatNC
 * - multiply a row vector and a matrix without
 *   checking size consistency.
 */
zVec zMulVecMatNC(zVec v1, zMat m, zVec v)
{
  zRawMulVecMat( zVecBuf(v1), zMatBuf(m),
    zMatRowSizeNC(m), zMatColSizeNC(m), zVecBuf(v) );
  return v;
}

/* zMulMatMatNC
 * - multiply two matrices without checking
 *   size consistency.
 */
zMat zMulMatMatNC(zMat m1, zMat m2, zMat m)
{
  zRawMulMatMat( zMatBuf(m1),zMatRowSizeNC(m1), zMatColSizeNC(m1),
    zMatBuf(m2), zMatRowSizeNC(m2), zMatColSizeNC(m2), zMatBuf(m) );
  return m;
}

/* zMulMatMatTNC
 * - multiply a matrix and transpose of a matrix
 *   ('m = m1 m2^T') without checking size consistency.
 */
zMat zMulMatMatTNC(zMat m1, zMat m2, zMat m)
{
  zRawMulMatMatT( zMatBuf(m1), zMatRowSizeNC(m1), zMatColSizeNC(m1),
    zMatBuf(m2), zMatRowSizeNC(m2), zMatColSizeNC(m2), zMatBuf(m) );
  return m;
}

/* zMulMatTMatNC
 * - multiply transpose of a matrix and a matrix
 *   'm = m1^T m2' without checking size consistency.
 */
zMat zMulMatTMatNC(zMat m1, zMat m2, zMat m)
{
  zRawMulMatTMat( zMatBuf(m1), zMatRowSizeNC(m1), zMatColSizeNC(m1),
    zMatBuf(m2), zMatRowSizeNC(m2), zMatColSizeNC(m2), zMatBuf(m) );
  return m;
}

/* zMulMatVec
 * - multiply a matrix and a column vector.
 */
zVec zMulMatVec(zMat m, zVec v1, zVec v)
{
  if( !zMatColVecSizeIsEqual( m, v1 ) || !zMatRowVecSizeIsEqual( m, v ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  return zMulMatVecNC( m, v1, v );
}

/* zMulVecMat
 * - multiply a row vector and a matrix.
 */
zVec zMulVecMat(zVec v1, zMat m, zVec v)
{
  if( !zMatRowVecSizeIsEqual( m, v1 ) || !zMatColVecSizeIsEqual( m, v ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  return zMulVecMatNC( v1, m, v );
}

/* zMulMatMat
 * - multiply two matrices.
 */
zMat zMulMatMat(zMat m1, zMat m2, zMat m)
{
  if( !zMatColRowSizeIsEqual( m1, m2 ) ||
      !zMatRowSizeIsEqual( m1, m ) || !zMatColSizeIsEqual( m2, m ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMulMatMatNC( m1, m2, m );
}

/* zMulMatMatT
 * - multiply a matrix and transpose of a matrix
 *   ('m = m1 m2^T').
 */
zMat zMulMatMatT(zMat m1, zMat m2, zMat m)
{
  if( !zMatColSizeIsEqual( m1, m2 ) || !zMatRowSizeIsEqual( m1, m ) ||
      !zMatRowColSizeIsEqual( m2, m ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMulMatMatTNC( m1, m2, m );
}

/* zMulMatTMat
 * - multiply transpose of a matrix and a matrix
 *   ('m = m1^T m2').
 */
zMat zMulMatTMat(zMat m1, zMat m2, zMat m)
{
  if( !zMatRowSizeIsEqual( m1, m2 ) ||
      !zMatColRowSizeIsEqual( m1, m ) ||
      !zMatColSizeIsEqual( m2, m ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMulMatTMatNC( m1, m2, m );
}

/* zMulMatVecDRC
 * - multiply a matrix and a column vector.
 *   the result is directly put into 'v'.
 */
zVec zMulMatVecDRC(zMat m, zVec v)
{
  zVec tmp;

  if( !( tmp = zVecAlloc( zVecSizeNC(v) ) ) ) return NULL;
  zMulMatVec( m, v, tmp );
  zVecCopyNC( tmp, v );
  zVecFree( tmp );
  return v;
}

/* zMulVecMatDRC
 * - multiply a row vector and a matrix.
 *   the result is directly put into 'v'.
 */
zVec zMulVecMatDRC(zVec v, zMat m)
{
  zVec tmp;

  if( !( tmp = zVecAlloc( zVecSizeNC(v) ) ) ) return NULL;
  zMulVecMat( v, m, tmp );
  zVecCopyNC( tmp, v );
  zVecFree( tmp );
  return v;
}

/* zMulMatMatDRC
 * - multiply two matrices.
 *   the result is directly put into 'm2'.
 *   note that 'm1' must be a square matrix.
 */
zMat zMulMatMatDRC(zMat m1, zMat m2)
{
  zMat tmp;

  if( !( tmp = zMatAlloc(zMatRowSizeNC(m2),zMatColSizeNC(m2)) ) )
    return NULL;
  zMulMatMat( m1, m2, tmp );
  zMatCopyNC( tmp, m2 );
  zMatFree( tmp );
  return m2;
}

/* zMulMatMatTDRC
 * - multiply a matrix and transpose of a matrix
 *   ('m = m1 m2^T').
 *   the result is directly put into 'm2'.
 */
zMat zMulMatMatTDRC(zMat m1, zMat m2)
{
  zMat tmp;

  if( !( tmp = zMatAlloc(zMatRowSizeNC(m2),zMatColSizeNC(m2)) ) )
    return NULL;
  if( zMulMatMatT( m1, m2, tmp ) )
    zMatCopyNC( tmp, m2 );
  zMatFree( tmp );
  return m2;
}

/* zMulMatTMatDRC
 * - multiply transpose of a matrix and a matrix
 *   ('m = m1^T m2').
 *   the result is directly put into 'm2'.
 */
zMat zMulMatTMatDRC(zMat m1, zMat m2)
{
  zMat tmp;

  if( !( tmp = zMatAlloc(zMatRowSizeNC(m2),zMatColSizeNC(m2)) ) )
    return NULL;
  if( zMulMatTMat( m1, m2, tmp ) )
    zMatCopyNC( tmp, m2 );
  zMatFree( tmp );
  return m2;
}

/* zMatQuadNC
 * - quadratic multiplication of matrices ('q = a diag{w} a^T')
 *   without checking size consistency.
 */
zMat zMatQuadNC(zMat a, zVec w, zMat q)
{
  register int i, j, k;
  double wa;

  zMatClear( q );
  for( k=0; k<zMatColSizeNC(a); k++ )
    for( i=0; i<zMatRowSizeNC(a); i++ ){
      wa = w ? zVecElem(w,k) * zMatElem(a,i,k) : zMatElem(a,i,k);
      for( j=i; j<zMatRowSizeNC(a); j++ )
        zMatElem(q,i,j) += wa * zMatElem(a,j,k);
    }
  for( i=0; i<zMatRowSizeNC(a); i++ )
    for( j=i; j<zMatRowSizeNC(a); j++ )
      zMatSetElem( q, j, i, zMatElem(q,i,j) );
  return q;
}

/* zMatQuad
 * - quadratic multiplication of matrices ('q = a diag{w} a^T').
 */
zMat zMatQuad(zMat a, zVec w, zMat q)
{
  if( w && !zMatColVecSizeIsEqual( a, w ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  if( !zMatRowSizeIsEqual( a, q ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  if( !zMatIsSqr( q ) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return NULL;
  }
  return zMatQuadNC( a, w, q );
}

/* zMatTQuadNC
 * - quadratic multiplication of matrices ('q = a^T diag{w} a')
 *   without checking size consistency.
 */
zMat zMatTQuadNC(zMat a, zVec w, zMat q)
{
  register int i, j, k;
  double wa;

  zMatClear( q );
  for( k=0; k<zMatRowSizeNC(a); k++ )
    for( i=0; i<zMatColSizeNC(a); i++ ){
      wa = w ? zVecElem(w,k) * zMatElem(a,k,i) : zMatElem(a,k,i);
      for( j=i; j<zMatColSizeNC(a); j++ )
        zMatElem(q,i,j) += wa * zMatElem(a,k,j);
    }
  for( i=0; i<zMatColSizeNC(a); i++ )
    for( j=i; j<zMatColSizeNC(a); j++ )
      zMatSetElem( q, j, i, zMatElem(q,i,j) );
  return q;
}

/* zMatTQuad
 * - quadratic multiplication of matrices ('q = a^T diag{w} a').
 */
zMat zMatTQuad(zMat a, zVec w, zMat q)
{
  if( w && !zMatRowVecSizeIsEqual( a, w ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  if( !zMatColRowSizeIsEqual( a, q ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  if( !zMatIsSqr( q ) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return NULL;
  }
  return zMatTQuadNC( a, w, q );
}

/* zMatReadFile
 * - input of matrix from file.
 */
zMat zMatReadFile(char filename[])
{
  FILE *fp;
  zMat m;

  if( !( fp = zOpenFile( filename, ZMATRIX_SUFFIX, "r" ) ) )
    return NULL;
  m = zMatFRead( fp );
  fclose( fp );
  return m;
}

/* zMatFRead
 * - input of matrix from file.
 */
zMat zMatFRead(FILE *fp)
{
  register unsigned i, j, row, col;
  zMat m;

  row = zFInt( fp );
  col = zFInt( fp );
  if( !( m = zMatAlloc( row, col ) ) ) return NULL;

  for( i=0; i<row; i++ )
    for( j=0; j<col; j++ )
      zMatSetElem( m, i, j, zFDouble( fp ) );
  return m;
}

/* zMatFWrite
 * - output of matrix to file.
 */
void zMatFWrite(FILE *fp, zMat m)
{
  register int i, j;

  if( !m )
    fprintf( fp, "(null matrix)\n" );
  else{
    fprintf( fp,"(%d, %d) {\n",
      zMatRowSizeNC(m), zMatColSizeNC(m) );
    for( i=0; i<zMatRowSizeNC(m); i++ ){
      for( j=0; j<zMatColSizeNC(m); j++ )
        fprintf( fp, " %.10g", zMatElem( m, i, j ) );
      fprintf( fp, "\n" );
    }
    fprintf( fp, "}\n" );
  }
}

/* zMatImg
 * - visualize matrix using one-charactor collage for debug.
 */
void zMatImg(zMat m)
{
  double max, min, d;
  const char pat_pos[] = ",x*M";
  const char pat_neg[] = ".oO@";
  register int i, j, c;

  max = fabs( zDataMax( zMatBuf(m), zMatRowSizeNC(m)*zMatColSizeNC(m), NULL ) );
  min = fabs( zDataMin( zMatBuf(m), zMatRowSizeNC(m)*zMatColSizeNC(m), NULL ) );
  d = zMax( max, min ) / 4;
  for( i=0; i<zMatRowSizeNC(m); i++ ){
    for( j=0; j<zMatColSizeNC(m); j++ ){
      if( zIsTiny( zMatElem(m,i,j) ) )
        printf( "  " );
      else{
        c = zRound( zMatElem(m,i,j) / d );
        c = zLimit( c, -3, 3 );
        printf( "%c ", c < 0 ? pat_neg[-c] : pat_pos[c] );
      }
    }
    zEndl();
  }
}
