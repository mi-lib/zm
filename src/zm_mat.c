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

/* set matrix components from value list. */
zMat _zMatSetElemList(zMat m, va_list args)
{
  register int i, j;

  for( i=0; i<zMatRowSizeNC(m); i++ )
    for( j=0; j<zMatColSizeNC(m); j++ )
      zMatSetElemNC( m, i, j, (double)va_arg( args, double ) );
  return m;
}

/* set matrix components from value list. */
zMat zMatSetElemList(zMat m, ... )
{
  va_list args;

  va_start( args, m );
  _zMatSetElemList( m, args );
  va_end( args );
  return m;
}

/* allocate memory for a matrix. */
zMat zMatAlloc(int row, int col)
{
  zMat m;

  if( !( m = zAlloc( zMatStruct, 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  if( !( zMatBufNC(m) = zAlloc( double, row*col ) ) ){
    ZALLOCERROR();
    zFree( m );
    return NULL;
  }
  zMatSetSizeNC( m, row, col );
  return m;
}

/* create a matrix from value list. */
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

/* free a matrix. */
void zMatFree(zMat m)
{
  if( m ){
    zFree( zMatBufNC(m) );
    free( m );
  }
}

/* free matrices at once. */
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

/* zero a matrix. */
zMat zMatZero(zMat m)
{
  zRawMatZero( zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
  return m;
}

/* touchup matrix. */
zMat zMatTouchup(zMat m)
{
  zRawMatTouchup( zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
  return m;
}

/* create identity matrix without checking size consistency. */
zMat zMatIdentNC(zMat m)
{
  zRawMatIdent( zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
  return m;
}

/* create diagonal matrix without checking size consistency. */
zMat zMatDiagNC(zMat m, zVec d)
{
  zRawMatDiag( zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m), zVecBufNC(d) );
  return m;
}

/* create identity matrix. */
zMat zMatIdent(zMat m)
{
  if( !zMatIsSqr( m ) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return NULL;
  }
  return zMatIdentNC( m );
}

/* create diagonal matrix. */
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

/* create a random matrix with a uniform range. */
zMat zMatRandUniform(zMat m, double min, double max)
{
  zRawMatRandUniform( zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m), min, max );
  return m;
}

/* create a random matrix with range matrices. */
zMat zMatRand(zMat m, zMat min, zMat max)
{
  zRawMatRand( zMatBufNC(m), zMatBufNC(min), zMatBufNC(max), zMatRowSizeNC(m), zMatColSizeNC(m) );
  return m;
}

/* copy a matrix without checking size consistency. */
zMat zMatCopyNC(zMat src, zMat dest)
{
  zRawMatCopy( zMatBufNC(src), zMatBufNC(dest), zMatRowSizeNC(dest), zMatColSizeNC(dest) );
  return dest;
}

/* copy a matrix. */
zMat zMatCopy(zMat src, zMat dest)
{
  if( !zMatSizeIsEqual( src, dest ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMatCopyNC( src, dest );
}

/* copy matrix from 2-dim array of double precision floating point values. */
zMat zMatCopyArray(double array[], int r, int c, zMat m)
{
  if( zMatRowSizeNC(m) != r || zMatColSizeNC(m) != c ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  zRawMatCopy( array, zMatBufNC(m), r, c );
  return m;
}

/* clone a matrix. */
zMat zMatClone(zMat src)
{
  zMat dest;

  if( ( dest = zMatAlloc( zMatRowSizeNC(src), zMatColSizeNC(src) ) ) )
    zMatCopyNC( src, dest );
  return dest;
}

/* create a matrix from an array of double precision floating point values. */
zMat zMatCloneArray(double array[], int r, int c)
{
  zMat m;

  if( ( m = zMatAlloc( r, c ) ) )
    zMatCopyArray( array, r, c, m );
  return m;
}

/* get submatrix without checking validity of size. */
zMat zMatGetNC(zMat src, int pr, int pc, zMat dest)
{
  zRawMatGet( zMatBufNC(src), zMatRowSizeNC(src), zMatColSizeNC(src), pr, pc, zMatBufNC(dest), zMatRowSizeNC(dest), zMatColSizeNC(dest) );
  return dest;
}

/* get submatrix. */
zMat zMatGet(zMat src, int pr, int pc, zMat dest)
{
  if( pr + zMatRowSizeNC(dest) > zMatRowSizeNC(src) ||
      pc + zMatColSizeNC(dest) > zMatColSizeNC(src) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMatGetNC( src, pr, pc, dest );
}

/* get transpose of a submatrix without checking validity of size. */
zMat zMatTGetNC(zMat src, int pr, int pc, zMat dest)
{
  zRawMatTGet( zMatBufNC(src), zMatRowSizeNC(src), zMatColSizeNC(src), pr, pc, zMatBufNC(dest), zMatRowSizeNC(dest), zMatColSizeNC(dest) );
  return dest;
}

/* get transpose of a submatrix. */
zMat zMatTGet(zMat src, int pr, int pc, zMat dest)
{
  if( pr + zMatColSizeNC(dest) > zMatRowSizeNC(src) ||
      pc + zMatRowSizeNC(dest) > zMatColSizeNC(src) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMatTGetNC( src, pr, pc, dest );
}

/* put submatrix without checking the size validity. */
zMat zMatPutNC(zMat dest, int pr, int pc, zMat src)
{
  zRawMatPut( zMatBufNC(dest), zMatRowSizeNC(dest), zMatColSizeNC(dest), pr, pc, zMatBufNC(src), zMatRowSizeNC(src), zMatColSizeNC(src) );
  return dest;
}

/* put submatrix. */
zMat zMatPut(zMat dest, int pr, int pc, zMat src)
{
  if( pr + zMatRowSizeNC(src) > zMatRowSizeNC(dest) ||
      pc + zMatColSizeNC(src) > zMatColSizeNC(dest) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMatPutNC( dest, pr, pc, src );
}

/* put transpose of a submatrix without checking validity of size. */
zMat zMatTPutNC(zMat dest, int pr, int pc, zMat src)
{
  zRawMatTPut( zMatBufNC(dest), zMatRowSizeNC(dest), zMatColSizeNC(dest), pr, pc, zMatBufNC(src), zMatRowSizeNC(src), zMatColSizeNC(src) );
  return dest;
}

/* put transpose of a submatrix. */
zMat zMatTPut(zMat dest, int pr, int pc, zMat src)
{
  if( pr + zMatColSizeNC(src) > zMatRowSizeNC(dest) ||
      pc + zMatRowSizeNC(src) > zMatColSizeNC(dest) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMatTPutNC( dest, pr, pc, src );
}

/* abstract row vector of matrix without checking size consistency. */
zVec zMatGetRowNC(zMat m, int row, zVec v)
{
  zRawMatGetRow( zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m), row, zVecBufNC(v) );
  return v;
}

/* abstract column vector of matrix without checking size consistency. */
zVec zMatGetColNC(zMat m, int col, zVec v)
{
  zRawMatGetCol( zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m), col, zVecBufNC(v) );
  return v;
}

/* abstract row vector of matrix. */
zVec zMatGetRow(zMat m, int row, zVec v)
{
  if( !zMatColVecSizeIsEqual( m, v ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  if( row >= zMatRowSizeNC(m) ){
    ZRUNERROR( ZM_ERR_INV_ROW );
    return NULL;
  }
  return zMatGetRowNC( m, row, v );
}

/* abstract column vector of matrix. */
zVec zMatGetCol(zMat m, int col, zVec v)
{
  if( !zMatRowVecSizeIsEqual( m, v ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  if( col >= zMatColSizeNC(m) ){
    ZRUNERROR( ZM_ERR_INV_COL );
    return NULL;
  }
  return zMatGetColNC( m, col, v );
}

/* put a row vector to matrix without checking size consistency. */
zMat zMatPutRowNC(zMat m, int row, zVec v)
{
  zRawMatPutRow( zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m), row, zVecBufNC(v) );
  return m;
}

/* put a column vector to matrix without checking size consistency. */
zMat zMatPutColNC(zMat m, int col, zVec v)
{
  zRawMatPutCol( zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m), col, zVecBufNC(v) );
  return m;
}

/* put a row vector to matrix. */
zMat zMatPutRow(zMat m, int row, zVec v)
{
  if( !zMatColVecSizeIsEqual( m, v ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  if( row >= zMatRowSizeNC(m) ){
    ZRUNERROR( ZM_ERR_INV_ROW );
    return NULL;
  }
  return zMatPutRowNC( m, row, v );
}

/* put a column vector to matrix. */
zMat zMatPutCol(zMat m, int col, zVec v)
{
  if( !zMatRowVecSizeIsEqual( m, v ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  if( col >= zMatColSizeNC(m) ){
    ZRUNERROR( ZM_ERR_INV_COL );
    return NULL;
  }
  return zMatPutColNC( m, col, v );
}

/* swap two matrix rows without checking size consistency. */
zMat zMatSwapRowNC(zMat m, int r1, int r2)
{
  zRawMatSwapRow( zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m), r1, r2 );
  return m;
}

/* swap two rows of a matrix. */
zMat zMatSwapRow(zMat m, int r1, int r2)
{
  if( r1 >= zMatRowSizeNC(m) || r2 >= zMatRowSizeNC(m) ){
    ZRUNERROR( ZM_ERR_INV_ROW );
    return NULL;
  }
  return zMatSwapRowNC( m, r1, r2 );
}

/* swap two columns of a matrix without checking size consistency. */
zMat zMatSwapColNC(zMat m, int c1, int c2)
{
  zRawMatSwapCol( zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m), c1, c2 );
  return m;
}

/* swap two matrix columns. */
zMat zMatSwapCol(zMat m, int c1, int c2)
{
  if( c1 >= zMatColSizeNC(m) || c2 >= zMatColSizeNC(m) ){
    ZRUNERROR( ZM_ERR_INV_COL );
    return NULL;
  }
  return zMatSwapColNC( m, c1, c2 );
}

/* shift diagonal values of a matrix. */
void zMatShift(zMat m, double shift)
{
  register int i;

  for( i=0; i<zMatRowSizeNC(m); i++ )
    zMatElemNC( m, i, i ) += shift;
}

/* see if two matrices are equal. */
bool zMatIsEqual(zMat m1, zMat m2)
{
  register int i, j;

  if( !zMatSizeIsEqual( m1, m2 ) ) return false;
  for( i=0; i<zMatRowSizeNC(m1); i++ )
    for( j=0; j<zMatColSizeNC(m1); j++ )
      if( !zIsTiny( zMatElemNC(m1,i,j)/zMatElemNC(m2,i,j) - 1.0 ) ) return false;
  return true;
}

/* test if a matrix is tiny. */
bool zMatIsTol(zMat m, double tol)
{
  return zRawMatIsTol( zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m), tol );
}

/* row regression of a matrix. */
zMat zMatRowReg(zMat m, int rank)
{
  if( rank < zMatRowSizeNC(m) ) zMatSetRowSizeNC( m, rank );
  return m;
}

/* column regression of a matrix. */
zMat zMatColReg(zMat m, int rank)
{
  register int i;
  double *sp, *dp;

  if( rank < zMatColSizeNC(m) ){
    for( sp=dp=zMatBufNC(m), i=0; i<zMatRowSizeNC(m); i++, sp+=zMatColSizeNC(m), dp+=rank )
      memmove( dp, sp, sizeof(double)*rank );
    zMatSetColSizeNC( m, rank );
  }
  return m;
}

/* add matrices without checking size consistency. */
zMat zMatAddNC(zMat m1, zMat m2, zMat m)
{
  zRawMatAdd( zMatBufNC(m1), zMatBufNC(m2), zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
  return m;
}

/* subtract a matrix from nother without checking size consistency. */
zMat zMatSubNC(zMat m1, zMat m2, zMat m)
{
  zRawMatSub( zMatBufNC(m1), zMatBufNC(m2), zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
  return m;
}

/* reverse a matrix without checking size consistency. */
zMat zMatRevNC(zMat m1, zMat m)
{
  zRawMatRev( zMatBufNC(m1), zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
  return m;
}

/* multiply a matrix by a scalar value without checking size consistency. */
zMat zMatMulNC(zMat m1, double k, zMat m)
{
  zRawMatMul( zMatBufNC(m1), k, zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
  return m;
}

/* divide a matrix by a scalar value rwithout checking size consistency. */
zMat zMatDivNC(zMat m1, double k, zMat m)
{
  zRawMatDiv( zMatBufNC(m1), k, zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
  return m;
}

/* concatenate a matrix with another without checking size consistency. */
zMat zMatCatNC(zMat m1, double k, zMat m2, zMat m)
{
  zRawMatCat( zMatBufNC(m1), k, zMatBufNC(m2), zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
  return m;
}

/* add matrices. */
zMat zMatAdd(zMat m1, zMat m2, zMat m)
{
  if( !zMatSizeIsEqual(m1,m2) || !zMatSizeIsEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMatAddNC( m1, m2, m );
}

/* substract a matrix from another. */
zMat zMatSub(zMat m1, zMat m2, zMat m)
{
  if( !zMatSizeIsEqual(m1,m2) || !zMatSizeIsEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMatSubNC( m1, m2, m );
}

/* reverse a matrix. */
zMat zMatRev(zMat m1, zMat m)
{
  if( !zMatSizeIsEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMatRevNC( m1, m );
}

/* multiply a matrix by a scalar value. */
zMat zMatMul(zMat m1, double k, zMat m)
{
  if( !zMatSizeIsEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMatMulNC( m1, k, m );
}

/* divide a matrix by a scalar value. */
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

/* concatenate a matrix with another. */
zMat zMatCat(zMat m1, double k, zMat m2, zMat m)
{
  if( !zMatSizeIsEqual(m1,m2) || !zMatSizeIsEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMatCatNC( m1, k, m2, m );
}

/* squared norm of a matrix. */
double zMatSqrNorm(zMat m)
{
  return zRawMatSqrNorm( zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
}

/* infinity norm of a matrix. */
double zMatInfNorm(zMat m)
{
  double *mp, rs, rsmax = 0;
  register int i, j;

  mp = zMatBufNC(m);
  for( i=0; i<zMatRowSizeNC(m); i++ ){
    for( rs=0, j=0; j<zMatColSizeNC(m); j++ )
      rs += fabs( *mp++ );
    if( rs > rsmax ) rsmax = rs;
  }
  return rsmax;
}

/* transpose a matrix without checking size consistency. */
zMat zMatTNC(zMat m, zMat tm)
{
  zRawMatT( zMatBufNC(m), zMatBufNC(tm), zMatRowSizeNC(tm), zMatColSizeNC(tm) );
  return tm;
}

/* transpose a matrix. */
zMat zMatT(zMat m, zMat tm)
{
  if( !zMatColRowSizeIsEqual( tm, m ) || !zMatColRowSizeIsEqual( m, tm ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMatTNC( m, tm );
}

/* transpose a matrix directly. */
zMat zMatTDRC(zMat m)
{
  int row, col;

  row = zMatRowSizeNC(m);
  col = zMatColSizeNC(m);
  zRawMatTDRC( zMatBufNC(m), row, col );
  zMatSetSizeNC( m, col, row );
  return m;
}

/* clone transpose of a matrix. */
zMat zMatTClone(zMat src)
{
  zMat dest;

  if( ( dest = zMatAlloc( zMatColSizeNC(src), zMatRowSizeNC(src) ) ) )
    zMatTNC( src, dest );
  return dest;
}

/* dyadic product of vectors without checking size consistency. */
zMat zVecDyadNC(zVec v1, zVec v2, zMat dyad)
{
  zRawVecDyad( zVecBufNC(v1), zVecSizeNC(v1), zVecBufNC(v2), zVecSizeNC(v2), zMatBufNC(dyad) );
  return dyad;
}

/* dyadic product of vectors. */
zMat zVecDyad(zVec v1, zVec v2, zMat dyad)
{
  if( zVecSizeNC(v1) != zMatRowSizeNC(dyad) || zVecSizeNC(v2) != zMatColSizeNC(dyad) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  return zVecDyadNC( v1, v2, dyad );
}

/* add dyadic product of vectors to a matrix without checking size consistency. */
zMat zMatAddDyadNC(zMat m, zVec v1, zVec v2)
{
  zRawMatAddDyad( zMatBufNC(m), zVecBufNC(v1), zVecSizeNC(v1), zVecBufNC(v2), zVecSizeNC(v2) );
  return m;
}

/* add dyadic product of vectors to a matrix. */
zMat zMatAddDyad(zMat m, zVec v1, zVec v2)
{
  if( zVecSizeNC(v1) != zMatRowSizeNC(m) || zVecSizeNC(v2) != zMatColSizeNC(m) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  return zMatAddDyadNC( m, v1, v2 );
}

/* subtract dyadic product of vectors from a matrix without checking size consistency. */
zMat zMatSubDyadNC(zMat m, zVec v1, zVec v2)
{
  zRawMatSubDyad( zMatBufNC(m), zVecBufNC(v1), zVecSizeNC(v1), zVecBufNC(v2), zVecSizeNC(v2) );
  return m;
}

/* subtract dyadic product of vectors from a matrix. */
zMat zMatSubDyad(zMat m, zVec v1, zVec v2)
{
  if( zVecSizeNC(v1) != zMatRowSizeNC(m) || zVecSizeNC(v2) != zMatColSizeNC(m) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  return zMatSubDyadNC( m, v1, v2 );
}

/* add dyadic product of vectors multiplied by a scalar value to a matrix without checking size consistency. */
zMat zMatCatDyadNC(zMat m, double k, zVec v1, zVec v2)
{
  zRawMatCatDyad( zMatBufNC(m), k, zVecBufNC(v1), zVecSizeNC(v1), zVecBufNC(v2), zVecSizeNC(v2) );
  return m;
}

/* add dyadic product of vectors multiplied by a scalar to a matrix. */
zMat zMatCatDyad(zMat m, double k, zVec v1, zVec v2)
{
  if( zVecSizeNC(v1) != zMatRowSizeNC(m) || zVecSizeNC(v2) != zMatColSizeNC(m) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  return zMatCatDyadNC( m, k, v1, v2 );
}

/* trace of a matrix without checking size consistency. */
double zMatTrNC(zMat m)
{
  return zRawMatTr( zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
}

/* trace of a matrix. */
double zMatTr(zMat m)
{
  if( !zMatIsSqr( m ) ){
    ZRUNERROR( ZM_ERR_NONSQR_MAT );
    return 0;
  }
  return zMatTrNC( m );
}

/* multiply a vector by a matrix from the left side without checking size consistency. */
zVec zMulMatVecNC(zMat m, zVec v1, zVec v)
{
  zRawMulMatVec( zMatBufNC(m), zVecBufNC(v1), zMatRowSizeNC(m), zMatColSizeNC(m), zVecBufNC(v) );
  return v;
}

/* multiply a vector by transpose of a matrix from the left side without checking size consistency. */
zVec zMulMatTVecNC(zMat m, zVec v1, zVec v)
{
  zRawMulMatTVec( zMatBufNC(m), zVecBufNC(v1), zMatRowSizeNC(m), zMatColSizeNC(m), zVecBufNC(v) );
  return v;
}

/* multiply two matrices without checking size consistency. */
zMat zMulMatMatNC(zMat m1, zMat m2, zMat m)
{
  zRawMulMatMat( zMatBufNC(m1),zMatRowSizeNC(m1), zMatColSizeNC(m1),
    zMatBufNC(m2), zMatRowSizeNC(m2), zMatColSizeNC(m2), zMatBufNC(m) );
  return m;
}

/* multiply a matrix and transpose of a matrix ('m = m1 m2^T') without checking size consistency. */
zMat zMulMatMatTNC(zMat m1, zMat m2, zMat m)
{
  zRawMulMatMatT( zMatBufNC(m1), zMatRowSizeNC(m1), zMatColSizeNC(m1),
    zMatBufNC(m2), zMatRowSizeNC(m2), zMatColSizeNC(m2), zMatBufNC(m) );
  return m;
}

/* multiply transpose of a matrix and a matrix 'm = m1^T m2' without checking size consistency. */
zMat zMulMatTMatNC(zMat m1, zMat m2, zMat m)
{
  zRawMulMatTMat( zMatBufNC(m1), zMatRowSizeNC(m1), zMatColSizeNC(m1),
    zMatBufNC(m2), zMatRowSizeNC(m2), zMatColSizeNC(m2), zMatBufNC(m) );
  return m;
}

/* multiply a vector by a matrix from the left side. */
zVec zMulMatVec(zMat m, zVec v1, zVec v)
{
  if( !zMatColVecSizeIsEqual( m, v1 ) || !zMatRowVecSizeIsEqual( m, v ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  return zMulMatVecNC( m, v1, v );
}

/* multiply a vector by transpose of a matrix from the left side. */
zVec zMulMatTVec(zMat m, zVec v1, zVec v)
{
  if( !zMatRowVecSizeIsEqual( m, v1 ) || !zMatColVecSizeIsEqual( m, v ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  return zMulMatTVecNC( m, v1, v );
}

/* multiply two matrices. */
zMat zMulMatMat(zMat m1, zMat m2, zMat m)
{
  if( !zMatColRowSizeIsEqual( m1, m2 ) ||
      !zMatRowSizeIsEqual( m1, m ) || !zMatColSizeIsEqual( m2, m ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMulMatMatNC( m1, m2, m );
}

/* multiply a matrix by transpose of another matrix from the right side. */
zMat zMulMatMatT(zMat m1, zMat m2, zMat m)
{
  if( !zMatColSizeIsEqual( m1, m2 ) || !zMatRowSizeIsEqual( m1, m ) ||
      !zMatRowColSizeIsEqual( m2, m ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMulMatMatTNC( m1, m2, m );
}

/* multiply a matrix by transpose of another matrix from the left side. */
zMat zMulMatTMat(zMat m1, zMat m2, zMat m)
{
  if( !zMatRowSizeIsEqual( m1, m2 ) ||
      !zMatColRowSizeIsEqual( m1, m ) || !zMatColSizeIsEqual( m2, m ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  return zMulMatTMatNC( m1, m2, m );
}

/* multiply a vector by a matrix directly.
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

/* multiply a vector by transpose of a matrix directly. */
zVec zMulMatTVecDRC(zMat m, zVec v)
{
  zVec tmp;

  if( !( tmp = zVecAlloc( zVecSizeNC(v) ) ) ) return NULL;
  zMulMatTVec( m, v, tmp );
  zVecCopyNC( tmp, v );
  zVecFree( tmp );
  return v;
}

/* quadratic multiplication of matrices ('q = a diag{w} a^T') without checking size consistency. */
zMat zMatQuadNC(zMat a, zVec w, zMat q)
{
  register int i, j, k;
  double wa;

  zMatZero( q );
  for( k=0; k<zMatColSizeNC(a); k++ )
    for( i=0; i<zMatRowSizeNC(a); i++ ){
      wa = w ? zVecElem(w,k) * zMatElemNC(a,i,k) : zMatElemNC(a,i,k);
      for( j=i; j<zMatRowSizeNC(a); j++ )
        zMatElemNC(q,i,j) += wa * zMatElemNC(a,j,k);
    }
  for( i=0; i<zMatRowSizeNC(a); i++ )
    for( j=i; j<zMatRowSizeNC(a); j++ )
      zMatSetElemNC( q, j, i, zMatElemNC(q,i,j) );
  return q;
}

/* quadratic multiplication of matrices ('q = a diag{w} a^T'). */
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

/* quadratic multiplication of matrices ('q = a^T diag{w} a') without checking size consistency. */
zMat zMatTQuadNC(zMat a, zVec w, zMat q)
{
  register int i, j, k;
  double wa;

  zMatZero( q );
  for( k=0; k<zMatRowSizeNC(a); k++ )
    for( i=0; i<zMatColSizeNC(a); i++ ){
      wa = w ? zVecElem(w,k) * zMatElemNC(a,k,i) : zMatElemNC(a,k,i);
      for( j=i; j<zMatColSizeNC(a); j++ )
        zMatElemNC(q,i,j) += wa * zMatElemNC(a,k,j);
    }
  for( i=0; i<zMatColSizeNC(a); i++ )
    for( j=i; j<zMatColSizeNC(a); j++ )
      zMatSetElemNC( q, j, i, zMatElemNC(q,i,j) );
  return q;
}

/* quadratic multiplication of matrices ('q = a^T diag{w} a'). */
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

/* ********************************************************** */
/* I/O
 * ********************************************************** */

/* read a matrix from a ZTK format processor. */
zMat zMatFromZTK(ZTK *ztk)
{
  register int i, j, row, col;
  zMat m;

  row = ZTKInt(ztk);
  col = ZTKInt(ztk);
  if( !( m = zMatAlloc( row, col ) ) ) return NULL;
  for( i=0; i<row; i++ )
    for( j=0; j<col; j++ )
      zMatSetElemNC( m, i, j, ZTKDouble(ztk) );
  return m;
}

/* scan information of a matrix from file. */
zMat zMatFScan(FILE *fp)
{
  register unsigned i, j, row, col;
  zMat m;

  row = zFInt( fp );
  col = zFInt( fp );
  if( !( m = zMatAlloc( row, col ) ) ) return NULL;
  for( i=0; i<row; i++ )
    for( j=0; j<col; j++ )
      zMatSetElemNC( m, i, j, zFDouble( fp ) );
  return m;
}

/* print information of a matrix to file. */
void zMatFPrint(FILE *fp, zMat m)
{
  register int i, j;

  if( !m )
    fprintf( fp, "(null matrix)\n" );
  else{
    fprintf( fp,"(%d, %d) {\n",
      zMatRowSizeNC(m), zMatColSizeNC(m) );
    for( i=0; i<zMatRowSizeNC(m); i++ ){
      for( j=0; j<zMatColSizeNC(m); j++ )
        fprintf( fp, " %.10g", zMatElemNC( m, i, j ) );
      fprintf( fp, "\n" );
    }
    fprintf( fp, "}\n" );
  }
}

/* visualize a matrix using one-charactor collage for debug. */
void zMatImg(zMat m)
{
  double max, min, d;
  const char pat_pos[] = ",x*M";
  const char pat_neg[] = ".oO@";
  register int i, j, c;

  max = fabs( zDataMax( zMatBufNC(m), zMatRowSizeNC(m)*zMatColSizeNC(m), NULL ) );
  min = fabs( zDataMin( zMatBufNC(m), zMatRowSizeNC(m)*zMatColSizeNC(m), NULL ) );
  d = zMax( max, min ) / 4;
  for( i=0; i<zMatRowSizeNC(m); i++ ){
    for( j=0; j<zMatColSizeNC(m); j++ ){
      if( zIsTiny( zMatElemNC(m,i,j) ) )
        printf( "  " );
      else{
        c = zRound( zMatElemNC(m,i,j) / d );
        c = zLimit( c, -3, 3 );
        printf( "%c ", c < 0 ? pat_neg[-c] : pat_pos[c] );
      }
    }
    zEndl();
  }
}
