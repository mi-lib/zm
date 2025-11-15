/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_vec_mat - vector and matrix class : matrix class.
 */

#include <zm/zm_mat.h>

/* set matrix components from value list. */
static zMat _zMatSetElemList(zMat m, va_list args)
{
  int i, j;

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
zMat zMatAlloc(int rowsize, int colsize)
{
  zMat m;

  if( !( m = zAlloc( zMatStruct, 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  zArray2Alloc( m, double, rowsize, colsize );
  if( !zMatBufNC(m) ){
    zFree( m );
    return NULL;
  }
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
void zMatFreeAtOnce(int num, ...)
{
  va_list arg;
  zMat m;
  int i;

  va_start( arg, num );
  for( i=0; i<num; i++ ){
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
zMat zMatTouchup(zMat m, double tol)
{
  zRawMatTouchup( zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m), tol );
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
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return NULL;
  }
  return zMatIdentNC( m );
}

/* create diagonal matrix. */
zMat zMatDiag(zMat m, zVec d)
{
  if( !zMatIsSqr( m ) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return NULL;
  }
  if( !zMatColVecSizeEqual( m, d ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
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
zMat zMatCopyNC(const zMat src, zMat dest)
{
  zRawMatCopy( zMatBufNC(src), zMatBufNC(dest), zMatRowSizeNC(dest), zMatColSizeNC(dest) );
  return dest;
}

/* copy a matrix. */
zMat zMatCopy(const zMat src, zMat dest)
{
  if( !zMatSizeEqual( src, dest ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  return zMatCopyNC( src, dest );
}

/* copy matrix from 2-dim array of double precision floating point values. */
zMat zMatCopyArray(const double array[], int r, int c, zMat m)
{
  if( zMatRowSizeNC(m) != r || zMatColSizeNC(m) != c ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  zRawMatCopy( array, zMatBufNC(m), r, c );
  return m;
}

/* clone a matrix. */
zMat zMatClone(const zMat src)
{
  zMat dest;

  if( !src ) return NULL;
  if( ( dest = zMatAlloc( zMatRowSizeNC(src), zMatColSizeNC(src) ) ) )
    zMatCopyNC( src, dest );
  return dest;
}

/* create a matrix from an array of double precision floating point values. */
zMat zMatCloneArray(const double array[], int r, int c)
{
  zMat m;

  if( ( m = zMatAlloc( r, c ) ) )
    zMatCopyArray( array, r, c, m );
  return m;
}

/* get submatrix without checking validity of size. */
zMat zMatGetNC(const zMat src, int pr, int pc, zMat dest)
{
  zRawMatGet( zMatBufNC(src), zMatRowSizeNC(src), zMatColSizeNC(src), pr, pc, zMatBufNC(dest), zMatRowSizeNC(dest), zMatColSizeNC(dest) );
  return dest;
}

/* get submatrix. */
zMat zMatGet(const zMat src, int pr, int pc, zMat dest)
{
  if( pr + zMatRowSizeNC(dest) > zMatRowSizeNC(src) ||
      pc + zMatColSizeNC(dest) > zMatColSizeNC(src) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  return zMatGetNC( src, pr, pc, dest );
}

/* get transpose of a submatrix without checking validity of size. */
zMat zMatTGetNC(const zMat src, int pr, int pc, zMat dest)
{
  zRawMatTGet( zMatBufNC(src), zMatRowSizeNC(src), zMatColSizeNC(src), pr, pc, zMatBufNC(dest), zMatRowSizeNC(dest), zMatColSizeNC(dest) );
  return dest;
}

/* get transpose of a submatrix. */
zMat zMatTGet(const zMat src, int pr, int pc, zMat dest)
{
  if( pr + zMatColSizeNC(dest) > zMatRowSizeNC(src) ||
      pc + zMatRowSizeNC(dest) > zMatColSizeNC(src) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  return zMatTGetNC( src, pr, pc, dest );
}

/* put submatrix without checking the size validity. */
zMat zMatPutNC(zMat dest, int pr, int pc, const zMat src)
{
  zRawMatPut( zMatBufNC(dest), zMatRowSizeNC(dest), zMatColSizeNC(dest), pr, pc, zMatBufNC(src), zMatRowSizeNC(src), zMatColSizeNC(src) );
  return dest;
}

/* put submatrix. */
zMat zMatPut(zMat dest, int pr, int pc, const zMat src)
{
  if( pr + zMatRowSizeNC(src) > zMatRowSizeNC(dest) ||
      pc + zMatColSizeNC(src) > zMatColSizeNC(dest) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  return zMatPutNC( dest, pr, pc, src );
}

/* put transpose of a submatrix without checking validity of size. */
zMat zMatTPutNC(zMat dest, int pr, int pc, const zMat src)
{
  zRawMatTPut( zMatBufNC(dest), zMatRowSizeNC(dest), zMatColSizeNC(dest), pr, pc, zMatBufNC(src), zMatRowSizeNC(src), zMatColSizeNC(src) );
  return dest;
}

/* put transpose of a submatrix. */
zMat zMatTPut(zMat dest, int pr, int pc, const zMat src)
{
  if( pr + zMatColSizeNC(src) > zMatRowSizeNC(dest) ||
      pc + zMatRowSizeNC(src) > zMatColSizeNC(dest) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  return zMatTPutNC( dest, pr, pc, src );
}

/* abstract row vector of matrix without checking size consistency. */
zVec zMatGetRowNC(const zMat m, int row, zVec v)
{
  zRawMatGetRow( zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m), row, zVecBufNC(v) );
  return v;
}

/* abstract column vector of matrix without checking size consistency. */
zVec zMatGetColNC(const zMat m, int col, zVec v)
{
  zRawMatGetCol( zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m), col, zVecBufNC(v) );
  return v;
}

/* abstract row vector of matrix. */
zVec zMatGetRow(const zMat m, int row, zVec v)
{
  if( !zMatColVecSizeEqual( m, v ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  if( row >= zMatRowSizeNC(m) ){
    ZRUNERROR( ZM_ERR_INVALID_ROW );
    return NULL;
  }
  return zMatGetRowNC( m, row, v );
}

/* abstract column vector of matrix. */
zVec zMatGetCol(const zMat m, int col, zVec v)
{
  if( !zMatRowVecSizeEqual( m, v ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  if( col >= zMatColSizeNC(m) ){
    ZRUNERROR( ZM_ERR_INVALID_COL );
    return NULL;
  }
  return zMatGetColNC( m, col, v );
}

/* put a row vector to matrix without checking size consistency. */
zMat zMatPutRowNC(zMat m, int row, const zVec v)
{
  zRawMatPutRow( zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m), row, zVecBufNC(v) );
  return m;
}

/* put a column vector to matrix without checking size consistency. */
zMat zMatPutColNC(zMat m, int col, const zVec v)
{
  zRawMatPutCol( zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m), col, zVecBufNC(v) );
  return m;
}

/* put a row vector to matrix. */
zMat zMatPutRow(zMat m, int row, const zVec v)
{
  if( !zMatColVecSizeEqual( m, v ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  if( row >= zMatRowSizeNC(m) ){
    ZRUNERROR( ZM_ERR_INVALID_ROW );
    return NULL;
  }
  return zMatPutRowNC( m, row, v );
}

/* put a column vector to matrix. */
zMat zMatPutCol(zMat m, int col, const zVec v)
{
  if( !zMatRowVecSizeEqual( m, v ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  if( col >= zMatColSizeNC(m) ){
    ZRUNERROR( ZM_ERR_INVALID_COL );
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
    ZRUNERROR( ZM_ERR_INVALID_ROW );
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
    ZRUNERROR( ZM_ERR_INVALID_COL );
    return NULL;
  }
  return zMatSwapColNC( m, c1, c2 );
}

/* shift diagonal values of a matrix. */
void zMatShift(zMat m, double shift)
{
  int i;

  for( i=0; i<zMatRowSizeNC(m); i++ )
    zMatElemNC( m, i, i ) += shift;
}

/* maximum of matrix elements. */
double zMatElemMax(const zMat m, int *im){ return _zMatElemMax( m, im ); }
/* minimum of vector elements. */
double zMatElemMin(const zMat m, int *im){ return _zMatElemMin( m, im ); }
/* absolute maximum of vector elements. */
double zMatElemAbsMax(const zMat m, int *im){ return _zMatElemAbsMax( m, im ); }
/* absolute minimum of vector elements. */
double zMatElemAbsMin(const zMat m, int *im){ return _zMatElemAbsMin( m, im ); }

/* check if two matrices are equal. */
bool zMatEqual(const zMat m1, const zMat m2, double tol)
{
  int i, j;

  if( !zMatSizeEqual( m1, m2 ) ) return false;
  for( i=0; i<zMatRowSizeNC(m1); i++ )
    for( j=0; j<zMatColSizeNC(m1); j++ )
      if( !zEqual( zMatElemNC(m1,i,j), zMatElemNC(m2,i,j), tol ) ) return false;
  return true;
}

/* check if two matrices exactly match with each other. */
bool zMatMatch(const zMat m1, const zMat m2)
{
  int i, j;

  if( !zMatSizeEqual( m1, m2 ) ) return false;
  for( i=0; i<zMatRowSizeNC(m1); i++ )
    for( j=0; j<zMatColSizeNC(m1); j++ )
      if( zMatElemNC(m1,i,j) != zMatElemNC(m2,i,j) ) return false;
  return true;
}

/* test if a matrix is tiny. */
bool zMatIsTol(const zMat m, double tol)
{
  return zRawMatIsTol( zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m), tol );
}

/* test if a matrix is diagonal. */
bool zMatIsDiag(zMat m, double tol)
{
  int i, j;

  for( i=0; i<zMatRowSizeNC(m); i++ )
    for( j=0; j<zMatColSizeNC(m); j++ ){
      if( j == i ) continue;
      if( !zIsTol( zMatElemNC(m,i,j), tol ) ) return false;
    }
  return true;
}

/* test if a matrix is the identity matrix. */
bool zMatIsIdent(zMat m, double tol)
{
  int i, j;

  for( i=0; i<zMatRowSizeNC(m); i++ )
    for( j=0; j<zMatColSizeNC(m); j++ ){
      if( j == i ){
        if( !zEqual( zMatElemNC(m,i,j), 1.0, tol ) ) return false;
      } else
        if( !zIsTol( zMatElemNC(m,i,j), tol ) ) return false;
    }
  return true;
}

/* check if a matrix is square and symmetric. */
bool zMatIsSymmetric(const zMat m)
{
  int i, j;

  if( !zMatIsSqr(m) ) return false;
  for( i=0; i<zMatRowSizeNC(m); i++ )
    for( j=i+1; j<zMatRowSizeNC(m); j++ )
      if( zMatElemNC(m,i,j) != zMatElemNC(m,j,i) ) return false;
  return true;
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
  int i;
  double *sp, *dp;

  if( rank < zMatColSizeNC(m) ){
    for( sp=dp=zMatBufNC(m), i=0; i<zMatRowSizeNC(m); i++, sp+=zMatColSizeNC(m), dp+=rank )
      memmove( dp, sp, sizeof(double)*rank );
    zMatSetColSizeNC( m, rank );
  }
  return m;
}

/* add matrices without checking size consistency. */
zMat zMatAddNC(const zMat m1, const zMat m2, zMat m)
{
  zRawMatAdd( zMatBufNC(m1), zMatBufNC(m2), zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
  return m;
}

/* subtract a matrix from nother without checking size consistency. */
zMat zMatSubNC(const zMat m1, const zMat m2, zMat m)
{
  zRawMatSub( zMatBufNC(m1), zMatBufNC(m2), zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
  return m;
}

/* reverse a matrix without checking size consistency. */
zMat zMatRevNC(const zMat m1, zMat m)
{
  zRawMatRev( zMatBufNC(m1), zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
  return m;
}

/* multiply a matrix by a scalar value without checking size consistency. */
zMat zMatMulNC(const zMat m1, double k, zMat m)
{
  zRawMatMul( zMatBufNC(m1), k, zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
  return m;
}

/* divide a matrix by a scalar value rwithout checking size consistency. */
zMat zMatDivNC(const zMat m1, double k, zMat m)
{
  zRawMatDiv( zMatBufNC(m1), k, zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
  return m;
}

/* concatenate a matrix with another without checking size consistency. */
zMat zMatCatNC(const zMat m1, double k, const zMat m2, zMat m)
{
  zRawMatCat( zMatBufNC(m1), k, zMatBufNC(m2), zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
  return m;
}

/* add matrices. */
zMat zMatAdd(const zMat m1, const zMat m2, zMat m)
{
  if( !zMatSizeEqual(m1,m2) || !zMatSizeEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  return zMatAddNC( m1, m2, m );
}

/* substract a matrix from another. */
zMat zMatSub(const zMat m1, const zMat m2, zMat m)
{
  if( !zMatSizeEqual(m1,m2) || !zMatSizeEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  return zMatSubNC( m1, m2, m );
}

/* reverse a matrix. */
zMat zMatRev(const zMat m1, zMat m)
{
  if( !zMatSizeEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  return zMatRevNC( m1, m );
}

/* multiply a matrix by a scalar value. */
zMat zMatMul(const zMat m1, double k, zMat m)
{
  if( !zMatSizeEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  return zMatMulNC( m1, k, m );
}

/* divide a matrix by a scalar value. */
zMat zMatDiv(const zMat m1, double k, zMat m)
{
  if( !zMatSizeEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  if( zIsTiny( k ) ){
    ZRUNERROR( ZM_ERR_ZERODIV );
    return NULL;
  }
  return zMatDivNC( m1, k, m );
}

/* concatenate a matrix with another. */
zMat zMatCat(const zMat m1, double k, const zMat m2, zMat m)
{
  if( !zMatSizeEqual(m1,m2) || !zMatSizeEqual(m1,m) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  return zMatCatNC( m1, k, m2, m );
}

/* squared norm of a matrix. */
double zMatSqrNorm(const zMat m)
{
  return zRawMatSqrNorm( zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
}

/* infinity norm of a matrix. */
double zMatInfNorm(const zMat m)
{
  double *mp, rs, rsmax = 0;
  int i, j;

  mp = zMatBufNC(m);
  for( i=0; i<zMatRowSizeNC(m); i++ ){
    for( rs=0, j=0; j<zMatColSizeNC(m); j++ )
      rs += fabs( *mp++ );
    if( rs > rsmax ) rsmax = rs;
  }
  return rsmax;
}

/* transpose a matrix without checking size consistency. */
zMat zMatTNC(const zMat m, zMat tm)
{
  zRawMatT( zMatBufNC(m), zMatBufNC(tm), zMatRowSizeNC(tm), zMatColSizeNC(tm) );
  return tm;
}

/* transpose a matrix. */
zMat zMatT(const zMat m, zMat tm)
{
  if( !zMatColRowSizeEqual( tm, m ) || !zMatColRowSizeEqual( m, tm ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
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
zMat zMatTClone(const zMat src)
{
  zMat dest;

  if( !src ) return NULL;
  if( ( dest = zMatAlloc( zMatColSizeNC(src), zMatRowSizeNC(src) ) ) )
    zMatTNC( src, dest );
  return dest;
}

/* dyadic product of vectors without checking size consistency. */
zMat zVecDyadNC(const zVec v1, const zVec v2, zMat dyad)
{
  zRawVecDyad( zVecBufNC(v1), zVecSizeNC(v1), zVecBufNC(v2), zVecSizeNC(v2), zMatBufNC(dyad) );
  return dyad;
}

/* dyadic product of vectors. */
zMat zVecDyad(const zVec v1, const zVec v2, zMat dyad)
{
  if( zVecSizeNC(v1) != zMatRowSizeNC(dyad) || zVecSizeNC(v2) != zMatColSizeNC(dyad) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  return zVecDyadNC( v1, v2, dyad );
}

/* add dyadic product of vectors to a matrix without checking size consistency. */
zMat zMatAddDyadNC(zMat m, const zVec v1, const zVec v2)
{
  zRawMatAddDyad( zMatBufNC(m), zVecBufNC(v1), zVecSizeNC(v1), zVecBufNC(v2), zVecSizeNC(v2) );
  return m;
}

/* add dyadic product of vectors to a matrix. */
zMat zMatAddDyad(zMat m, const zVec v1, const zVec v2)
{
  if( zVecSizeNC(v1) != zMatRowSizeNC(m) || zVecSizeNC(v2) != zMatColSizeNC(m) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  return zMatAddDyadNC( m, v1, v2 );
}

/* subtract dyadic product of vectors from a matrix without checking size consistency. */
zMat zMatSubDyadNC(zMat m, const zVec v1, const zVec v2)
{
  zRawMatSubDyad( zMatBufNC(m), zVecBufNC(v1), zVecSizeNC(v1), zVecBufNC(v2), zVecSizeNC(v2) );
  return m;
}

/* subtract dyadic product of vectors from a matrix. */
zMat zMatSubDyad(zMat m, const zVec v1, const zVec v2)
{
  if( zVecSizeNC(v1) != zMatRowSizeNC(m) || zVecSizeNC(v2) != zMatColSizeNC(m) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  return zMatSubDyadNC( m, v1, v2 );
}

/* add dyadic product of vectors multiplied by a scalar value to a matrix without checking size consistency. */
zMat zMatCatDyadNC(zMat m, double k, const zVec v1, const zVec v2)
{
  zRawMatCatDyad( zMatBufNC(m), k, zVecBufNC(v1), zVecSizeNC(v1), zVecBufNC(v2), zVecSizeNC(v2) );
  return m;
}

/* add dyadic product of vectors multiplied by a scalar to a matrix. */
zMat zMatCatDyad(zMat m, double k, const zVec v1, const zVec v2)
{
  if( zVecSizeNC(v1) != zMatRowSizeNC(m) || zVecSizeNC(v2) != zMatColSizeNC(m) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  return zMatCatDyadNC( m, k, v1, v2 );
}

/* trace of a matrix without checking size consistency. */
double zMatTraceNC(const zMat m)
{
  return zRawMatTrace( zMatBufNC(m), zMatRowSizeNC(m), zMatColSizeNC(m) );
}

/* trace of a matrix. */
double zMatTrace(const zMat m)
{
  if( !zMatIsSqr( m ) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return 0;
  }
  return zMatTraceNC( m );
}

/* multiply a vector by a matrix from the left side without checking size consistency. */
zVec zMulMatVecNC(const zMat m, const zVec v1, zVec v)
{
  zRawMulMatVec( zMatBufNC(m), zVecBufNC(v1), zMatRowSizeNC(m), zMatColSizeNC(m), zVecBufNC(v) );
  return v;
}

/* multiply a vector by transpose of a matrix from the left side without checking size consistency. */
zVec zMulMatTVecNC(const zMat m, const zVec v1, zVec v)
{
  zRawMulMatTVec( zMatBufNC(m), zVecBufNC(v1), zMatRowSizeNC(m), zMatColSizeNC(m), zVecBufNC(v) );
  return v;
}

/* multiply two matrices without checking size consistency. */
zMat zMulMatMatNC(const zMat m1, const zMat m2, zMat m)
{
  zRawMulMatMat( zMatBufNC(m1),zMatRowSizeNC(m1), zMatColSizeNC(m1),
    zMatBufNC(m2), zMatRowSizeNC(m2), zMatColSizeNC(m2), zMatBufNC(m) );
  return m;
}

/* multiply a matrix and transpose of a matrix ('m = m1 m2^T') without checking size consistency. */
zMat zMulMatMatTNC(const zMat m1, const zMat m2, zMat m)
{
  zRawMulMatMatT( zMatBufNC(m1), zMatRowSizeNC(m1), zMatColSizeNC(m1),
    zMatBufNC(m2), zMatRowSizeNC(m2), zMatColSizeNC(m2), zMatBufNC(m) );
  return m;
}

/* multiply transpose of a matrix and a matrix 'm = m1^T m2' without checking size consistency. */
zMat zMulMatTMatNC(const zMat m1, const zMat m2, zMat m)
{
  zRawMulMatTMat( zMatBufNC(m1), zMatRowSizeNC(m1), zMatColSizeNC(m1),
    zMatBufNC(m2), zMatRowSizeNC(m2), zMatColSizeNC(m2), zMatBufNC(m) );
  return m;
}

/* multiply a vector by a matrix from the left side. */
zVec zMulMatVec(const zMat m, const zVec v1, zVec v)
{
  if( !zMatColVecSizeEqual( m, v1 ) || !zMatRowVecSizeEqual( m, v ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  return zMulMatVecNC( m, v1, v );
}

/* multiply a vector by transpose of a matrix from the left side. */
zVec zMulMatTVec(const zMat m, const zVec v1, zVec v)
{
  if( !zMatRowVecSizeEqual( m, v1 ) || !zMatColVecSizeEqual( m, v ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  return zMulMatTVecNC( m, v1, v );
}

/* multiply two matrices. */
zMat zMulMatMat(const zMat m1, const zMat m2, zMat m)
{
  if( !zMatColRowSizeEqual( m1, m2 ) ||
      !zMatRowSizeEqual( m1, m ) || !zMatColSizeEqual( m2, m ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  return zMulMatMatNC( m1, m2, m );
}

/* multiply a matrix by transpose of another matrix from the right side. */
zMat zMulMatMatT(const zMat m1, const zMat m2, zMat m)
{
  if( !zMatColSizeEqual( m1, m2 ) || !zMatRowSizeEqual( m1, m ) ||
      !zMatRowColSizeEqual( m2, m ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  return zMulMatMatTNC( m1, m2, m );
}

/* multiply a matrix by transpose of another matrix from the left side. */
zMat zMulMatTMat(const zMat m1, const zMat m2, zMat m)
{
  if( !zMatRowSizeEqual( m1, m2 ) ||
      !zMatColRowSizeEqual( m1, m ) || !zMatColSizeEqual( m2, m ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  return zMulMatTMatNC( m1, m2, m );
}

/* multiply a vector by a matrix directly.
 */
zVec zMulMatVecDRC(const zMat m, zVec v)
{
  zVec tmp;

  if( !( tmp = zVecAlloc( zVecSizeNC(v) ) ) ) return NULL;
  zMulMatVec( m, v, tmp );
  zVecCopyNC( tmp, v );
  zVecFree( tmp );
  return v;
}

/* multiply a vector by transpose of a matrix directly. */
zVec zMulMatTVecDRC(const zMat m, zVec v)
{
  zVec tmp;

  if( !( tmp = zVecAlloc( zVecSizeNC(v) ) ) ) return NULL;
  zMulMatTVec( m, v, tmp );
  zVecCopyNC( tmp, v );
  zVecFree( tmp );
  return v;
}

/* quadratic multiplication of matrices ('q = a diag{w} a^T') without checking size consistency. */
zMat zMatQuadNC(const zMat a, const zVec w, zMat q)
{
  int i, j, k;
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
zMat zMatQuad(const zMat a, const zVec w, zMat q)
{
  if( w && !zMatColVecSizeEqual( a, w ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  if( !zMatRowSizeEqual( a, q ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  if( !zMatIsSqr( q ) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return NULL;
  }
  return zMatQuadNC( a, w, q );
}

/* quadratic multiplication of matrices ('q = a^T diag{w} a') without checking size consistency. */
zMat zMatTQuadNC(const zMat a, const zVec w, zMat q)
{
  int i, j, k;
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
zMat zMatTQuad(const zMat a, const zVec w, zMat q)
{
  if( w && !zMatRowVecSizeEqual( a, w ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  if( !zMatColRowSizeEqual( a, q ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  if( !zMatIsSqr( q ) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return NULL;
  }
  return zMatTQuadNC( a, w, q );
}

/* quadratic multiplication of matrices ('m = a q a^T') without checking size consistency. */
zMat zMulMatMatMatTNC(const zMat a, const zMat q, zMat m)
{
  int i, j, k;
  double y;

  zMatZero( m );
  for( i=0; i<zMatRowSizeNC(q); i++ )
    for( j=0; j<zMatRowSizeNC(a); j++ ){
      y = 0;
      for( k=0; k<zMatColSizeNC(q); k++ )
        y += zMatElemNC(q,i,k) * zMatElemNC(a,j,k);
      for( k=0; k<zMatRowSizeNC(a); k++ )
        zMatElemNC(m,k,j) += zMatElemNC(a,k,i) * y;
    }
  return m;
}

/* quadratic multiplication of matrices ('m = a q a^T'). */
zMat zMulMatMatMatT(const zMat a, const zMat q, zMat m)
{
  if( !zMatIsSqr( q ) || !zMatIsSqr( m ) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return NULL;
  }
  if( !zMatColRowSizeEqual( a, q ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  return zMulMatMatMatTNC( a, q, m );
}

/* quadratic multiplication of matrices ('m = a^T q a') without checking size consistency. */
zMat zMulMatTMatMatNC(const zMat a, const zMat q, zMat m)
{
  int i, j, k;
  double y;

  zMatZero( m );
  for( i=0; i<zMatRowSizeNC(q); i++ )
    for( j=0; j<zMatColSizeNC(a); j++ ){
      y = 0;
      for( k=0; k<zMatColSizeNC(q); k++ )
        y += zMatElemNC(q,i,k) * zMatElemNC(a,k,j);
      for( k=0; k<zMatColSizeNC(a); k++ )
        zMatElemNC(m,k,j) += zMatElemNC(a,i,k) * y;
    }
  return m;
}

/* quadratic multiplication of matrices ('m = a^T q a'). */
zMat zMulMatTMatMat(const zMat a, const zMat q, zMat m)
{
  if( !zMatIsSqr( q ) || !zMatIsSqr( m ) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return NULL;
  }
  if( !zMatRowSizeEqual( a, q ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  return zMulMatTMatMatNC( a, q, m );
}

/* ********************************************************** */
/* I/O
 * ********************************************************** */

/* read a matrix from a ZTK format processor. */
zMat zMatFromZTK(ZTK *ztk)
{
  int i, j, row, col;
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
  int i, j;
  int row, col;
  zMat m;

  if( !zFInt( fp, &row ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZENOTFOUND );
    return NULL;
  }
  if( !zFInt( fp, &col ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZENOTFOUND );
    return NULL;
  }
  if( !( m = zMatAlloc( row, col ) ) ) return NULL;
  for( i=0; i<row; i++ )
    for( j=0; j<col; j++ )
      if( !zFDouble( fp, &zMatElemNC(m,i,j) ) ){
        ZRUNERROR( ZM_WARN_MAT_SIZEMISMATCH, i*row+j, row*col );
        break;
      }
  return m;
}

/* print information of a matrix to file. */
void zMatFPrint(FILE *fp, const zMat m)
{
  int i, j;

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
void zMatImg(const zMat m)
{
  double d;
  const char pat_pos[] = ",x*M";
  const char pat_neg[] = ".oO@";
  int i, j;
  int c;

  d = 0.25 * zMatElemAbsMax( m, NULL );
  for( i=0; i<zMatRowSizeNC(m); i++ ){
    for( j=0; j<zMatColSizeNC(m); j++ ){
      if( zIsTiny( zMatElemNC(m,i,j) ) )
        printf( "  " );
      else{
        c = zRound( zMatElemNC(m,i,j) / d );
        c = _zLimit( c, -3, 3 );
        printf( "%c ", c < 0 ? pat_neg[-c] : pat_pos[c] );
      }
    }
    zEndl();
  }
}
