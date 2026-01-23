/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_raw_mat - raw vector and matrix : matrix.
 */

#include <zm/zm_raw.h>

/* zero a raw matrix. */
void zRawMatZero(double *m, int colcapacity, int rowsize, int colsize)
{
  for( ; --rowsize>=0; m+=colcapacity ) zRawVecZero( m, colsize );
}

/* touch-up a raw matrix. */
void zRawMatTouchup(double *m, int colcapacity, int rowsize, int colsize, double tol)
{
  for( ; --rowsize>=0; m+=colcapacity ) zRawVecTouchup( m, colsize, tol );
}

/* copy a raw matrix. */
void zRawMatCopy(const double *src, int srccolcapacity, double *dest, int destcolcapacity, int rowsize, int colsize)
{
  for( ; --rowsize>=0; src+=srccolcapacity, dest+=destcolcapacity ) zRawVecCopy( src, dest, colsize );
}

/* check if two raw matrices are equal. */
bool zRawMatEqual(const double *m1, int colcapacity1, const double *m2, int colcapacity2, int rowsize, int colsize, double tol)
{
  for( ; --rowsize>=0; m1+=colcapacity1, m2+=colcapacity2 )
    if( !zRawVecEqual( m1, m2, colsize, tol ) ) return false;
  return true;
}

/* check if two raw matrices exactly match with each other. */
bool zRawMatMatch(const double *m1, int colcapacity1, const double *m2, int colcapacity2, int rowsize, int colsize)
{
  for( ; --rowsize>=0; m1+=colcapacity1, m2+=colcapacity2 )
    if( !zRawVecMatch( m1, m2, colsize ) ) return false;
  return true;
}

/* check if all components of a raw matrix are smaller than a tolerance. */
bool zRawMatIsTol(const double *m, int colcapacity, int rowsize, int colsize, double tol)
{
  for( ; --rowsize>=0; m+=colcapacity )
    if( !zRawVecIsTol( m, colsize, tol ) ) return false;
  return true;
}

/* create a raw identity matrix. */
void zRawMatIdent(double *m, int colcapacity, int rowsize, int colsize)
{
  int size;

  zRawMatZero( m, colcapacity, rowsize, colsize );
  size = _zMin( rowsize, colsize );
  for( ; --size>=0; m+=colcapacity+1 ) *m = 1;
}

/* create a raw diagonal matrix. */
void zRawMatDiag(double *m, int colcapacity, double *d, int rowsize, int colsize)
{
  int size;

  zRawMatZero( m, colcapacity, rowsize, colsize );
  size = _zMin( rowsize, colsize );
  for( ; --size>=0; m+=colcapacity+1, d++ ) *m = *d;
}

/* create a random raw matrix. */
void zRawMatRandUniform(double *m, int colcapacity, double min, double max, int rowsize, int colsize)
{
  for( ; --rowsize>=0; m+=colcapacity ) zRawVecRandUniform( m, colsize, min, max );
}

/* create a random raw matrix. */
void zRawMatRand(double *m, int colcapacity, double *matmin, int matmincolcapacity, double *matmax, int matmaxcolcapacity, int rowsize, int colsize)
{
  for( ; --rowsize>=0; m+=colcapacity, matmin+=matmincolcapacity, matmax+=matmaxcolcapacity )
    zRawVecRand( m, matmin, matmax, colsize );
}

/* abstract a row vector from a raw matrix. */
void zRawMatGetRow(const double *m, int colcapacity, int rowsize, int colsize, int row, double *rowvec)
{
  zRawVecGet( m+(row*colcapacity), 0, rowvec, colsize );
}

/* abstract a column vector from a raw matrix. */
void zRawMatGetCol(const double *m, int colcapacity, int rowsize, int colsize, int col, double *colvec)
{
  for( m+=col; --rowsize>=0; m+=colcapacity, colvec++ ) *colvec = *m;
}

/* put a row vector into a raw matrix. */
void zRawMatPutRow(double *m, int colcapacity, int rowsize, int colsize, int row, const double *rowvec)
{
  zRawVecPut( m+(row*colcapacity), 0, rowvec, colsize );
}

/* put a column vector into a raw matrix. */
void zRawMatPutCol(double *m, int colcapacity, int rowsize, int colsize, int col, const double *colvec)
{
  for( m+=col; --rowsize>=0; m+=colcapacity, colvec++ ) *m = *colvec;
}

/* swap two row vectors in a raw matrix. */
void zRawMatSwapRow(double *m, int colcapacity, int rowsize, int colsize, int row1, int row2)
{
  double *m1, *m2, tmp;

  m1 = m + row1 * colcapacity;
  m2 = m + row2 * colcapacity;
  for( ; --colsize>=0; m1++, m2++ ){
    tmp = *m1;
    *m1 = *m2;
    *m2 = tmp;
  }
}

/* swap two column vectors in a raw matrix. */
void zRawMatSwapCol(double *m, int colcapacity, int rowsize, int colsize, int col1, int col2)
{
  double tmp;

  for( ; --rowsize>=0; m+=colcapacity ){
    tmp = *( m + col1 );
    *( m + col1 ) = *( m + col2 );
    *( m + col2 ) = tmp;
  }
}

/* get a submatrix from a raw matrix. */
void zRawMatGet(const double *src, int srccolcapacity, int srcrowsize, int srccolsize, int r, int c, double *dest, int destcolcapacity, int destrowsize, int destcolsize)
{
  if( r + destrowsize > srcrowsize ){
    ZRUNWARN( ZM_WARN_MAT_DEST_ROWOVERSIZED, destrowsize, destcolsize, srcrowsize, srccolsize, r, c );
    destrowsize = srcrowsize - r;
  }
  if( c + destcolsize > srccolsize ){
    ZRUNWARN( ZM_WARN_MAT_DEST_COLUMNOVERSIZED, destrowsize, destcolsize, srcrowsize, srccolsize, r, c );
    destcolsize = srccolsize - c;
  }
  for( src+=r*srccolcapacity+c; --destrowsize>=0; src+=srccolcapacity, dest+=destcolcapacity )
    zRawVecGet( src, 0, dest, destcolsize );
}

/* put a submatrix into another raw matrix. */
void zRawMatPut(double *dest, int destcolcapacity, int destrowsize, int destcolsize, int r, int c, const double *src, int srccolcapacity, int srcrowsize, int srccolsize)
{
  if( r + srcrowsize > destrowsize ){
    ZRUNWARN( ZM_WARN_MAT_SRC_ROWOVERSIZED, srcrowsize, srccolsize, destrowsize, destcolsize, r, c );
    srcrowsize = destrowsize - r;
  }
  if( c + srccolsize > destcolsize ){
    ZRUNWARN( ZM_WARN_MAT_SRC_COLUMNOVERSIZED, srcrowsize, srccolsize, destrowsize, destcolsize, r, c );
    srccolsize = destcolsize - c;
  }
  for( dest+=r*destcolcapacity+c; --srcrowsize>=0; dest+=destcolcapacity, src+=srccolcapacity )
    zRawVecPut( dest, 0, src, srccolsize );
}

/* get transpose of a submatrix from a raw matrix. */
void zRawMatTGet(const double *src, int srccolcapacity, int srcrowsize, int srccolsize, int r, int c, double *dest, int destcolcapacity, int destrowsize, int destcolsize)
{
  if( r + destcolsize > srcrowsize ){
    ZRUNWARN( ZM_WARN_MAT_DEST_ROWOVERSIZED, destrowsize, destcolsize, srcrowsize, srccolsize, r, c );
    destcolsize = srcrowsize - r;
  }
  if( c + destrowsize > srccolsize ){
    ZRUNWARN( ZM_WARN_MAT_DEST_COLUMNOVERSIZED, destrowsize, destcolsize, srcrowsize, srccolsize, r, c );
    destrowsize = srccolsize - c;
  }
  for( src+=r*srccolcapacity+c; --destrowsize>=0; src++, dest+=destcolcapacity )
    zRawMatGetCol( src, srccolcapacity, destcolsize, destrowsize, 0, dest );
}

/* put transpose of a submatrix into another matrix. */
void zRawMatTPut(double *dest, int destcolcapacity, int destrowsize, int destcolsize, int r, int c, const double *src, int srccolcapacity, int srcrowsize, int srccolsize)
{
  if( r + srccolsize > destrowsize ){
    ZRUNWARN( ZM_WARN_MAT_SRC_ROWOVERSIZED, srcrowsize, srccolsize, destrowsize, destcolsize, r, c );
    srccolsize = destrowsize - r;
  }
  if( c + srcrowsize > destcolsize ){
    ZRUNWARN( ZM_WARN_MAT_SRC_COLUMNOVERSIZED, srcrowsize, srccolsize, destrowsize, destcolsize, r, c );
    srcrowsize = destcolsize - c;
  }
  for( dest+=r*destcolcapacity+c; --srcrowsize>=0; src+=srccolcapacity, dest++ )
    zRawMatPutCol( dest, destcolcapacity, srccolsize, destcolsize, 0, src );
}

/* add two raw matrices. */
void zRawMatAdd(const double *m1, int colcapacity1, const double *m2, int colcapacity2, double *m, int colcapacity3, int rowsize, int colsize)
{
  for( ; --rowsize>=0; m1+=colcapacity1, m2+=colcapacity2, m+=colcapacity3 )
    zRawVecAdd( m1, m2, m, colsize );
}

/* subtract a raw matrix from another. */
void zRawMatSub(const double *m1, int colcapacity1, const double *m2, int colcapacity2, double *m, int colcapacity3, int rowsize, int colsize)
{
  for( ; --rowsize>=0; m1+=colcapacity1, m2+=colcapacity2, m+=colcapacity3 )
    zRawVecSub( m1, m2, m, colsize );
}

/* reverse a raw matrix. */
void zRawMatRev(const double *m1, int colcapacity1, double *m, int colcapacity2, int rowsize, int colsize)
{
  for( ; --rowsize>=0; m1+=colcapacity1, m+=colcapacity2 )
    zRawVecRev( m1, m, colsize );
}

/* multiply a raw matrix by a real number. */
void zRawMatMul(const double *m1, int colcapacity1, double k, double *m, int colcapacity2, int rowsize, int colsize)
{
  for( ; --rowsize>=0; m1+=colcapacity1, m+=colcapacity2 )
    zRawVecMul( m1, k, m, colsize );
}

/* divide a raw matrix by a real number. */
void zRawMatDiv(const double *m1, int colcapacity1, double k, double *m, int colcapacity2, int rowsize, int colsize)
{
  zRawMatMul( m1, colcapacity1, 1.0/k, m, colcapacity2, rowsize, colsize );
}

/* concatenate a raw matrix by another multiplied by a real number. */
void zRawMatCat(const double *m1, int colcapacity1, double k, const double *m2, int colcapacity2, double *m, int colcapacity3, int rowsize, int colsize)
{
  for( ; --rowsize>=0; m1+=colcapacity1, m2+=colcapacity2, m+=colcapacity3 )
    zRawVecCat( m1, k, m2, m, colsize );
}

/* add a raw matrix directly to another. */
void zRawMatAddDRC(double *m1, int colcapacity1, const double *m2, int colcapacity2, int rowsize, int colsize)
{
  for( ; --rowsize>=0; m1+=colcapacity1, m2+=colcapacity2 )
    zRawVecAddDRC( m1, m2, colsize );
}

/* subtract a raw matrix directly from another. */
void zRawMatSubDRC(double *m1, int colcapacity1, const double *m2, int colcapacity2, int rowsize, int colsize)
{
  for( ; --rowsize>=0; m1+=colcapacity1, m2+=colcapacity2 )
    zRawVecSubDRC( m1, m2, colsize );
}

/* reverse a raw matrix directly. */
void zRawMatRevDRC(double *m, int colcapacity, int rowsize, int colsize)
{
  for( ; --rowsize>=0; m+=colcapacity )
    zRawVecRevDRC( m, colsize );
}

/* multiply a raw matrix directly by a real number. */
void zRawMatMulDRC(double *m, int colcapacity, double k, int rowsize, int colsize)
{
  for( ; --rowsize>=0; m+=colcapacity )
    zRawVecMulDRC( m, k, colsize );
}

/* divide a raw matrix directly by a real number. */
void zRawMatDivDRC(double *m, int colcapacity, double k, int rowsize, int colsize)
{
  zRawMatMulDRC( m, colcapacity, 1.0/k, rowsize, colsize );
}

/* concatenate a raw matrix multiplied by a real number directly. */
void zRawMatCatDRC(double *m1, int colcapacity1, double k, const double *m2, int colcapacity2, int rowsize, int colsize)
{
  for( ; --rowsize>=0; m1+=colcapacity1, m2+=colcapacity2 )
    zRawVecCatDRC( m1, k, m2, colsize );
}

/* add a raw vector to a column vector of a raw matrix directly. */
void zRawMatColAddDRC(double *m, int colcapacity, const double *colvec, int rowsize, int colsize, int col)
{
  for( m+=col; --rowsize>=0; m+=colcapacity, colvec++ ) *m += *colvec;
}

/* subtract a raw vector from a column vector of a raw matrix directly. */
void zRawMatColSubDRC(double *m, int colcapacity, const double *colvec, int rowsize, int colsize, int col)
{
  for( m+=col; --rowsize>=0; m+=colcapacity, colvec++ ) *m -= *colvec;
}

/* multiply a column vector of a raw matrix by a scalar value directly. */
void zRawMatColMulDRC(double *m, int colcapacity, double k, int rowsize, int colsize, int col)
{
  for( m+=col; --rowsize>=0; m+=colcapacity ) *m *= k;
}

/* concatenate a raw vector multiplied by a scalar value to a column vector of a raw matrix directly. */
void zRawMatColCatDRC(double *m, int colcapacity, double k, const double *colvec, int rowsize, int colsize, int col)
{
  for( m+=col; --rowsize>=0; m+=colcapacity, colvec++ ) *m += k * *colvec;
}

/* inner product of a column vector of a raw matrix and another raw vector. */
double zRawMatColInnerProd(const double *m, int colcapacity, const double *v, int rowsize, int colsize, int col)
{
  double s=0, s_prev=0, c, q=0, r;

  for( m+=col; --rowsize>=0; m+=colcapacity ){
    s = s_prev + ( c = *m * *v++ );
    r = s - s_prev;
    q += c - r;
    s_prev = s;
  }
  return s + q;
}

/* squared norm of a raw matrix. */
double zRawMatSqrNorm(const double *m, int colcapacity, int rowsize, int colsize)
{
  double norm = 0;

  for( ; --rowsize>=0; m+=colcapacity )
    norm += zRawVecSqrNorm( m, colsize );
  return norm;
}

/* transpose a raw matrix. */
void zRawMatT(const double *m, int colcapacity1, double *tm, int colcapacity2, int rowsize, int colsize)
{
  for( ; --colsize>=0; m++, tm+=colcapacity2 )
    zRawMatGetCol( m, colcapacity1, rowsize, colsize, 0, tm );
}

/* find the previous index to move a component for in-place transpose of a raw matrix. */
static int _zRawMatTDRCPrevIndex(int index, int rowcapacity, int colcapacity)
{
  return ( index % rowcapacity ) * colcapacity + ( index / rowcapacity );
}

/* directly transpose a raw matrix in-place. */
void zRawMatTDRC(double *mat, int rowcapacity, int colcapacity)
{
  int i, size;
  int index, prev_index;
  double tmp;

  if( ( size = rowcapacity * colcapacity ) <= 1 ) return;
  for( i=0; i<size; i++ ){
    prev_index = _zRawMatTDRCPrevIndex( i, rowcapacity, colcapacity );
    while( prev_index > i )
      prev_index = _zRawMatTDRCPrevIndex( prev_index, rowcapacity, colcapacity );
    if( prev_index != i ) continue;
    tmp = mat[i];
    prev_index = _zRawMatTDRCPrevIndex( ( index = i ), rowcapacity, colcapacity );
    while( prev_index != i ){
      mat[index] = mat[prev_index];
      prev_index = _zRawMatTDRCPrevIndex( ( index = prev_index ), rowcapacity, colcapacity );
    }
    mat[index] = tmp;
  }
}

/* calculate the trace of a raw matrix. */
double zRawMatTrace(const double *m, int colcapacity, int rowsize, int colsize)
{
  double trace = 0, trace_prev = 0, q = 0, r;
  int n;

  n = _zMin( rowsize, colsize );
  for( ; --n>=0; m+=colcapacity+1 ){
    trace = trace_prev + *m;
    r = trace - trace_prev;
    q += *m - r;
    trace_prev = trace;
  }
  return trace + q;
}

/* multiply a raw vector by a raw matrix. */
void zRawMulMatVec(const double *m, int colcapacity, const double *v, int rowsize, int colsize, double *mv)
{
  for( ; --rowsize>=0; m+=colcapacity, mv++ )
    *mv = zRawVecInnerProd( m, v, colsize );
}

/* multiply a raw vector by transpose of a raw matrix. */
void zRawMulMatTVec(const double *m, int colcapacity, const double *v, int rowsize, int colsize, double *mv)
{
  int i;

  for( i=0; i<colsize; i++ )
    *( mv + i ) = zRawMatColInnerProd( m, colcapacity, v, rowsize, colsize, i );
}

/* multiply a raw matrix by another. */
void zRawMulMatMat(const double *m1, int colcapacity1, int rowsize1, int colsize1, const double *m2, int colcapacity2, int rowsize2, int colsize2, double *m, int colcapacity3)
{
  /* rowsize2 must be equal to colsize1. */
  for( ; --rowsize1>=0; m1+=colcapacity1, m+=colcapacity3 )
    zRawMulMatTVec( m2, colcapacity2, m1, rowsize2, colsize2, m );
}

/* multiply a raw matrix by transpose of another. */
void zRawMulMatMatT(const double *m1, int colcapacity1, int rowsize1, int colsize1, const double *m2, int colcapacity2, int rowsize2, int colsize2, double *m, int colcapacity3)
{
  /* colsize2 must be equal to colsize1. */
  for( ; --rowsize1>=0; m1+=colcapacity1, m+=colcapacity3 )
    zRawMulMatVec( m2, colcapacity2, m1, rowsize2, colsize2, m );
}

/* multiply transpose of a raw matrix by another raw matrix. */
void zRawMulMatTMat(const double *m1, int colcapacity1, int rowsize1, int colsize1, const double *m2, int colcapacity2, int rowsize2, int colsize2, double *m, int colcapacity3)
{
  const double *mp1, *mp2;
  double *mp;
  int i, j;

  /* rowsize2 must be equal to rowsize1. */
  for( ; --colsize1>=0; m1++, m+=colcapacity3 )
    for( i=0, mp=m; i<colsize2; i++, mp++ )
      for( j=0, mp1=m1, mp2=m2+i, *mp=0; j<rowsize1; j++, mp1+=colcapacity1, mp2+=colcapacity2 )
        *mp += *mp1 * *mp2;
}

/* dyadic product of two raw vectors. */
void zRawVecDyad(const double *v1, int size1, const double *v2, int size2, double *dyad, int colcapacity)
{
  const double *vp;
  double *mp;
  int i;

  for( ; --size1>=0; v1++, dyad+=colcapacity )
    for( mp=dyad, vp=v2, i=0; i<size2; i++, vp++, mp++ )
      *mp = *v1 * *vp;
}

/* add dyadic product of two raw vectors to a raw matrix. */
void zRawMatAddDyad(double *m, int colcapacity, const double *v1, int size1, const double *v2, int size2)
{
  const double *vp;
  double *mp;
  int i;

  for( ; --size1>=0; v1++, m+=colcapacity )
    for( mp=m, vp=v2, i=0; i<size2; i++, vp++, mp++ )
      *mp += *v1 * *vp;
}

/* subtract dyadic product of two raw vectors from a raw matrix. */
void zRawMatSubDyad(double *m, int colcapacity, const double *v1, int size1, const double *v2, int size2)
{
  const double *vp;
  double *mp;
  int i;

  for( ; --size1>=0; v1++, m+=colcapacity )
    for( mp=m, vp=v2, i=0; i<size2; i++, vp++, mp++ )
      *mp -= *v1 * *vp;
}

/* concatenate a raw matrix with dyadic product of two raw vectors multiplied by a scalar value. */
void zRawMatCatDyad(double *m, int colcapacity, double k, const double *v1, int size1, const double *v2, int size2)
{
  const double *vp;
  double *mp;
  int i;

  for( ; --size1>=0; v1++, m+=colcapacity )
    for( mp=m, vp=v2, i=0; i<size2; i++, vp++, mp++ )
      *mp += k * *v1 * *vp;
}

/* print a raw matrix out to a file. */
void zRawMatFPrint(FILE *fp, const double *m, int rowcapacity, int colcapacity, int rowsize, int colsize)
{
  fprintf( fp, "(%d", rowsize );
  if( rowcapacity != rowsize )
    fprintf( fp, "/%d", rowcapacity );
  fprintf( fp, ", %d", colsize );
  if( colcapacity != colsize )
    fprintf( fp, "/%d", colcapacity );
  fprintf( fp, ") {\n" );
  for( ; --rowsize>=0; m+=colcapacity )
    zRawVecValueFPrint( fp, m, colsize );
  fprintf( fp, "}\n" );
}
