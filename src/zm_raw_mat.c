/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_raw_mat - raw vector and matrix : matrix.
 */

#include <zm/zm_raw.h>

/* create a raw identity matrix. */
void zRawMatIdent(double *m, int rowsize, int colsize)
{
  int i, size;

  zRawMatZero( m, rowsize, colsize );
  size = _zMin( rowsize, colsize );
  for( i=0; i<size; i++, m+=colsize+1 ) *m = 1;
}

/* create a raw diagonal matrix. */
void zRawMatDiag(double *m, int rowsize, int colsize, double *d)
{
  int i, size;

  zRawMatZero( m, rowsize, colsize );
  size = _zMin( rowsize, colsize );
  for( i=0; i<size; i++, m+=colsize+1, d++ ) *m = *d;
}

/* get a submatrix from a raw matrix. */
void zRawMatGet(const double *src, int sr, int sc, int pr, int pc, double *dest, int dr, int dc)
{
  src += pr*sc;
  for( ; dr>0; src+=sc, dest+=dc, dr-- )
    zRawVecGet( src, pc, dest, dc );
}

/* get transpose of a submatrix from a raw matrix. */
void zRawMatTGet(const double *src, int sr, int sc, int pr, int pc, double *dest, int dr, int dc)
{
  for( src+=pr*sc+pc; dr>0; src++, dest+=dc, dr-- )
    zRawMatGetCol( src, dc, sc, 0, dest );
}

/* put a submatrix into another raw matrix. */
void zRawMatPut(double *dest, int dr, int dc, int pr, int pc, const double *src, int sr, int sc)
{
  dest += pr*dc;
  for( ; sr>0; dest+=dc, src+=sc, sr-- )
    zRawVecPut( dest, pc, src, sc );
}

/* put transpose of a submatrix into another matrix. */
void zRawMatTPut(double *dest, int dr, int dc, int pr, int pc, const double *src, int sr, int sc)
{
  dest += pr*dc + pc;
  for( ; sr>0; src+=sc, dest++, sr-- )
    zRawMatPutCol( dest, sc, dc, 0, src );
}

/* abstract a row vector from a raw matrix. */
void zRawMatGetRow(const double *m, int rowsize, int colsize, int sr, double *v)
{
  zRawVecGet( m, colsize*sr, v, colsize );
}

/* abstract a column vector from a raw matrix. */
void zRawMatGetCol(const double *m, int rowsize, int colsize, int sc, double *v)
{
  for( m+=sc; rowsize>0; rowsize--, v++, m+=colsize ) *v = *m;
}

/* put a row vector into a raw matrix. */
void zRawMatPutRow(double *m, int rowsize, int colsize, int dr, const double *v)
{
  zRawVecPut( m, colsize*dr, v, colsize );
}

/* put a column vector into a raw matrix. */
void zRawMatPutCol(double *m, int rowsize, int colsize, int dc, const double *v)
{
  for( m+=dc; rowsize>0; rowsize--, v++, m+=colsize ) *m = *v;
}

/* swap two row vectors in a raw matrix. */
void zRawMatSwapRow(double *m, int rowsize, int colsize, int r1, int r2)
{
  int i, _r1, _r2;

  if( r1 < r2 ){
    _r1 = r1; _r2 = r2;
  } else{
    _r1 = r2; _r2 = r1;
  }
  for( i=(_r2-_r1)*colsize, m+=colsize*_r1; colsize>0; colsize--, m++ )
    zRawVecSwap( m, 0, i );
}

/* swap two column vectors in a raw matrix. */
void zRawMatSwapCol(double *m, int rowsize, int colsize, int c1, int c2)
{
  int i, _c1, _c2;

  if( c1 < c2 ){
    _c1 = c1; _c2 = c2;
  } else{
    _c1 = c2; _c2 = c1;
  }
  for( i=_c2-_c1, m+=_c1; rowsize>0; rowsize--, m+=colsize )
    zRawVecSwap( m, 0, i );
}

/* add a raw vector to a column vector of a raw matrix directly. */
double *zRawMatColAddDRC(double *m, const double *colvec, int rowsize, int colsize, int col)
{
  for( m+=col; rowsize>0; rowsize--, m+=colsize, colvec++ )
    *m += *colvec;
  return m;
}

/* subtract a raw vector from a column vector of a raw matrix directly. */
double *zRawMatColSubDRC(double *m, const double *colvec, int rowsize, int colsize, int col)
{
  for( m+=col; rowsize>0; rowsize--, m+=colsize, colvec++ )
    *m -= *colvec;
  return m;
}

/* multiply a column vector of a raw matrix by a scalar value directly. */
double *zRawMatColMulDRC(double *m, double k, int rowsize, int colsize, int col)
{
  for( m+=col; rowsize>0; rowsize--, m+=colsize )
    *m *= k;
  return m;
}

/* concatenate a raw vector multiplied by a scalar value to a column vector of a raw matrix directly. */
double *zRawMatColCatDRC(double *m, double k, const double *colvec, int rowsize, int colsize, int col)
{
  for( m+=col; rowsize>0; rowsize--, m+=colsize, colvec++ )
    *m += k * *colvec;
  return m;
}

/* inner product of a column vector of a raw matrix and another raw vector. */
double zRawMatColInnerProd(const double *m, const double *v, int rowsize, int colsize, int col)
{
  double s=0, s_prev=0, c, q=0, r;

  for( m+=col; rowsize-->0; m+=colsize ){
    s = s_prev + ( c = *m * *v++ );
    r = s - s_prev;
    q += c - r;
    s_prev = s;
  }
  return s + q;
}

/* transpose a raw matrix. */
void zRawMatT(const double *m, double *tm, int rowsize, int colsize)
{
  int i;

  /* tm : ('rowsize' x 'colsize') */
  for( i=0; i<rowsize; i++, tm+=colsize )
    zRawMatGetCol( m, colsize, rowsize, i, tm );
}

/* directly transpose a raw matrix. */
void zRawMatTDRC(double *mat, int rowsize, int colsize)
{
  int r, c;
  long offset, offset_next, size;
  bool *check;
  double *mp, *cp, tmp;

  size = rowsize * colsize;
  if( !( check = zAlloc( bool, size ) ) ){
    ZALLOCERROR();
    return;
  }
  for( mp=mat+1; mp<mat+size; mp++ ){
    tmp = *mp;
    if( check[(offset=mp-mat)] ) continue;
    for( cp=mp; ; cp=mat+(offset=offset_next) ){
      check[offset] = true;
      r = offset / rowsize;
      c = offset - r * rowsize;
      if( check[( offset_next = colsize * c + r /* transpose */)] ) break;
      *cp = mat[offset_next];
    }
    *cp = tmp;
  }
  free( check );
}

/* calculate the trace of a raw matrix. */
double zRawMatTrace(const double *m, int rowsize, int colsize)
{
  double result = 0;
  int n;

  n = _zMin( rowsize, colsize );
  for( ; n>0 ; n--, m+=colsize+1 ) result += *m;
  return result;
}

/* multiply a raw vector by a raw matrix. */
void zRawMulMatVec(const double *m, const double *v1, int rowsize, int colsize, double *v)
{
  for( ; rowsize>0; rowsize--, v++, m+=colsize )
    *v = zRawVecInnerProd( m, v1, colsize );
}

/* multiply a raw vector by transpose of a raw matrix. */
void zRawMulMatTVec(const double *m, const double *v1, int rowsize, int colsize, double *v)
{
  int i, j;
  const double *vp, *mp;

  for( i=0; i<colsize; i++, v++, m++ ){
    *v = 0;
    for( vp=v1, mp=m, j=0; j<rowsize; j++, vp++, mp+=colsize )
      *v += (*mp) * (*vp);
  }
}

/* multiply a raw matrix by another. */
void zRawMulMatMat(const double *m1, int rowsize1, int colsize1, const double *m2, int rowsize2, int colsize2, double *m)
{
  int i;

  for( i=0; i<rowsize1; i++, m1+=colsize1, m+=colsize2 )
    zRawMulMatTVec( m2, m1, colsize1, colsize2, m );
}

/* multiply a raw matrix by transpose of another. */
void zRawMulMatMatT(const double *m1, int rowsize1, int colsize1, const double *m2, int rowsize2, int colsize2, double *m)
{
  int i;
  const double *mp;

  for( ; rowsize1>0; rowsize1--, m1+=colsize1 )
    for( mp=m2, i=0; i<rowsize2; i++, mp+=colsize1 )
      *m++ = zRawVecInnerProd( m1, mp, colsize1 );
}

/* multiply transpose of a raw matrix by another raw matrix. */
void zRawMulMatTMat(const double *m1, int rowsize1, int colsize1, const double *m2, int rowsize2, int colsize2, double *m)
{
  int i, j, k;
  const double *mp, *mp1, *mp2;

  for( i=0; i<colsize1; i++, m1++ )
    for( mp=m2, j=0; j<colsize2; j++, mp++, m++ ){
      *m = 0;
      for( mp1=m1, mp2=mp, k=0; k<rowsize1; k++, mp1+=colsize1, mp2+=colsize2 )
        *m += (*mp1) * (*mp2);
    }
}

/* dyadic product of two raw vectors. */
void zRawVecDyad(const double *v1, int size1, const double *v2, int size2, double *dyad)
{
  int i, j;
  const double *vp;

  for( i=0; i<size1; i++, v1++ )
    for( j=0, vp=v2; j<size2; j++, vp++ )
      *dyad++ = *v1 * *vp;
}

/* add dyadic product of two raw vectors to a raw matrix. */
void zRawMatAddDyad(double *m, const double *v1, int size1, const double *v2, int size2)
{
  int i, j;
  const double *vp;

  for( i=0; i<size1; i++, v1++ )
    for( j=0, vp=v2; j<size2; j++, vp++ )
      *m++ += *v1 * *vp;
}

/* subtract dyadic product of two raw vectors from a raw matrix. */
void zRawMatSubDyad(double *m, const double *v1, int size1, const double *v2, int size2)
{
  int i, j;
  const double *vp;

  for( i=0; i<size1; i++, v1++ )
    for( j=0, vp=v2; j<size2; j++, vp++ )
      *m++ -= *v1 * *vp;
}

/* concatenate a raw matrix with dyadic product of two raw vectors multiplied by a scalar value. */
void zRawMatCatDyad(double *m, double k, const double *v1, int size1, const double *v2, int size2)
{
  int i, j;
  const double *vp;

  for( i=0; i<size1; i++, v1++ )
    for( j=0, vp=v2; j<size2; j++, vp++ )
      *m++ += k * *v1 * *vp;
}

/* print a raw matrix out to a file. */
void zRawMatFPrint(FILE *fp, const double *m, int rowsize, int colsize)
{
  for( ; rowsize>0; rowsize--, m+=colsize )
    zRawVecFPrint( fp, m, colsize );
}
