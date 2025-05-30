/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_raw_mat - raw vector and matrix : matrix.
 */

#include <zm/zm_raw.h>

/* create a raw identity matrix. */
void zRawMatIdent(double *m, int row, int col)
{
  zRawMatZero( m, row, col );
  for( ; row>0; row--, m+=col+1 ) *m = 1;
}

/* create a raw diagonal matrix. */
void zRawMatDiag(double *m, int row, int col, double *d)
{
  zRawMatZero( m, row, col );
  for( ; row>0; row--, m+=col+1, d++ ) *m = *d;
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
void zRawMatGetRow(const double *m, int row, int col, int sr, double *v)
{
  zRawVecGet( m, col*sr, v, col );
}

/* abstract a column vector from a raw matrix. */
void zRawMatGetCol(const double *m, int row, int col, int sc, double *v)
{
  for( m+=sc; row>0; row--, v++, m+=col ) *v = *m;
}

/* put a row vector into a raw matrix. */
void zRawMatPutRow(double *m, int row, int col, int dr, const double *v)
{
  zRawVecPut( m, col*dr, v, col );
}

/* put a column vector into a raw matrix. */
void zRawMatPutCol(double *m, int row, int col, int dc, const double *v)
{
  for( m+=dc; row>0; row--, v++, m+=col ) *m = *v;
}

/* swap two row vectors in a raw matrix. */
void zRawMatSwapRow(double *m, int row, int col, int r1, int r2)
{
  int i, _r1, _r2;

  if( r1 < r2 ){
    _r1 = r1; _r2 = r2;
  } else{
    _r1 = r2; _r2 = r1;
  }
  for( i=(_r2-_r1)*col, m+=col*_r1; col>0; col--, m++ )
    zRawVecSwap( m, 0, i );
}

/* swap two column vectors in a raw matrix. */
void zRawMatSwapCol(double *m, int row, int col, int c1, int c2)
{
  int i, _c1, _c2;

  if( c1 < c2 ){
    _c1 = c1; _c2 = c2;
  } else{
    _c1 = c2; _c2 = c1;
  }
  for( i=_c2-_c1, m+=_c1; row>0; row--, m+=col )
    zRawVecSwap( m, 0, i );
}

/* transpose a raw matrix. */
void zRawMatT(const double *m, double *tm, int row, int col)
{
  int i;

  /* tm : ('row' x 'col') */
  for( i=0; i<row; i++, tm+=col )
    zRawMatGetCol( m, col, row, i, tm );
}

/* directly transpose a raw matrix. */
void zRawMatTDRC(double *mat, int row, int col)
{
  int r, c;
  long offset, offset_next, size;
  bool *check;
  double *mp, *cp, tmp;

  size = row * col;
  if( !( check = zAlloc( bool, size ) ) ){
    ZALLOCERROR();
    return;
  }
  for( mp=mat+1; mp<mat+size; mp++ ){
    tmp = *mp;
    if( check[(offset=mp-mat)] ) continue;
    for( cp=mp; ; cp=mat+(offset=offset_next) ){
      check[offset] = true;
      r = offset / row;
      c = offset - r * row;
      if( check[( offset_next = col * c + r /* transpose */)] ) break;
      *cp = mat[offset_next];
    }
    *cp = tmp;
  }
  free( check );
}

/* calculate the trace of a raw matrix. */
double zRawMatTrace(const double *m, int row, int col)
{
  double result = 0;
  int n;

  n = _zMin( row, col );
  for( ; n>0 ; n--, m+=col+1 ) result += *m;
  return result;
}

/* multiply a raw vector by a raw matrix. */
void zRawMulMatVec(const double *m, const double *v1, int row, int col, double *v)
{
  for( ; row>0; row--, v++, m+=col )
    *v = zRawVecInnerProd( m, v1, col );
}

/* multiply a raw vector by transpose of a raw matrix. */
void zRawMulMatTVec(const double *m, const double *v1, int row, int col, double *v)
{
  int i, j;
  const double *vp, *mp;

  for( i=0; i<col; i++, v++, m++ ){
    *v = 0;
    for( vp=v1, mp=m, j=0; j<row; j++, vp++, mp+=col )
      *v += (*mp) * (*vp);
  }
}

/* multiply a raw matrix by another. */
void zRawMulMatMat(const double *m1, int r1, int c1, const double *m2, int r2, int c2, double *m)
{
  int i;

  for( i=0; i<r1; i++, m1+=c1, m+=c2 )
    zRawMulMatTVec( m2, m1, c1, c2, m );
}

/* multiply a raw matrix by transpose of another. */
void zRawMulMatMatT(const double *m1, int r1, int c1, const double *m2, int r2, int c2, double *m)
{
  int i;
  const double *mp;

  for( ; r1>0; r1--, m1+=c1 )
    for( mp=m2, i=0; i<r2; i++, mp+=c1 )
      *m++ = zRawVecInnerProd( m1, mp, c1 );
}

/* multiply transpose of a raw matrix by another raw matrix. */
void zRawMulMatTMat(const double *m1, int r1, int c1, const double *m2, int r2, int c2, double *m)
{
  int i, j, k;
  const double *mp, *mp1, *mp2;

  for( i=0; i<c1; i++, m1++ )
    for( mp=m2, j=0; j<c2; j++, mp++, m++ ){
      *m = 0;
      for( mp1=m1, mp2=mp, k=0; k<r1; k++, mp1+=c1, mp2+=c2 )
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
void zRawMatFPrint(FILE *fp, const double *m, int row, int col)
{
  for( ; row>0; row--, m+=col )
    zRawVecFPrint( fp, m, col );
}
