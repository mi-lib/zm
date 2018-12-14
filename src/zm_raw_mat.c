/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_raw_mat - raw vector and matrix : matrix.
 */

#include <zm/zm_raw.h>

/* zRawMatIdent
 * - create a raw identity matrix.
 */
void zRawMatIdent(double *m, int row, int col)
{
  zRawMatClear( m, row, col );
  for( ; row>0; row--, m+=col+1 ) *m = 1;
}

/* zRawMatDiag
 * - create a raw diagonal matrix.
 */
void zRawMatDiag(double *m, int row, int col, double *d)
{
  zRawMatClear( m, row, col );
  for( ; row>0; row--, m+=col+1, d++ ) *m = *d;
}

/* zRawMatGet
 * - get submatrix.
 */
void zRawMatGet(double *src, int sr, int sc, int pr, int pc, double *dest, int dr, int dc)
{
  src += pr*sc;
  for( ; dr>0; src+=sc, dest+=dc, dr-- )
    zRawVecGet( src, pc, dest, dc );
}

/* zRawMatTGet
 * - get transpose of a submatrix.
 */
void zRawMatTGet(double *src, int sr, int sc, int pr, int pc, double *dest, int dr, int dc)
{
  for( src+=pr*sc+pc; dr>0; src++, dest+=dc, dr-- )
    zRawMatGetCol( src, dc, sc, 0, dest );
}

/* zRawMatPut
 * - put submatrix.
 */
void zRawMatPut(double *dest, int dr, int dc, int pr, int pc, double *src, int sr, int sc)
{
  dest += pr*dc;
  for( ; sr>0; dest+=dc, src+=sc, sr-- )
    zRawVecPut( dest, pc, src, sc );
}

/* zRawMatTPut
 * - put transpose of a submatrix.
 */
void zRawMatTPut(double *dest, int dr, int dc, int pr, int pc, double *src, int sr, int sc)
{
  dest += pr*dc + pc;
  for( ; sr>0; src+=sc, dest++, sr-- )
    zRawMatSetCol( dest, sc, dc, 0, src );
}

/* zRawMatGetRow
 * - abstruction of a row vector from a raw matrix.
 */
void zRawMatGetRow(double *m, int row, int col, int sr, double *v)
{
  zRawVecGet( m, col*sr, v, col );
}

/* zRawMatGetCol
 * - abstruction of a column vector from a raw matrix.
 */
void zRawMatGetCol(double *m, int row, int col, int sc, double *v)
{
  for( m+=sc; row>0; row--, v++, m+=col ) *v = *m;
}

/* zRawMatSetRow
 * - set of a row vector into a raw matrix.
 */
void zRawMatSetRow(double *m, int row, int col, int dr, double *v)
{
  zRawVecPut( m, col*dr, v, col );
}

/* zRawMatSetCol
 * - set of a column vector into a raw matrix.
 */
void zRawMatSetCol(double *m, int row, int col, int dc, double *v)
{
  for( m+=dc; row>0; row--, v++, m+=col ) *m = *v;
}

void zRawMatSwapRow(double *m, int row, int col, int r1, int r2)
{
  register int i, _r1, _r2;

  if( r1 < r2 ){
    _r1 = r1; _r2 = r2;
  } else{
    _r1 = r2; _r2 = r1;
  }
  for( i=(_r2-_r1)*col, m+=col*_r1; col>0; col--, m++ )
    zRawVecSwap( m, 0, i );
}

void zRawMatSwapCol(double *m, int row, int col, int c1, int c2)
{
  register int i, _c1, _c2;

  if( c1 < c2 ){
    _c1 = c1; _c2 = c2;
  } else{
    _c1 = c2; _c2 = c1;
  }
  for( i=_c2-_c1, m+=_c1; row>0; row--, m+=col )
    zRawVecSwap( m, 0, i );
}

/* zRawMatT
 * - transpose of a raw array matrix.
 */
void zRawMatT(double *m, double *tm, int row, int col)
{
  register int i;

  /* tm : ('row' x 'col') */
  for( i=0; i<row; i++, tm+=col )
    zRawMatGetCol( m, col, row, i, tm );
}

/* zRawMatTDST
 * - transpose of a raw array matrix (destructive).
 */
void zRawMatTDST(double *mat, int row, int col)
{
  register int r, c;
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

/* zRawMatTr
 * - calculation of a trace value of a raw array as a matrix.
 */
double zRawMatTr(double *m, int row, int col)
{
  double result = 0;
  register int n;

  n = zMin( row, col );
  for( ; n>0 ; n--, m+=col+1 ) result += *m;
  return result;
}

/* zRawMulMatVec
 * - multiplication of raw arrays as a matrix and a vector.
 */
void zRawMulMatVec(double *m, double *v1, int row, int col, double *v)
{
  for( ; row>0; row--, v++, m+=col )
    *v = zRawVecInnerProd( m, v1, col );
}

/* zRawMulVecMat
 * - multiplication of raw arrays as a row vector and a matrix.
 */
void zRawMulVecMat(double *v1, double *m, int row, int col, double *v)
{
  register int i, j;
  double *vp, *mp;

  zRawVecClear( v, col );
  for( i=0; i<col; i++, v++, m++ )
    for( vp=v1, mp=m, j=0; j<row; j++, vp++, mp+=col )
      *v += (*mp) * (*vp);
}

/* zRawMulMatMat
 * - multiplication of raw arrays as matrices.
 */
void zRawMulMatMat(double *m1, int r1, int c1, double *m2, int r2, int c2, double *m)
{
  register int i;

  for( i=0; i<r1; i++, m1+=c1, m+=c2 )
    zRawMulVecMat( m1, m2, c1, c2, m );
}

/* zRawMulMatMatT
 * - multiplication of raw arrays as a matrix and transpose of a matrix.
 */
void zRawMulMatMatT(double *m1, int r1, int c1, double *m2, int r2, int c2, double *m)
{
  register int i;
  double *mp;

  for( ; r1>0; r1--, m1+=c1 )
    for( mp=m2, i=0; i<r2; i++, mp+=c1 )
      *m++ = zRawVecInnerProd( m1, mp, c1 );
}

/* zRawMulMatMatT
 * - multiplication of raw arrays as transpose of a matrix and a matrix.
 */
void zRawMulMatTMat(double *m1, int r1, int c1, double *m2, int r2, int c2, double *m)
{
  register int i, j, k;
  double *mp, *mp1, *mp2;

  zRawMatClear( m, c1, c2 );
  for( i=0; i<c1; i++, m1++ )
    for( mp=m2, j=0; j<c2; j++, mp++, m++ )
      for( mp1=m1, mp2=mp, k=0; k<r1; k++, mp1+=c1, mp2+=c2 )
        *m += (*mp1) * (*mp2);
}

/* zRawVecDyad
 * - dyad of raw vectors.
 */
void zRawVecDyad(double *v1, int size1, double *v2, int size2, double *dyad)
{
  register int i, j;
  double *vp;

  for( i=0; i<size1; i++, v1++ )
    for( j=0, vp=v2; j<size2; j++, vp++ )
      *dyad++ = *v1 * *vp;
}

/* zRawMatAddDyad
 * - add dyad of raw vectors to matrix.
 */
void zRawMatAddDyad(double *m, double *v1, int size1, double *v2, int size2)
{
  register int i, j;
  double *vp;

  for( i=0; i<size1; i++, v1++ )
    for( j=0, vp=v2; j<size2; j++, vp++ )
      *m++ += *v1 * *vp;
}

/* zRawMatSubDyad
 * - subtract dyad of raw vectors to matrix.
 */
void zRawMatSubDyad(double *m, double *v1, int size1, double *v2, int size2)
{
  register int i, j;
  double *vp;

  for( i=0; i<size1; i++, v1++ )
    for( j=0, vp=v2; j<size2; j++, vp++ )
      *m++ -= *v1 * *vp;
}

/* zRawMatCatDyad
 * - add scalar-multiplied dyad of raw vectors to matrix.
 */
void zRawMatCatDyad(double *m, double k, double *v1, int size1, double *v2, int size2)
{
  register int i, j;
  double *vp;

  for( i=0; i<size1; i++, v1++ )
    for( j=0, vp=v2; j<size2; j++, vp++ )
      *m++ += k * *v1 * *vp;
}

/* zRawMatFWrite
 * - output raw matrix.
 */
void zRawMatFWrite(FILE *fp, double *m, int row, int col)
{
  for( ; row>0; row--, m+=col )
    zRawVecFWrite( fp, m, col );
}
