/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_raw_mat - raw vector and matrix : matrix.
 */

#ifndef __ZM_RAW_MAT_H__
#define __ZM_RAW_MAT_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* METHOD:
 * zRawMatClear, zRawMatTouchup, zRawMatCopy
 * - clear, touchup and copy raw matrix.
 *
 * 'zRawMatClear()' clears a matrix 'm', setting all
 * components for zeros.
 * #
 * 'zRawMatTouchup()' touches up 'm', namely, replace
 * all components which are less than zTOL for zeros.
 * #
 * 'zRawMatCopy()' copies a matrix 'src' to 'dest'.
 * 'row' and 'col' are the row and column size of matrices,
 * respectively.
 * #
 * 'row' and 'col' are the row and column size of 'm',
 * respectively.
 * [RETURN VALUE]
 * Neither 'zRawMatClear()' nor 'zRawMatTouchup()' return
 * any values.
 * #
 * 'zRawMatCopy()' returns a pointer to 'dest'.
 */
#define zRawMatClear(m,r,c)      zRawVecClear( m, (r)*(c) )
#define zRawMatTouchup(m,r,c)    zRawVecTouchup( m, (r)*(c) )
#define zRawMatCopy(src,dst,r,c) zRawVecCopy( src, dst, (r)*(c) )

/* METHOD:
 * zRawMatIdent, zRawMatDiag, zRawMatRand
 * - identity, diagonal and random raw matrices.
 *
 * 'zRawMatIdent()' creates an identity matrix 'm'.
 * #
 * 'zRawMatDiag()' creates a diagonal matrix 'm', referring
 * a raw-value array 'd' as the diagonal values.
 * 'm' does not have to be a squared matrix.
 * #
 * 'zRawMatRand()' sets all the components of 'm' randomly
 * within the range from 'min' to 'max'.
 * #
 * 'row' and 'col' are the row and column size of 'm',
 * respectively.
 * [RETURN VALUE]
 * These functions return no values.
 */
__EXPORT void zRawMatIdent(double *m, int row, int col);
__EXPORT void zRawMatDiag(double *m, int row, int col, double *d);
#define zRawMatRandUniform(m,r,c,min,max) zRawVecRandUniform( m, (r)*(c), min, max )
#define zRawMatRand(m,min,max,r,c) zRawVecRand( m, min, max, (r)*(c) )

/* METHOD:
 * zRawMatGet, zRawMatPut
 * - partial copy of matrices.
 *
 * 'zRawMatGet()' gets a submatrix of 'src' with the
 * size of 'dr' x 'dc' from ('pr', 'pc) to 'dest'.
 * 'dr' and 'dc' is also the size of 'dest'. The size
 * of 'src' is 'sr' x 'sc'.
 * 'zRawMatPut()' puts 'dest' to 'src' as a submatrix
 * with the size of 'sr' x 'sc' at ('pr', 'pc).
 * 'dr' and 'dc' is the size of 'dest'.
 * [RETURN VALUES]
 * These functions return no values.
 */
__EXPORT void zRawMatGet(double *src, int sr, int sc, int pr, int pc, double *dest, int dr, int dc);
__EXPORT void zRawMatPut(double *dest, int dr, int dc, int pr, int pc, double *src, int sr, int sc);
__EXPORT void zRawMatTGet(double *src, int sr, int sc, int pr, int pc, double *dest, int dr, int dc);
__EXPORT void zRawMatTPut(double *dest, int dr, int dc, int pr, int pc, double *src, int sr, int sc);

/* METHOD:
 * zRawMatGetRow, zRawMatGetCol, zRawMatSetRow, zRawMatSetCol,
 * zRawMatSwapRow, zRawMatSwapCol
 * - abstraction, set and swap of row/column vector from raw matrix.
 *
 * 'zRawMatGetRow()' abstracts a row vector at 'sr' from a
 * matrix 'm'. 'zRawMatGetCol()' abstracts a column vector
 * at 'sc' from 'm'. Both put the result into 'v'.
 * 'zRawMatSetRow()' sets a row vector 'v' to 'm' at 'dr'.
 * 'zRawMatSetCol()' sets a column vector 'v' to 'm'
 * at 'dc'.
 * #
 * 'zRawMatSwapRow()' swaps 'r1'th row and 'r2'th row of
 * 'm'. 'zRawMatSwapCol()' swaps 'c1'th column and 'c2'th
 * column of 'm'.
 * #
 * 'row' and 'col' are the row and column size of matrices,
 * respectively.
 * [RETURN VALUE]
 * These functions return no values.
 * [NOTES]
 * 'zRawMatSwapRow()' and 'zRawMatSwapCol()' are only
 * available in the user space.
 */
__EXPORT void zRawMatGetRow(double *m, int row, int col, int sr, double *v);
__EXPORT void zRawMatGetCol(double *m, int row, int col, int sc, double *v);
__EXPORT void zRawMatSetRow(double *m, int row, int col, int dr, double *v);
__EXPORT void zRawMatSetCol(double *m, int row, int col, int dc, double *v);
__EXPORT void zRawMatSwapRow(double *m, int row, int col, int r1, int r2);
__EXPORT void zRawMatSwapCol(double *m, int row, int col, int c1, int c2);

/* METHOD:
 * zRawMatIsTol, zRawMatIsTiny - see if matrix is tiny.
 *
 * 'zRawMatIsTol()' returns the true value if every components
 * of 'm' are smaller than 'e', or the false value, otherwise.
 * 'zRawMatIsTiny()' compares each component of 'm' with
 * zTOL(defined in "zm_misc.h"), and returns the same result
 * with 'zRawMatIsTol()'.
 * #
 * 'row' and 'col' are the row and column size of matrices,
 * respectively.
 * [RETURN VALUE]
 * 'zRawMatIsTol()', 'zRawMatIsTiny()' return results
 * as boolean values.
 */
#define zRawMatIsTol(m,r,c,t) zRawVecIsTol( m, (r)*(c), t )
#define zRawMatIsTiny(m,r,c)  zRawMatIsTol( m, r, c, zTOL )

/* METHOD:
 * zRawMatAdd, zRawMatSub, zRawMatRev,
 * zRawMatMul, zRawMatDiv, zRawMatCat,
 * zRawMatAddDRC, zRawMatSubDRC, zRawMatRevDRC,
 * zRawMatMulDRC, zRawMatDivDRC, zRawMatCatDRC
 * - basic arithmetics for matrix.
 *
 * 'zRawMatAdd()' adds two matrices 'm1' and 'm2'.
 * 'zRawMatSub()' subtracts 'm2' from 'm1'.
 * 'zRawMatRev()' reverses the sign of 'm1'.
 * 'zRawMatMul()' multiplies 'm1' by a scalar value 'k'.
 * 'zRawMatDiv()' divides 'm1' by 'k'.
 * 'zRawMatCat()' concatenates 'm1', adding multiplied
 * 'm2' by 'k'.
 * These functions put the result into 'm'.
 * #
 * 'zRawMatAddDRC()' adds 'm2' to 'm1' directly.
 * 'zRawMatSubDRC()' subtracts 'm2' from 'm1' directly.
 * 'zRawMatRevDRC()' reverses the sign of 'm' directly.
 * 'zRawMatMulDRC()' multiplies 'm' by a scalar value
 * 'k' directly.
 * 'zRawMatDivDRC()' divides 'm' by 'k' directly.
 * 'zRawMatCatDRC()' concatenates 'm1', adding multiplied
 * 'm2' by 'k'.
 * #
 * 'row' and 'col' are the row and column size of matrices,
 * respectively.
 * [RETURN VALUE]
 * These functions return no values.
 */
#define zRawMatAdd(m1,m2,m,r,c)    zRawVecAdd( m1, m2, m, (r)*(c) )
#define zRawMatSub(m1,m2,m,r,c)    zRawVecSub( m1, m2, m, (r)*(c) )
#define zRawMatRev(m1,m,r,c)       zRawVecRev( m1, m, (r)*(c) )
#define zRawMatMul(m1,k,m,r,c)     zRawVecMul( m1, k, m, (r)*(c) )
#define zRawMatDiv(m1,k,m,r,c)     zRawVecDiv( m1, k, m, (r)*(c) )
#define zRawMatCat(m1,k,m2,m,r,c)  zRawVecCat( m1, k, m2, m, (r)*(c) )

#define zRawMatAddDRC(m1,m2,r,c)   zRawVecAddDRC( m1, m2, (r)*(c) )
#define zRawMatSubDRC(m1,m2,r,c)   zRawVecSubDRC( m1, m2, (r)*(c) )
#define zRawMatRevDRC(m,r,c)       zRawVecRevDRC( m, (r)*(c) )
#define zRawMatMulDRC(m,k,r,c)     zRawVecMulDRC( m, k, (r)*(c) )
#define zRawMatDivDRC(m,k,r,c)     zRawVecDivDRC( m, k, (r)*(c) )
#define zRawMatCatDRC(m1,k,m2,r,c) zRawVecCatDRC( m1, k, m2, (r)*(c) )

/* METHOD:
 * zRawMatSqrNorm, zRawMatNorm - calculation of matrix norm.
 *
 * 'zRawMatSqrNorm()' calculates the squared norm of 'm',
 * and 'zRawMatNorm()' calculates the norm of 'm'.
 * #
 * 'row' and 'col' are the row and column size of matrices,
 * respectively.
 * [RETURN VALUE]
 * These functions return the norm calculated.
 */
#define zRawMatSqrNorm(m,r,c) zRawVecSqrNorm( m, (r)*(c) )
#define zRawMatNorm(m,r,c)    zRawVecNorm( m, (r)*(c) )

/* METHOD:
 * zRawMatT, zRawMatTDST, zRawMatTr
 * - transpose of matrix.
 *
 * 'zRawMatT()' creates a transpose matrix of 'm'
 * and puts it into 'tm'. 'row' and 'col' are of 'tm'.
 * #
 * 'zRawMatTDST()' destructively modifies 'm' to
 * the transpose of itself. 'row' and 'col' are of 'm'.
 * #
 * 'zRawMatTr()' calculates the trace value of a
 * matrix 'm'.
 * #
 * 'row' and 'col' are the row and column size of matrices,
 * respectively.
 * [RETURN VALUES]
 *  'zRawMatT()' and 'zRawMatTDST()' return no values.
 */
__EXPORT void zRawMatT(double *m, double *tm, int row, int col);
__EXPORT void zRawMatTDST(double *m, int row, int col);
__EXPORT double  zRawMatTr(double *m, int row, int col);

/* METHOD:
 * zRawMulMatVec, zRawMulVecMat,
 * zRawMulMatMat, zRawMulMatMatT, zRawMulMatTMat
 * - multiplication of double-precision floating-point value
 *   vector and matrix.
 *
 * These multiplication functions are prepared for allowing
 * one to treat raw arrays as vectors or matrices.
 * #
 * 'zRawMulMatVec()' multiplies a matrix 'm' by a vector
 * 'v1'. 'zRawMulVecMat()' multiplies a row vector 'v1'
 * by 'm'. These two functions put the result into 'v'.
 * #
 * 'zRawMulMatMat()' multiplies a matrix 'm1' by another
 * matrix 'm2'. 'm1' is 'r1' x 'c1' matrix, and 'm2' is
 * 'c1' x 'c2' matrix.
 * 'zRawMulMatMatT()' multiplies 'm1' by the
 * transpose of 'm2'. 'm1' is 'r1' x 'c1' matrix, and 'm2' is
 * 'r2' x 'c1' matrix.
 * 'zRawMulMatTMat()' multiplies the transpose
 * of 'm1' by 'm2'. 'm1' is 'r1' x 'c1' matrix, and 'm2' is
 * 'r1' x 'c2' matrix.
 * These three functions put the result into 'm'.
 * [RETURN VALUE]
 * These functions return no values.
 */
__EXPORT void zRawMulMatVec(double *m, double *v1, int row, int col, double *v);
__EXPORT void zRawMulVecMat(double *v1, double *m, int row, int col, double *v);
__EXPORT void zRawMulMatMat(double *m1, int r1, int c1, double *m2, int r2, int c2, double *m);
__EXPORT void zRawMulMatMatT(double *m1, int r1, int c1, double *m2, int r2, int c2, double *m);
__EXPORT void zRawMulMatTMat(double *m1, int r1, int c1, double *m2, int r2, int c2, double *m);

/* METHOD:
 * zRawVecDyad, zRawMatAddDyad, zRawMatSubDyad, zRawMatCatDyad
 * - dyadic product of double-precision floating-point value vectors.
 *
 * 'zRawVecDyad()' calculates dyadic product of
 * two vector 'v1' and 'v2', namely, 'v1' 'v2'^T.
 * The result will be put into 'dyad'.
 * #
 * 'zRawMatAddDyad()' adds dyadic products of 'v1'
 * and 'v2' to matrix 'm'.
 * 'zRawMatSubDyad()' subtracts dyadic products of 'v1'
 * and 'v2' to matrix 'm'.
 * 'zRawMatCatDyad()' adds a multiplied dyadic
 * products of 'v1' and 'v2' by 'k' to 'm'.
 * #
 * For those three functions, 'size1' and 'size2'
 * are the size of 'v1' and 'v2', respectively.
 * [RETURN VALUE]
 * 'zRawVecDyad()', 'zRawMatAddDyad()', 'zRawMatSubDyad()'
 * and 'zRawMatCatDyad()' return no value.
 */
__EXPORT void zRawVecDyad(double *v1, int size1, double *v2, int size2, double *dyad);
__EXPORT void zRawMatAddDyad(double *m, double *v1, int size1, double *v2, int size2);
__EXPORT void zRawMatSubDyad(double *m, double *v1, int size1, double *v2, int size2);
__EXPORT void zRawMatCatDyad(double *m, double k, double *v1, int size1, double *v2, int size2);

/* METHOD:
 * zRawMatFWrite, zRawMatWrite
 * - output raw matrix.
 *
 * 'zRawMatFWrite()' outputs all components of a given 2-D raw
 * array of floating-point values pointed  by 'm' with a row size
 * 'row' and a column size 'col' to the current position of file
 * 'fp'.
 * 'zRawMatWrite()' outputs all components of 'm' simply to
 * the standard output.
 * [RETURN VALUE]
 * Neither 'zRawMatFWrite()' nor 'zRawMatWrite()' return
 * any values.
 */
__EXPORT void zRawMatFWrite(FILE *fp, double *m, int row, int col);
#define zRawMatWrite(m,r,c) zRawMatFWrite( stdout, m, r, c )

__END_DECLS

#endif /* __ZM_RAW_MAT_H__ */
