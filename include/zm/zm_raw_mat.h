/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_raw_mat - raw vector and matrix : matrix.
 */

#ifndef __ZM_RAW_MAT_H__
#define __ZM_RAW_MAT_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief zero a raw matrix.
 *
 * zRawMatZero() sets all components of a raw matrix \a m for zeros.
 *
 * zRawMatTouchup() replaces all components less than zTOL (defined in zm_misc.h) of a raw matrix
 * \a m for zeros.
 *
 * \a r and \a c are the row and column sizes of the matrices, respectively.
 * \return
 * zRawMatZero() returns a pointer \a m.
 * zRawMatTouchup() returns no value.
 */
#define zRawMatZero(m,r,c)    zRawVecZero( m, (r)*(c) )
#define zRawMatTouchup(m,r,c) zRawVecTouchup( m, (r)*(c) )

/*! \brief copy a raw matrix.
 *
 * zRawMatCopy() copies a raw matrix \a src to another \a dest. \a r and \a c are the row and
 * column sizes of the matrices, respectively.
 * \return
 * zRawMatCopy() returns a pointer \a dest.
 */
#define zRawMatCopy(src,dest,r,c) zRawVecCopy( src, dest, (r)*(c) )

/*! \brief make identity, diagonal and random raw matrices.
 *
 * zRawMatIdent() creates an identity matrix \a m.
 *
 * zRawMatDiag() creates a diagonal matrix \a m, referring a raw-value array \a d as the diagonal values.
 * \a m does not need to be a squared matrix.
 *
 * zRawMatRand() sets all components of \a m randomly within the range from \a min to \a max.
 * \return
 * zRawMatIdent(), zRawMatDiag()
 */
__ZM_EXPORT void zRawMatIdent(double *m, int row, int col);
__ZM_EXPORT void zRawMatDiag(double *m, int row, int col, double *d);
#define zRawMatRandUniform(m,r,c,min,max) zRawVecRandUniform( m, (r)*(c), min, max )
#define zRawMatRand(m,min,max,r,c) zRawVecRand( m, min, max, (r)*(c) )

/*! \brief partially copy a matrix.
 *
 * zRawMatGet() gets a submatrix of \a src from (\a pr, \a pc) to (\a pr + \a dr - 1, \a pc + \a dc - 1)
 * and puts it into another matrix \a dest.
 * The size of \a dest must be \a dr x \a dc.
 * The size of \a src is specified as \a sr x \a sc.
 *
 * zRawMatTGet() gets a submatrix of \a src from (\a pr, \a pc) to (\a pr + \a dc - 1, \a pc + \a dr - 1)
 * and puts its transpose into another matrix \a dest.
 * The size of \a dest must be \a dr x \a dc.
 * The size of \a src is specified as \a sr x \a sc.
 *
 * zRawMatPut() puts a \a sr x \a sc matrix \a src into another matrix \a dest from (\a pr, \a pc).
 * The size of \a dest is specified as \a dr x \a dc.
 *
 * zRawMatTPut() puts the transpose of a \a sr x \a sc matrix \a src into another matrix \a dest
 * from (\a pr, \a pc). The size of \a dest is specified as \a dr x \a dc.
 * \return
 * These functions do not retur any values.
 */
__ZM_EXPORT void zRawMatGet(const double *src, int sr, int sc, int pr, int pc, double *dest, int dr, int dc);
__ZM_EXPORT void zRawMatPut(double *dest, int dr, int dc, int pr, int pc, const double *src, int sr, int sc);
__ZM_EXPORT void zRawMatTGet(const double *src, int sr, int sc, int pr, int pc, double *dest, int dr, int dc);
__ZM_EXPORT void zRawMatTPut(double *dest, int dr, int dc, int pr, int pc, const double *src, int sr, int sc);

/*! \brief abstract, set and swap row/column vectors from a raw matrix.
 *
 * zRawMatGetRow() abstracts the \a sr'th row vector of a matrix \a m. zRawMatGetCol() abstracts
 * the \a sc th column vector of a matrix \a m. Both functions put the result into \a v.
 * zRawMatPutRow() sets a vector \a v to the \a dr th row of a matrix \a m.
 * zRawMatPutCol() sets a vector \a v to the \a dc th column of a matrix \a m.
 *
 * zRawMatSwapRow() swaps the \a r1 th row and the \a r2 th row of a matrix \a m. zRawMatSwapCol()
 * swaps the \a c1 th column and \a c2 th column of \a m. \a row and \a col are the row and column
 * size of \a m, respectively.
 * \return
 * These functions do not retur any values.
 */
__ZM_EXPORT void zRawMatGetRow(const double *m, int row, int col, int sr, double *v);
__ZM_EXPORT void zRawMatGetCol(const double *m, int row, int col, int sc, double *v);
__ZM_EXPORT void zRawMatPutRow(double *m, int row, int col, int dr, const double *v);
__ZM_EXPORT void zRawMatPutCol(double *m, int row, int col, int dc, const double *v);
__ZM_EXPORT void zRawMatSwapRow(double *m, int row, int col, int r1, int r2);
__ZM_EXPORT void zRawMatSwapCol(double *m, int row, int col, int c1, int c2);

/*! \brief check if a matrix is tiny.
 *
 * zRawMatIsTol() returns the true value if all elements of a matrix \a m are smaller than \a e,
 * or the false value, otherwise.
 * zRawMatIsTiny() returns the true value if all elements of \a m are smaller than zTOL
 * (defined in zm_misc.h), or the false value, otherwise.
 *
 * \a row and \a col are the row and column size of \a m, respectively.
 * \return
 * zRawMatIsTol() and zRawMatIsTiny() return results as boolean values.
 */
#define zRawMatIsTol(m,r,c,t) zRawVecIsTol( m, (r)*(c), t )
#define zRawMatIsTiny(m,r,c)  zRawMatIsTol( m, r, c, zTOL )

/*! \brief basic arithmetics for matrix.
 *
 * zRawMatAdd() adds two matrices \a m1 and \a m2.
 * zRawMatSub() subtracts \a m2 from \a m1.
 * zRawMatRev() reverses the sign of \a m1.
 * zRawMatMul() multiplies \a m1 by a scalar value \a k.
 * zRawMatDiv() divides \a m1 by \a k.
 * zRawMatCat() concatenates \a m1 with \a m2 multiplied by \a k.
 * For all the above functions, the result is put into \a m.
 *
 * zRawMatAddDRC() adds \a m2 to \a m1 directly.
 * zRawMatSubDRC() subtracts \a m2 from \a m1 directly.
 * zRawMatRevDRC() reverses the sign of \a m1 directly.
 * zRawMatMulDRC() multiplies \a m1 by a scalar value \a k directly.
 * zRawMatDivDRC() divides \a m1 by \a k directly.
 * zRawMatCatDRC() concatenates \a m1 with \a m2 multiplied by \a k directly.
 *
 * \a row and \a col are the row and column size of the matrices, respectively.
 * \return
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

/*! \brief calculate the norm of a raw matrix.
 *
 * \return
 * zRawMatSqrNorm() returns the squared norm of a raw matrix \a m.
 * zRawMatNorm() returns the norm of \a m.
 *
 * \a row and \a col are the row and column size of \a m, respectively.
 */
#define zRawMatSqrNorm(m,r,c) zRawVecSqrNorm( m, (r)*(c) )
#define zRawMatNorm(m,r,c)    zRawVecNorm( m, (r)*(c) )

/*! \brief transpose a raw matrix.
 *
 * zRawMatT() transposes a raw matrix \a m. The result is put into \a tm.
 * The size of \a tm must be \a row x \a col.
 *
 * zRawMatTDRC() directly transposes \a m. The size of \a m must be \a row x \a col.
 *
 * zRawMatTr() calculates the trace value of \a m. The size of \a m must be \a row x \a col.
 * \return
 * zRawMatT() and zRawMatTDRC() return no values.
 * zRawMatTr() returns the value calculated.
 */
__ZM_EXPORT void zRawMatT(const double *m, double *tm, int row, int col);
__ZM_EXPORT void zRawMatTDRC(double *m, int row, int col);
__ZM_EXPORT double zRawMatTr(const double *m, int row, int col);

/*! \brief multiply a raw vector by a raw matrix.
 *
 * zRawMulMatVec() multiplies a raw vector \a v1 by a raw matrix \a m.
 * zRawMulMatTVec() multiplies a raw vector \a v1 by the transpose of a raw matrix \a m.
 * For the both functions, the result is put into \a v.
 *
 * zRawMulMatMat() multiplies a raw matrix \a m1 by another raw matrix \a m2 from the right side.
 * The sizes of \a m1 and \a m2 must be \a r1 x \a c1 and \a c1 x \a c2, respectively.
 * zRawMulMatMatT() multiplies \a m1 by the transpose of \a m2 from the right side. The sizes of
 * \a m1 and \a m2 must be \a r1 x \a c1 and \a r2 x \a c1, respectively.
 * zRawMulMatTMat() multiplies \a m2 by the transpose of \a m1 from the left side. The sizes of
 * \a m1 and \a m2 are \a r1 x \a c1 and \a r1 x \a c2.
 * For those functions, the result is put into \a m.
 * \return
 * These functions return no values.
 */
__ZM_EXPORT void zRawMulMatVec(const double *m, const double *v1, int row, int col, double *v);
__ZM_EXPORT void zRawMulMatTVec(const double *m, const double *v1, int row, int col, double *v);
__ZM_EXPORT void zRawMulMatMat(const double *m1, int r1, int c1, const double *m2, int r2, int c2, double *m);
__ZM_EXPORT void zRawMulMatMatT(const double *m1, int r1, int c1, const double *m2, int r2, int c2, double *m);
__ZM_EXPORT void zRawMulMatTMat(const double *m1, int r1, int c1, const double *m2, int r2, int c2, double *m);

/*! \brief dyadic product of raw vectors.
 *
 * zRawVecDyad() calculates the dyadic product of two raw vectors \a v1 and \a v2. Namely, \a v1 \a v2^T.
 * The result is put into \a dyad.
 *
 * zRawMatAddDyad() adds the dyadic products of \a v1 and \a v2 to a raw matrix \a m.
 * zRawMatSubDyad() subtracts the dyadic products of \a v1 and \a v2 from a raw matrix \a m.
 * zRawMatCatDyad() concatenates \a m with the dyadic product of \a v1 and \a v2 multiplied by a scalar
 * value \a k.
 * For those functions, \a size1 and \a size2 are the sizes of \a v1 and \a v2, respectively.
 * \return
 * These functions return no value.
 */
__ZM_EXPORT void zRawVecDyad(const double *v1, int size1, const double *v2, int size2, double *dyad);
__ZM_EXPORT void zRawMatAddDyad(double *m, const double *v1, int size1, const double *v2, int size2);
__ZM_EXPORT void zRawMatSubDyad(double *m, const double *v1, int size1, const double *v2, int size2);
__ZM_EXPORT void zRawMatCatDyad(double *m, double k, const double *v1, int size1, const double *v2, int size2);

/*! \brief print a raw matrix.
 *
 * zRawMatFPrint() prints a raw matrix \a m to the current position of a file \a fp. The size of
 * \a m is specified as \a row x \a col.
 *
 * zRawMatPrint() prints \a m out to the standard output.
 * \return
 * These functions return no value.
 */
__ZM_EXPORT void zRawMatFPrint(FILE *fp, const double *m, int row, int col);
#define zRawMatPrint(m,r,c) zRawMatFPrint( stdout, m, r, c )

__END_DECLS

#endif /* __ZM_RAW_MAT_H__ */
