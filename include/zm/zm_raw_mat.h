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
 * zRawMatTouchup() replaces all components less than \a tol of a raw matrix
 * \a m for zeros.
 *
 * \a rowsize and \a colsize are the row and column sizes of the matrices, respectively.
 * \return
 * zRawMatZero() returns a pointer \a m.
 * zRawMatTouchup() returns no value.
 */
#define zRawMatZero(m,rowsize,colsize)        zRawVecZero( m, (rowsize)*(colsize) )
#define zRawMatTouchup(m,rowsize,colsize,tol) zRawVecTouchup( m, (rowsize)*(colsize), tol )

/*! \brief copy a raw matrix.
 *
 * zRawMatCopy() copies a raw matrix \a src to another \a dest. \a rowsize and \a colsize are the row and
 * column sizes of the matrices, respectively.
 * \return
 * zRawMatCopy() returns a pointer \a dest.
 */
#define zRawMatCopy(src,dest,rowsize,colsize) zRawVecCopy( src, dest, (rowsize)*(colsize) )

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
__ZM_EXPORT void zRawMatIdent(double *m, int rowsize, int colsize);
__ZM_EXPORT void zRawMatDiag(double *m, int rowsize, int colsize, double *d);
#define zRawMatRandUniform(m,rowsize,colsize,min,max) zRawVecRandUniform( m, (rowsize)*(colsize), min, max )
#define zRawMatRand(m,min,max,rowsize,colsize) zRawVecRand( m, min, max, (rowsize)*(colsize) )

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
 * swaps the \a c1 th column and \a c2 th column of \a m. \a rowsize and \a colsize are the row and column
 * size of \a m, respectively.
 * \return
 * These functions do not retur any values.
 */
__ZM_EXPORT void zRawMatGetRow(const double *m, int rowsize, int colsize, int sr, double *v);
__ZM_EXPORT void zRawMatGetCol(const double *m, int rowsize, int colsize, int sc, double *v);
__ZM_EXPORT void zRawMatPutRow(double *m, int rowsize, int colsize, int dr, const double *v);
__ZM_EXPORT void zRawMatPutCol(double *m, int rowsize, int colsize, int dc, const double *v);
__ZM_EXPORT void zRawMatSwapRow(double *m, int rowsize, int colsize, int r1, int r2);
__ZM_EXPORT void zRawMatSwapCol(double *m, int rowsize, int colsize, int c1, int c2);

/*! \brief check if two raw matrices are equal.
 *
 * zRawMatEqual() checks if two raw matrices \a m1 and \a m2 are equal. \a tol is the tolerance to regard
 * two values are the same. \a rowsize and \a colsize are the row and column sizes of the two matrices.
 * \return
 * zRawMatEqual() returns the true value if differences of all corresponding components of the two
 * matrices \a m1 and \a m2 are less than or equal to \a tol. Otherwise, it returns the false value.
 */
#define zRawMatEqual(m1,m2,rowsize,colsize,tol) zRawVecEqual( m1, m2, (rowsize)*(colsize), tol )
/*! \brief check if two matrices exactly match with each other.
 *
 * zRawMatMatch() checks if two raw matrices \a m1 and \a m2 match with each other.
 * \a rowsize and \a colsize are the row and column sizes of the two matrices.
 * \return
 * zRawMatMatch() returns the true value if all corresponding components of the two matrices \a m1 and
 * \a m2 match with each other. Otherwise, it returns the false value.
 */
#define zRawMatMatch(m1,m2,rowsize,colsize) zRawVecMatch( m1, m2, (rowsize)*(colsize) )

/*! \brief check if a matrix is tiny.
 *
 * zRawMatIsTol() returns the true value if all elements of a matrix \a m are smaller than \a e,
 * or the false value, otherwise.
 * zRawMatIsTiny() returns the true value if all elements of \a m are smaller than zTOL
 * (defined in zm_misc.h), or the false value, otherwise.
 *
 * \a rowsize and \a colsize are the row and column size of \a m, respectively.
 * \return
 * zRawMatIsTol() and zRawMatIsTiny() return results as boolean values.
 */
#define zRawMatIsTol(m,rowsize,colsize,t) zRawVecIsTol( m, (rowsize)*(colsize), t )
#define zRawMatIsTiny(m,rowsize,colsize)  zRawMatIsTol( m, rowsize, colsize, zTOL )

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
 * \a rowsize and \a colsize are the row and column size of the matrices, respectively.
 * \return
 * These functions return no values.
 */
#define zRawMatAdd(m1,m2,m,rowsize,colsize)    zRawVecAdd( m1, m2, m, (rowsize)*(colsize) )
#define zRawMatSub(m1,m2,m,rowsize,colsize)    zRawVecSub( m1, m2, m, (rowsize)*(colsize) )
#define zRawMatRev(m1,m,rowsize,colsize)       zRawVecRev( m1, m, (rowsize)*(colsize) )
#define zRawMatMul(m1,k,m,rowsize,colsize)     zRawVecMul( m1, k, m, (rowsize)*(colsize) )
#define zRawMatDiv(m1,k,m,rowsize,colsize)     zRawVecDiv( m1, k, m, (rowsize)*(colsize) )
#define zRawMatCat(m1,k,m2,m,rowsize,colsize)  zRawVecCat( m1, k, m2, m, (rowsize)*(colsize) )

#define zRawMatAddDRC(m1,m2,rowsize,colsize)   zRawVecAddDRC( m1, m2, (rowsize)*(colsize) )
#define zRawMatSubDRC(m1,m2,rowsize,colsize)   zRawVecSubDRC( m1, m2, (rowsize)*(colsize) )
#define zRawMatRevDRC(m,rowsize,colsize)       zRawVecRevDRC( m, (rowsize)*(colsize) )
#define zRawMatMulDRC(m,k,rowsize,colsize)     zRawVecMulDRC( m, k, (rowsize)*(colsize) )
#define zRawMatDivDRC(m,k,rowsize,colsize)     zRawVecDivDRC( m, k, (rowsize)*(colsize) )
#define zRawMatCatDRC(m1,k,m2,rowsize,colsize) zRawVecCatDRC( m1, k, m2, (rowsize)*(colsize) )

/*! \brief add a raw vector to a column vector of a raw matrix directly.
 *
 * zRawMatColAddDRC() adds a raw vector \a colvec to the \a col th column vector of a raw matrix \a m
 * directly.
 * \a rowsize and \a colsize are the row and column sizes of \a m, where \a rowsize is the same with
 * that of \a v.
 * \return
 * zRawMatColAddDRC() returns a pointer \a m.
 */
__ZM_EXPORT double *zRawMatColAddDRC(double *m, const double *colvec, int rowsize, int colsize, int col);
/*! \brief subtract a raw vector from a column vector of a raw matrix directly.
 *
 * zRawMatColSubDRC() subtracts a raw vector \a colvec from the \a col th column vector of a raw
 * matrix \a m directly.
 * \a rowsize and \a colsize are the row and column sizes of \a m, where \a rowsize is the same with
 * that of \a v.
 * \return
 * zRawMatColSubDRC() returns a pointer \a m.
 */
__ZM_EXPORT double *zRawMatColSubDRC(double *m, const double *colvec, int rowsize, int colsize, int col);
/*! \brief multiply a column vector of a raw matrix by a scalar value directly.
 *
 * zRawMatColMulDRC() multiplies the \a col th column vector of a raw matrix \a m by a scalar value
 * \a k directly.
 * \a rowsize and \a colsize are the row and column sizes of \a m.
 * \return
 * zRawMatColMulDRC() returns a pointer \a m.
 */
__ZM_EXPORT double *zRawMatColMulDRC(double *m, double k, int rowsize, int colsize, int col);
/*! \brief concatenate a raw vector multiplied by a scalar value to a column vector of a raw matrix directly.
 *
 * zRawMatColCatDRC() concatenates a raw vector \a colvec multiplied by a scalar value \a k to the \a col
 * th column vector of a raw matrix \a m directly.
 * \a rowsize and \a colsize are the row and column sizes of \a m, where \a rowsize is the same with
 * that of \a v.
 * \return
 * zRawMatColCatDRC() returns a pointer \a m.
 */
__ZM_EXPORT double *zRawMatColCatDRC(double *m, double k, const double *colvec, int rowsize, int colsize, int col);

/*! \brief inner product of a column vector of a raw matrix and another raw vector.
 *
 * zRawMatColInnerProd() calculates the inner product of the \a col th column vector of a raw matrix \a m
 * and another raw vector \a v.
 * \a rowsize and \a colsize are the row and column sizes of \a m, where \a rowsize is the same with
 * that of \a v.
 * \return
 * zRawMatColInnerProd() returns the calculated inner product.
 */
__ZM_EXPORT double zRawMatColInnerProd(const double *m, const double *v, int rowsize, int colsize, int col);

/*! \brief calculate the norm of a raw matrix.
 *
 * \return
 * zRawMatSqrNorm() returns the squared norm of a raw matrix \a m.
 * zRawMatNorm() returns the norm of \a m.
 *
 * \a rowsize and \a colsize are the row and column size of \a m, respectively.
 */
#define zRawMatSqrNorm(m,rowsize,colsize) zRawVecSqrNorm( m, (rowsize)*(colsize) )
#define zRawMatNorm(m,rowsize,colsize)    zRawVecNorm( m, (rowsize)*(colsize) )

/*! \brief transpose a raw matrix.
 *
 * zRawMatT() transposes a raw matrix \a m. The result is put into \a tm.
 * The size of \a tm must be \a rowsize x \a colsize.
 *
 * zRawMatTDRC() directly transposes \a m. The size of \a m must be \a rowsize x \a colsize.
 *
 * zRawMatTrace() calculates the trace value of \a m. The size of \a m must be \a rowsize x \a colsize.
 * \return
 * zRawMatT() and zRawMatTDRC() return no values.
 * zRawMatTrace() returns the value calculated.
 */
__ZM_EXPORT void zRawMatT(const double *m, double *tm, int rowsize, int colsize);
__ZM_EXPORT void zRawMatTDRC(double *m, int rowsize, int colsize);
__ZM_EXPORT double zRawMatTrace(const double *m, int rowsize, int colsize);

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
__ZM_EXPORT void zRawMulMatVec(const double *m, const double *v1, int rowsize, int colsize, double *v);
__ZM_EXPORT void zRawMulMatTVec(const double *m, const double *v1, int rowsize, int colsize, double *v);
__ZM_EXPORT void zRawMulMatMat(const double *m1, int rowsize1, int colsize1, const double *m2, int rowsize2, int colsize2, double *m);
__ZM_EXPORT void zRawMulMatMatT(const double *m1, int rowsize1, int colsize1, const double *m2, int rowsize2, int colsize2, double *m);
__ZM_EXPORT void zRawMulMatTMat(const double *m1, int rowsize1, int colsize1, const double *m2, int rowsize2, int colsize2, double *m);

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
 * \a m is specified as \a rowsize x \a colsize.
 *
 * zRawMatPrint() prints \a m out to the standard output.
 * \return
 * These functions return no value.
 */
__ZM_EXPORT void zRawMatFPrint(FILE *fp, const double *m, int rowsize, int colsize);
#define zRawMatPrint(m,rowsize,colsize) zRawMatFPrint( stdout, m, rowsize, colsize )

__END_DECLS

#endif /* __ZM_RAW_MAT_H__ */
