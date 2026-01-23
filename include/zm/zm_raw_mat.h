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
 * zRawMatZero() sets all components of a raw matrix \a m for zeroes.
 *
 * zRawMatTouchup() replaces all components less than \a tol of a raw matrix \a m for zeroes.
 *
 * \a colcapacity is the maximum column size allocated for \a m.
 * \a rowsize and \a colsize are the row and column sizes of the matrices, respectively.
 * \return
 * zRawMatZero() and zRawMatTouchup() do not return any value.
 */
__ZM_EXPORT void zRawMatZero(double *m, int colcapacity, int rowsize, int colsize);
__ZM_EXPORT void zRawMatTouchup(double *m, int colcapacity, int rowsize, int colsize, double tol);

/*! \brief copy a raw matrix.
 *
 * zRawMatCopy() copies a raw matrix \a src to another \a dest.
 * \a srccolcapacity and \a destcolcapacity are the maximum column sizes allocated for \a src and \a dest,
 * respectively.
 * \a rowsize and \a colsize are the row and column sizes of the matrices, respectively.
 * \return
 * zRawMatCopy() does not return any value.
 */
__ZM_EXPORT void zRawMatCopy(const double *src, int srccolcapacity, double *dest, int destcolcapacity, int rowsize, int colsize);

/*! \brief check if two raw matrices are equal.
 *
 * zRawMatEqual() checks if two raw matrices \a m1 and \a m2 are equal. \a tol is the tolerance to regard
 * two values are the same.
 *
 * \a colcapacity1 and \a colcapacity2 are the maximum column sizes allocated for \a m1 and \a m2,
 * respectively.
 * \a rowsize and \a colsize are the row and column sizes of the two matrices, respectively.
 * \return
 * zRawMatEqual() returns the true value if differences of all corresponding components of the two
 * matrices \a m1 and \a m2 are less than or equal to \a tol. Otherwise, it returns the false value.
 */
__ZM_EXPORT bool zRawMatEqual(const double *m1, int colcapacity1, const double *m2, int colcapacity2, int rowsize, int colsize, double tol);

/*! \brief check if two raw matrices are exactly the same with each other.
 *
 * zRawMatMatch() checks if two raw matrices \a m1 and \a m2 match with each other.
 *
 * \a colcapacity1 and \a colcapacity2 are the maximum column sizes allocated for \a m1 and \a m2,
 * respectively.
 * \a rowsize and \a colsize are the row and column sizes of the two matrices, respectively.
 * \return
 * zRawMatMatch() returns the true value if all corresponding components of the two matrices \a m1 and
 * \a m2 are exactly the same with each other. Otherwise, it returns the false value.
 */
__ZM_EXPORT bool zRawMatMatch(const double *m1, int colcapacity1, const double *m2, int colcapacity2, int rowsize, int colsize);

/*! \brief check if a matrix is tiny.
 *
 * zRawMatIsTol() checks if all components of a raw matrix \a m are smaller than a tolerance \a tol.
 * zRawMatIsTiny() checks if all components of \a m are smaller than zTOL defined in zm_misc.h.
 *
 * \a colcapacity is the maximum column size allocated for \a m.
 * \a rowsize and \a colsize are the row and column sizes of \a m, respectively.
 * \return
 * zRawMatIsTol() returns the true value if all components of \a m are smaller than \a tol. Otherwise, it
 * returns the false value.
 * zRawMatIsTiny() returns the true value if all components of \a m are smaller than zTOL. Otherwise, it
 * returns the false value.
 */
__ZM_EXPORT bool zRawMatIsTol(const double *m, int colcapacity, int rowsize, int colsize, double tol);
#define zRawMatIsTiny(m,colcapacity,rowsize,colsize) zRawMatIsTol( m, colcapacity, rowsize, colsize, zTOL )

/*! \brief create raw identity and diagonal matrices.
 *
 * zRawMatIdent() creates a raw identity matrix \a m.
 *
 * zRawMatDiag() creates a diagonal matrix \a m, referring a raw-value array \a d as the diagonal values.
 *
 * \a colcapacity is the maximum column size allocated for \a m.
 * \a rowsize and \a colsize are row and column sizes of \a m, respectively. Note that \a m does not need
 * to be a square matrix.
 *
 * \a d is an array of diagonal values. It must have more number of values than the minimum of \a rowsize
 * and \a colsize.
 * \return
 * zRawMatIdent() and zRawMatDiag() do not return any value.
 */
__ZM_EXPORT void zRawMatIdent(double *m, int colcapacity, int rowsize, int colsize);
__ZM_EXPORT void zRawMatDiag(double *m, int colcapacity, double *d, int rowsize, int colsize);

/*! \brief create a raw matrix at random.
 *
 * zRawMatRandUniform() sets all components of \a m randomly within the range from \a min to \a max.
 * zRawMatRand() sets all components of \a m randomly. Each component of \a m is within the range of
 * corresponding components of \a matmin and \a matmax.
 *
 * \a colcapacity is the maximum column size allocated for \a m.
 * \a rowsize and \a colsize are row and column sizes of \a m, respectively.
 * \return
 * zRawMatRandUniform() and zRawMatRand() do not return any value.
 */
__ZM_EXPORT void zRawMatRandUniform(double *m, int colcapacity, double min, double max, int rowsize, int colsize);
__ZM_EXPORT void zRawMatRand(double *m, int colcapacity, double *matmin, int matmincolcapacity, double *matmax, int matmaxcolcapacity, int rowsize, int colsize);

/*! \brief get/put a row/column vector in a raw matrix.
 *
 * zRawMatGetRow() gets the \a row th row vector of a raw matrix \a m into \a rowvec.
 * zRawMatGetCol() gets the \a col th column vector of a raw matrix \a m into \a colvec.
 *
 * zRawMatPutRow() puts a raw vector \a rowvec on the \a row th row of a raw matrix \a m.
 * zRawMatPutCol() puts a raw vector \a colvec on the \a col th column of \a m.
 *
 * \a colcapacity is the maximum column size allocated for \a m.
 * \a rowsize and \a colsize are row and column sizes of \a m.
 * \return
 * zRawMatGetRow(), zRawMatGetCol(), zRawMatPutRow() and zRawMatPutCol() do not return any value.
 */
__ZM_EXPORT void zRawMatGetRow(const double *m, int colcapacity, int rowsize, int colsize, int row, double *rowvec);
__ZM_EXPORT void zRawMatGetCol(const double *m, int colcapacity, int rowsize, int colsize, int col, double *colvec);
__ZM_EXPORT void zRawMatPutRow(double *m, int colcapacity, int rowsize, int colsize, int row, const double *rowvec);
__ZM_EXPORT void zRawMatPutCol(double *m, int colcapacity, int rowsize, int colsize, int col, const double *colvec);

/*! \brief swap two row/column vectors of a raw matrix.
 *
 * zRawMatSwapRow() swaps \a row1 th and \a row2 th rows of a raw matrix \a m.
 * zRawMatSwapCol() swaps \a col1 th and \a col2 th columns of \a m.
 *
 * \a colcapacity is the maximum column size allocated for \a m.
 * \a rowsize and \a colsize are row and column sizes of \a m.
 * \return
 * zRawMatSwapRow() and zRawMatSwapCol() do not return any value.
 */
__ZM_EXPORT void zRawMatSwapRow(double *m, int colcapacity, int rowsize, int colsize, int row1, int row2);
__ZM_EXPORT void zRawMatSwapCol(double *m, int colcapacity, int rowsize, int colsize, int col1, int col2);

/*! \brief partially copy a raw matrix.
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
 * These functions do not return any values.
 */
__ZM_EXPORT void zRawMatGet(const double *src, int srccolcapacity, int srcrowsize, int srccolsize, int r, int c, double *dest, int destcolcapacity, int destrowsize, int destcolsize);
__ZM_EXPORT void zRawMatPut(double *dest, int destcolcapacity, int destrowsize, int destcolsize, int r, int c, const double *src, int srccolcapacity, int srcrowsize, int srccolsize);
__ZM_EXPORT void zRawMatTGet(const double *src, int srccolcapacity, int srcrowsize, int srccolsize, int r, int c, double *dest, int destcolcapacity, int destrowsize, int destcolsize);
__ZM_EXPORT void zRawMatTPut(double *dest, int destcolcapacity, int destrowsize, int destcolsize, int r, int c, const double *src, int srccolcapacity, int srcrowsize, int srccolsize);

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
__ZM_EXPORT void zRawMatAdd(const double *m1, int colcapacity1, const double *m2, int colcapacity2, double *m, int colcapacity3, int rowsize, int colsize);
__ZM_EXPORT void zRawMatSub(const double *m1, int colcapacity1, const double *m2, int colcapacity2, double *m, int colcapacity3, int rowsize, int colsize);
__ZM_EXPORT void zRawMatRev(const double *m1, int colcapacity1, double *m, int colcapacity2, int rowsize, int colsize);
__ZM_EXPORT void zRawMatMul(const double *m1, int colcapacity1, double k, double *m, int colcapacity2, int rowsize, int colsize);
__ZM_EXPORT void zRawMatDiv(const double *m1, int colcapacity1, double k, double *m, int colcapacity2, int rowsize, int colsize);
__ZM_EXPORT void zRawMatCat(const double *m1, int colcapacity1, double k, const double *m2, int colcapacity2, double *m, int colcapacity3, int rowsize, int colsize);

__ZM_EXPORT void zRawMatAddDRC(double *m1, int colcapacity1, const double *m2, int colcapacity2, int rowsize, int colsize);
__ZM_EXPORT void zRawMatSubDRC(double *m1, int colcapacity1, const double *m2, int colcapacity2, int rowsize, int colsize);
__ZM_EXPORT void zRawMatRevDRC(double *m, int colcapacity, int rowsize, int colsize);
__ZM_EXPORT void zRawMatMulDRC(double *m, int colcapacity, double k, int rowsize, int colsize);
__ZM_EXPORT void zRawMatDivDRC(double *m, int colcapacity, double k, int rowsize, int colsize);
__ZM_EXPORT void zRawMatCatDRC(double *m1, int colcapacity1, double k, const double *m2, int colcapacity2, int rowsize, int colsize);

/*! \brief add a raw vector to a column vector of a raw matrix directly.
 *
 * zRawMatColAddDRC() adds a raw vector \a colvec to the \a col th column vector of a raw matrix \a m
 * directly.
 * \a rowsize and \a colsize are the row and column sizes of \a m, where \a rowsize is the same with
 * that of \a v.
 * \return
 * zRawMatColAddDRC() returns a pointer \a m.
 */
__ZM_EXPORT void zRawMatColAddDRC(double *m, int colcapacity, const double *colvec, int rowsize, int colsize, int col);
/*! \brief subtract a raw vector from a column vector of a raw matrix directly.
 *
 * zRawMatColSubDRC() subtracts a raw vector \a colvec from the \a col th column vector of a raw
 * matrix \a m directly.
 * \a rowsize and \a colsize are the row and column sizes of \a m, where \a rowsize is the same with
 * that of \a v.
 * \return
 * zRawMatColSubDRC() returns a pointer \a m.
 */
__ZM_EXPORT void zRawMatColSubDRC(double *m, int colcapacity, const double *colvec, int rowsize, int colsize, int col);
/*! \brief multiply a column vector of a raw matrix by a scalar value directly.
 *
 * zRawMatColMulDRC() multiplies the \a col th column vector of a raw matrix \a m by a scalar value
 * \a k directly.
 * \a rowsize and \a colsize are the row and column sizes of \a m.
 * \return
 * zRawMatColMulDRC() returns a pointer \a m.
 */
__ZM_EXPORT void zRawMatColMulDRC(double *m, int colcapacity, double k, int rowsize, int colsize, int col);
/*! \brief concatenate a raw vector multiplied by a scalar value to a column vector of a raw matrix directly.
 *
 * zRawMatColCatDRC() concatenates a raw vector \a colvec multiplied by a scalar value \a k to the \a col
 * th column vector of a raw matrix \a m directly.
 * \a rowsize and \a colsize are the row and column sizes of \a m, where \a rowsize is the same with
 * that of \a v.
 * \return
 * zRawMatColCatDRC() returns a pointer \a m.
 */
__ZM_EXPORT void zRawMatColCatDRC(double *m, int colcapacity, double k, const double *colvec, int rowsize, int colsize, int col);

/*! \brief inner product of a column vector of a raw matrix and another raw vector.
 *
 * zRawMatColInnerProd() calculates the inner product of the \a col th column vector of a raw matrix \a m
 * and another raw vector \a v.
 * \a rowsize and \a colsize are the row and column sizes of \a m, where \a rowsize is the same with
 * that of \a v.
 * \return
 * zRawMatColInnerProd() returns the calculated inner product.
 */
__ZM_EXPORT double zRawMatColInnerProd(const double *m, int colcapacity, const double *v, int rowsize, int colsize, int col);

/*! \brief calculate the norm of a raw matrix.
 *
 * \return
 * zRawMatSqrNorm() returns the squared norm of a raw matrix \a m.
 * zRawMatNorm() returns the norm of \a m.
 *
 * \a rowsize and \a colsize are the row and column size of \a m, respectively.
 */
__ZM_EXPORT double zRawMatSqrNorm(const double *m, int colcapacity, int rowsize, int colsize);
#define zRawMatNorm(m,colcapacity,rowsize,colsize) sqrt( zRawMatSqrNorm( m, colcapacity, rowsize, colsize ) )

/*! \brief transpose a raw matrix.
 *
 * zRawMatT() transposes a raw matrix \a m. The result is put into \a tm.
 * The size of \a tm must be \a rowsize x \a colsize.
 * \a colcapacity1 and \a colcapacity2 are the maximum column sizes allocated for \a m and \a tm, respectively.
 *
 * zRawMatTDRC() directly transposes \a m in-place. It swaps all components allocated for \a m.
 * The row-capacity and column-capacity are supposed to be swapped. Hence, the resulted transposed matrix
 * should be referred as (colsize x rowsize) matrix of (\a colcapacity x \a rowcapacity) array.
 * \return
 * zRawMatT() and zRawMatTDRC() do not return any value.
 */
__ZM_EXPORT void zRawMatT(const double *m, int colcapacity1, double *tm, int colcapacity2, int rowsize, int colsize);
__ZM_EXPORT void zRawMatTDRC(double *mat, int rowcapacity, int colcapacity);

/*! \brief trace of a raw matrix.
 *
 * zRawMatTrace() calculates the trace value of \a m.
 * \a colcapacity is the maximum column size allocated for \a m.
 * \a rowsize and \a colsize are row and column sizes of \a m, respectively.
 * \return
 * zRawMatTrace() returns the calculated trace value.
 */
__ZM_EXPORT double zRawMatTrace(const double *m, int colcapacity, int rowsize, int colsize);

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
__ZM_EXPORT void zRawMulMatVec(const double *m, int colcapacity, const double *v, int rowsize, int colsize, double *mv);
__ZM_EXPORT void zRawMulMatTVec(const double *m, int colcapacity, const double *v, int rowsize, int colsize, double *mv);
__ZM_EXPORT void zRawMulMatMat(const double *m1, int colcapacity1, int rowsize1, int colsize1, const double *m2, int colcapacity2, int rowsize2, int colsize2, double *m, int colcapacity3);
__ZM_EXPORT void zRawMulMatMatT(const double *m1, int colcapacity1, int rowsize1, int colsize1, const double *m2, int colcapacity2, int rowsize2, int colsize2, double *m, int colcapacity3);
__ZM_EXPORT void zRawMulMatTMat(const double *m1, int colcapacity1, int rowsize1, int colsize1, const double *m2, int colcapacity2, int rowsize2, int colsize2, double *m, int colcapacity3);

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
__ZM_EXPORT void zRawVecDyad(const double *v1, int size1, const double *v2, int size2, double *dyad, int colcapacity);
__ZM_EXPORT void zRawMatAddDyad(double *m, int colcapacity, const double *v1, int size1, const double *v2, int size2);
__ZM_EXPORT void zRawMatSubDyad(double *m, int colcapacity, const double *v1, int size1, const double *v2, int size2);
__ZM_EXPORT void zRawMatCatDyad(double *m, int colcapacity, double k, const double *v1, int size1, const double *v2, int size2);

/*! \brief print a raw matrix.
 *
 * zRawMatFPrint() prints a raw matrix \a m to the current position of a file \a fp. The size of
 * \a m is specified as \a rowsize x \a colsize.
 *
 * zRawMatPrint() prints \a m out to the standard output.
 * \return
 * These functions return no value.
 */
__ZM_EXPORT void zRawMatFPrint(FILE *fp, const double *m, int rowcapacity, int colcapacity, int rowsize, int colsize);
#define zRawMatPrint(m,rowcapacity,colcapacity,rowsize,colsize) zRawMatFPrint( stdout, m, rowcapacity, colcapacity, rowsize, colsize )

__END_DECLS

#endif /* __ZM_RAW_MAT_H__ */
