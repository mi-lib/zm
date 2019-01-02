/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_mat - matrix class.
 */

#ifndef __ZM_MAT_H__
#define __ZM_MAT_H__

#include <zm/zm_vec.h>

__BEGIN_DECLS

/* ********************************************************** */
/*! \brief double-precision floating-point value matrix class.
 *//* ******************************************************* */
zArray2Class( zMatStruct, double );
typedef zMatStruct * zMat;

/*! \brief row size of a matrix. */
#define zMatRowSizeNC(m)        zArray2RowSize(m)
/*! \brief column size of a matrix. */
#define zMatColSizeNC(m)        zArray2ColSize(m)
/*! \brief row size of a matrix.
 * \retval the row size of a matrix if \a m is not null.
 * \retval 0 if \a m is the null pointer.
 */
#define zMatRowSize(m)          ( (m) ? zMatRowSizeNC(m) : 0 )
/*! \brief column size of a matrix.
 * \retval the column size of a matrix if \a m is not null.
 * \retval 0 if \a m is the null pointer.
 */
#define zMatColSize(m)          ( (m) ? zMatColSizeNC(m) : 0 )
/*! \brief set the row size of a matrix. */
#define zMatSetRowSize(m,r)     ( zMatRowSizeNC(m) = (r) )
/*! \brief set the column size of a matrix. */
#define zMatSetColSize(m,c)     ( zMatColSizeNC(m) = (c) )

#define zMatSetSize(m,r,c) do{\
  zMatSetRowSize(m,r);\
  zMatSetColSize(m,c);\
} while(0)
#define zMatRowSizeIsEqual(m1,m2) \
  ( zMatRowSizeNC(m1) == zMatRowSizeNC(m2) )
#define zMatColSizeIsEqual(m1,m2) \
  ( zMatColSizeNC(m1) == zMatColSizeNC(m2) )
#define zMatSizeIsEqual(m1,m2) \
  ( zMatRowSizeIsEqual(m1,m2) && zMatColSizeIsEqual(m1,m2) )
#define zMatColVecSizeIsEqual(m,v) \
  ( zMatColSizeNC(m) == zVecSizeNC(v) )
#define zMatRowVecSizeIsEqual(m,v) \
  ( zMatRowSizeNC(m) == zVecSizeNC(v) )
#define zMatRowColSizeIsEqual(m1,m2) \
  ( zMatRowSizeNC(m1) == zMatColSizeNC(m2) )
#define zMatColRowSizeIsEqual(m1,m2) zMatRowColSizeIsEqual(m2,m1)
#define zMatIsSqr(m) zMatRowColSizeIsEqual(m,m)

/*! \brief pointer to the array buffer of of double-precision floating-point values in a matrix. */
#define zMatBuf(m)      zArray2Buf(m)
/*! \brief pointer to the \a r th row array buffer of of double-precision floating-point values in a matrix. */
#define zMatRowBuf(m,r) ( zMatBuf(m) + (r)*zMatColSizeNC(m) )

/*! \brief check if the specified row and column of a matrix is valid. */
#define zMatPosIsValid(m,r,c) ( (r) >= 0 && (r) < zMatRowSizeNC(m) && (c) >= 0 && (c) < zMatColSizeNC(m) )

/*! \brief get an element at the specified row and column of a matrix without checking the size. */
#define zMatElemNC(m,r,c)    zMatRowBuf(m,r)[c]
/*! \brief get an element at the specified row and column of a matrix. */
#define zMatElem(m,r,c)      ( zMatPosIsValid(m,r,c) ? zMatRowBuf(m,r)[c] : 0 )
/*! \brief set an alement at the specified row and column of a matrix without checking the size. */
#define zMatSetElemNC(m,r,c,e) ( zMatElemNC(m,r,c) = (e) )
/*! \brief set an alement at the specified row and column of a matrix. */
#define zMatSetElem(m,r,c,e) ( zMatPosIsValid(m,r,c) ? zMatSetElemNC(m,r,c,e) : 0 )

/*! \brief set elements of a matrix for values in the argument list.
 *
 * zMatSetElemList() sets all elements of a matrix \a m for the
 * values given by the argument list.
 * \return
 * zVecSetElemList() returns a pointer \a m.
 */
__EXPORT zMat zMatSetElemList(zMat m, ... );

/*! \brief create, destroy and cleanup a matrix.
 *
 * zMatAlloc() allocates memory for a new matrix with the size
 * \a row times \a col.
 * zMatAllocSqr() allocates memory for a square matrix with
 * the size \a size times \a size.
 *
 * zMatCreateList() creates a new matrix from the values given
 * by the argument list. \a row and \a col are the sizes of
 * row and column, respectively.
 *
 * zMatFree() frees the matrix \a m.
 *
 * zMatFreeAO() frees multiple matrices given by the argument
 * list at once. \a n is the number of matrices to be freed.
 *
 * zMatClear() sets all elements of the matrix \a m for zero.
 *
 * zMatTouchup() replaces all elements which are less than zTOL
 * for zeros.
 * \return
 * zMatAlloc() and zMatAllocSqr() return a pointer to the newly
 * allocated memory.
 *
 * zMatFree() and zMatFreeAO() return no values.
 *
 * zMatClear() and zMatTouchup() return a pointer \a m.
 * \notes
 * Because of a bug in glibc, the following call does not
 * work as expected.
 *   v = zVecCreateList( 3, 1, 2, 3 );
 * It should be written as follows.
 *   v = zVecCreateList( 3, 1.0, 2.0, 3.0 );
 */
__EXPORT zMat zMatAlloc(int row, int col);
#define zMatAllocSqr(s) zMatAlloc( (s), (s) )
__EXPORT zMat zMatCreateList(int row, int col, ...);
__EXPORT void zMatFree(zMat m);
__EXPORT void zMatFreeAO(int n, ...);
__EXPORT zMat zMatClear(zMat m);
__EXPORT zMat zMatTouchup(zMat m);

/*! \brief identity matrix, diagonal matrix and random matrix.
 *
 * zMatIdentNC() and zMatIdent() make an identity matrix:
 *  | 1.0 0.0 .  .  .     |
 *  | 0.0 1.0             |
 *  | .       .           |
 *  | .          .        |
 *  | .             . 0.0 |
 *  | .           0.0 1.0 |
 *
 * zMatDiagNC() and zMatDiag() make a diagonal matrix, where
 * the elements are given by a vector \a d:
 *  | d_1 0.0 .  .  .     |
 *  | 0.0 d_2             |
 *  | .       .           |
 *  | .          .        |
 *  | .             . 0.0 |
 *  | .           0.0 d_n |
 *
 * zMatRand() makes a random matrix, all elements of which are
 * within the range from \a min to \a max.
 *
 * zMatIdentNC(), zMatDiagNC() and zMatRandNC() do the same
 * operation but without checking the size consistency.
 * \return
 * zMatIdentNC(), zMatDiagNC(), zMatIdent(), zMatDiag() and
 * zMatRand() return a pointer \a m.
 * \notes
 * Since zMatIdentNC() and zMatDiagNC() does not check the size
 * consistency, anything might happen if \a m is not a square
 * matrix.
 * If it is not urgent and you are not hasty, you'd better use
 * zMatIdent() and zMatDiag() for safety.
 */
__EXPORT zMat zMatIdentNC(zMat m);
__EXPORT zMat zMatDiagNC(zMat m, zVec d);
__EXPORT zMat zMatIdent(zMat m);
__EXPORT zMat zMatDiag(zMat m, zVec d);
__EXPORT zMat zMatRandUniform(zMat m, double min, double max);
__EXPORT zMat zMatRand(zMat m, zMat min, zMat max);

/*! \brief copy and clone a matrix.
 *
 * zMatCopy() copies the matrix \a src to the other \a dest.
 * zMatCopyNC() also copies \a src to the other \a dest but
 * without checking the size consistency.
 *
 * zMatCopyArray() copies elements of an array with the size
 * \a row times \a col to \a m.
 *
 * zMatClone() creates a clone of \a m.
 * zMatCloneArray() creates a clone of \a array with the size \a s.
 * \return
 * zMatCopyNC() returns a pointer \a dest.
 *
 * zMatCopy() returns a pointer \a dest, or the null pointer
 * if the size of \a src and \a dest are not consistent.
 *
 * zMatCopyArray() returns a pointer \a m.
 *
 * zMatClone() and zMatCloneArray() returns a pointer to the
 * newly created matrix.
 * \notes
 * Since zMatCopyNC() does not check the size consistency,
 * anything might happen if the sizes of \a src and \a dest are
 * inconsistent. If it is not urgent and you are not hasty, you'd
 * better use zMatCopy() for safety.
 */
__EXPORT zMat zMatCopyNC(zMat src, zMat dest);
__EXPORT zMat zMatCopy(zMat src, zMat dest);
__EXPORT zMat zMatCopyArray(double array[], int r, int c, zMat m);
__EXPORT zMat zMatClone(zMat src);
__EXPORT zMat zMatCloneArray(double array[], int r, int c);

/*! \brief partially copy a matrix.
 *
 * zMatGet() gets a submatrix of \a src from (\a pr, \a pc)
 * to \a dest, while zMatPut() puts \a src to \a dest as a
 * submatrix at (\a pr, \a pc).
 *
 * It is expected that \a dest for zMatGet() (or \a src for
 * zMatPut()) has larger size than \a pr + the row size of
 * \a src (or \a dest) times \a pc + the column size of \a src
 * (or \a dest).
 * \return
 * zMatGetNC() and zMatPutNC() always return a pointer \a dest
 * without checking the size consistency between \a src and
 * \a dest, while zMatGet() and zMatPut() return the null
 * pointer if the sizes of \a src and \a dest are inconsistent.
 * \sa
 * zRawMatGet, zRawMatPut
 */
__EXPORT zMat zMatGetNC(zMat src, int pr, int pc, zMat dest);
__EXPORT zMat zMatGet(zMat src, int pr, int pc, zMat dest);
__EXPORT zMat zMatTGetNC(zMat src, int pr, int pc, zMat dest);
__EXPORT zMat zMatTGet(zMat src, int pr, int pc, zMat dest);
__EXPORT zMat zMatPutNC(zMat dest, int pr, int pc, zMat src);
__EXPORT zMat zMatPut(zMat dest, int pr, int pc, zMat src);
__EXPORT zMat zMatTPutNC(zMat dest, int pr, int pc, zMat src);
__EXPORT zMat zMatTPut(zMat dest, int pr, int pc, zMat src);

/*! \brief abstract, put and swap row/column vector of a matrix.
 *
 * zMatRowNC() abstracts the \a row'th row vector of a matrix
 * \a m and puts it into a vector \a v without checking the size.
 * zMatRow() abstracts the \a row'th row vector of \a m and puts
 * it into \a v.
 *
 * zMatColNC() abstracts the \a col'th column vector of \a m
 * and puts it into \a v without checking the size.
 * zMatCol() abstracts the \a col'th column vector of \a m and
 * puts it into \a v.
 *
 * zMatPutRowNC() puts the \a row'th row vector of \a m for
 * \a v without checking the size.
 * zMatPutRow() puts the \a row'th row vector of \a m for \a v.
 *
 * zMatPutColNC() puts the \a col'th column vector of \a m for
 * \a v without checking the size.
 * zMatPutCol() puts the \a col'th column vector of \a m for \a v.
 *
 * zMatSwapRowNC() swaps \a r1'th row and \a r2'th row of \a m
 * without checking the size.
 * zMatSwapRow() swaps \a r1'th row and \a r2'th row of \a m.
 * zMatSwapColNC() swaps \a c1'th column and \a c2'th column of
 * \a m without checking the size.
 * zMatSwapCol() swaps \a c1'th column and \a c2'th column of
 * \a m.
 * \return
 * zMatRowNC(), zMatColNC(), zMatRow() and zMatCol() return a
 * pointer to the abstracted vector.
 *
 * zMatPutRowNC(), zMatPutColNC(), zMatPutRow(), zMatPutCol(),
 * zMatSwapRow() and zMatSwapCol() return a pointer \a m.
 * \notes
 * If it is not urgent and you are not hasty, you'd better not
 * use NC functions for safety.
 */
__EXPORT zVec zMatGetRowNC(zMat m, int row, zVec v);
__EXPORT zVec zMatGetColNC(zMat m, int col, zVec v);
__EXPORT zVec zMatGetRow(zMat m, int row, zVec v);
__EXPORT zVec zMatGetCol(zMat m, int col, zVec v);
__EXPORT zMat zMatPutRowNC(zMat m, int row, zVec v);
__EXPORT zMat zMatPutColNC(zMat m, int col, zVec v);
__EXPORT zMat zMatPutRow(zMat m, int row, zVec v);
__EXPORT zMat zMatPutCol(zMat m, int col, zVec v);
__EXPORT zMat zMatSwapRowNC(zMat m, int r1, int r2);
__EXPORT zMat zMatSwapColNC(zMat m, int c1, int c2);
__EXPORT zMat zMatSwapRow(zMat m, int r1, int r2);
__EXPORT zMat zMatSwapCol(zMat m, int c1, int c2);

/*! \brief shift diagonal values of a matrix.
 *
 * zMatShift() add the specified value \a shift to diagonal
 * elements of a matrix \a m as offsets.
 * \return
 * zMatShift() returns no value.
 */
__EXPORT void zMatShift(zMat m, double shift);

/*! \brief check if two matrices are equal.
 *
 * zMatIsEqual() checks if given two matrices \a m1 and \a m2
 * are equal.
 * \return
 * zMatIsEqual() returns the true value if \a m1 and \a m2
 * are equal, or the false value otherwise.
 */
__EXPORT bool zMatIsEqual(zMat m1, zMat m2);

/*! \brief check if a matrix is tiny.
 *
 * zMatIsTol() checks if all elements of a matrix \a m is
 * smaller than the tolerance \a tol.
 *
 * zMatIsTiny() checks if all elements of \a m is smaller
 * than zTOL, which is defined in zm_misc.h.
 * \return
 * zMatIsTol() and zMatIsTiny() return the result as a
 * boolean value.
 */
__EXPORT bool zMatIsTol(zMat m, double tol);
#define zMatIsTiny(m) zMatIsTol( (m), zTOL )

/*! \brief matrix regression.
 *
 * zMatRowReg() regresses the row size of a matrix \a m,
 * namely, if \a rank is less than the row size of \a m,
 * it regresses \a m in row direction.
 *
 * zMatColReg() regresses the column size of \a m, namely,
 * if \a rank is less than the column size of \a m, it
 * regresses \a m in column direction.
 *
 * Those functions directly modify \a m.
 * \return
 * zMatRowReg() and zMatColReg() return a pointer \a m.
 */
__EXPORT zMat zMatRowReg(zMat m, int rank);
__EXPORT zMat zMatColReg(zMat m, int rank);

/*! \brief basic arithmetics for matrix.
 *
 * zMatAddNC() and zMatAdd() add two matrices \a m1 and \a m2.
 * The result is put into \a m.
 *
 * zMatSubNC() and zMatSub() subtract \a m2 from \a m1.
 * The result is put into \a m.
 *
 * zMatRevNC() and zMatRev() reverse \a m1. The result is put
 * into \a m.
 *
 * zMatMulNC() and zMatMul() multiply \a m1 by a scalar value
 * \a k. The result is put into \a m.
 *
 * zMatDivNC() and zMatDiv() divide \a m1 by \a k. The result
 * is put into \a m.
 *
 * zMatCatNC() and zMatCat() concatenate \a m1 by \a m2
 * multiplied by \a k. The result is put into \a m.
 *
 * zMatAddNCDRC() and zMatAddDRC() directly add \a m2 to \a m1.
 *
 * zMatSubNCDRC() and zMatSubDRC() directly subtract \a m2
 * from \a m1.
 *
 * zMatRevNCDRC() and zMatRevDRC() directly reverse \a m.
 *
 * zMatMulNCDRC() and zMatMulDRC() directly multiply \a m by
 * \a k.
 *
 * zMatDivNCDRC() and zMatDivDRC() directly divide \a m by \a k.
 *
 * zMatCatNCDRC() and zMatCatDRC() directly concatenate \a m1
 * by \a m2 multiplied by \a k.
 * \return
 * These functions return a pointer to the result.
 * \notes
 * NC-type functions calculate without checking the size
 * consistency. If it is not urgent and you are not hasty,
 * you'd better not use them.
 */
__EXPORT zMat zMatAddNC(zMat m1, zMat m2, zMat m);
__EXPORT zMat zMatSubNC(zMat m1, zMat m2, zMat m);
__EXPORT zMat zMatRevNC(zMat m1, zMat m);
__EXPORT zMat zMatMulNC(zMat m1, double k, zMat m);
__EXPORT zMat zMatDivNC(zMat m1, double k, zMat m);
__EXPORT zMat zMatCatNC(zMat m1, double k, zMat m2, zMat m);

#define zMatAddNCDRC(m1,m2)   zMatAddNC( (m1), (m2), (m1) )
#define zMatSubNCDRC(m1,m2)   zMatSubNC( (m1), (m2), (m1) )
#define zMatRevNCDRC(m)       zMatRevNC( (m), (m) )
#define zMatMulNCDRC(m,k)     zMatMulNC( (m), (k) , (m) )
#define zMatDivNCDRC(m,k)     zMatDivNC( (m), (k) , (m) )
#define zMatCatNCDRC(m1,k,m2) zMatCatNC( (m1), (k) , (m2), (m1) )

__EXPORT zMat zMatAdd(zMat m1, zMat m2, zMat m);
__EXPORT zMat zMatSub(zMat m1, zMat m2, zMat m);
__EXPORT zMat zMatRev(zMat m1, zMat m);
__EXPORT zMat zMatMul(zMat m1, double k, zMat m);
__EXPORT zMat zMatDiv(zMat m1, double k, zMat m);
__EXPORT zMat zMatCat(zMat m1, double k, zMat m2, zMat m);

#define zMatAddDRC(m1,m2)     zMatAdd( (m1), (m2), (m1) )
#define zMatSubDRC(m1,m2)     zMatSub( (m1), (m2), (m1) )
#define zMatRevDRC(m)         zMatRev( (m), (m) )
#define zMatMulDRC(m,k)       zMatMul( (m), (k), (m) )
#define zMatDivDRC(m,k)       zMatDiv( (m), (k), (m) )
#define zMatCatDRC(m1,k,m2)   zMatCat( (m1), (k), (m2), (m1) )

/*! \brief calculate the norm of a matrix.
 *
 * \return
 * zMatSqrNorm() returns the squared norm of a matrix \a m.
 * \return
 * zMatNorm() returns the norm of \a m.
 */
__EXPORT double zMatSqrNorm(zMat m);
#define zMatNorm(m) sqrt( zMatSqrNorm(m) )

__EXPORT double zMatInfNorm(zMat m);

/*! \brief transpose a matrix.
 *
 * zMatTNC() and zMatT() transpose a matrix \a m.
 * The result is put into \a tm.
 *
 * zMatTDRC() directly transposes \a m.
 * \return
 * zMatTNC() and zMatT() return a pointer \a tm.
 * zMatTDRC() returns a pointer \a m.
 * \notes
 * zMatTNC() does not check the size consistency.
 * If it is not urgent and you are not hasty, you'd better
 * not use it.
 */
__EXPORT zMat zMatTNC(zMat m, zMat tm);
__EXPORT zMat zMatT(zMat m, zMat tm);
__EXPORT zMat zMatTDRC(zMat m);
__EXPORT zMat zMatTClone(zMat src);

/*! \brief dyadic product of two vectors.
 *
 * zVecDyadNC() and zVecDyad() calculate the dyadic product
 * of two vectors \a v1 and \a v2, namely, \a v1 \a v2^T.
 * The result is put into \a dyad.
 *
 * zMatAddDyadNC() and zMatAddDyad() add the dyadic product
 * of \a v1 and \a v2 to \a m.
 * zMatSubDyadNC() and zMatSubDyad() subtract the dyad
 * product of \a v1 and \a v2 from \a m.
 * zMatCatDyadNC() and zMatCatDyad() add the dyadic product
 * of \a v1 and \a v2 multiplied by a scalar value \a k to
 * \a m.
 * \return
 * zVecDyadNC() and zVecDyad() return a pointer \a dyad.
 *
 * zMatAddDyad(), zMatSubDyad() and zMatCatDyad() return
 * a pointer \a m.
 * \notes
 * NC-type functions do not check the size consistency
 * between the vectors. If it is not urgent and you are
 * not hasty, you'd better not use them.
 */
__EXPORT zMat zVecDyadNC(zVec v1, zVec v2, zMat dyad);
__EXPORT zMat zVecDyad(zVec v1, zVec v2, zMat dyad);
__EXPORT zMat zMatAddDyadNC(zMat m, zVec v1, zVec v2);
__EXPORT zMat zMatAddDyad(zMat m, zVec v1, zVec v2);
__EXPORT zMat zMatSubDyadNC(zMat m, zVec v1, zVec v2);
__EXPORT zMat zMatSubDyad(zMat m, zVec v1, zVec v2);
__EXPORT zMat zMatCatDyadNC(zMat m, double k, zVec v1, zVec v2);
__EXPORT zMat zMatCatDyad(zMat m, double k, zVec v1, zVec v2);

/*! \brief trace of a matrix.
 *
 * \return
 * zMatTrNC() and zMatTr() return the trace value of a matrix
 * \a m, i.e., the sum of diagonal components.
 * zMatTr() returns 0 if \a m is not square.
 * \notes
 * \a m must be a square matrix.
 * zMatTrNC() does not check if \a m is square.
 */
__EXPORT double zMatTrNC(zMat m);
__EXPORT double zMatTr(zMat m);

/*! \brief multiplication of a matrix and a vector, or of two matrices.
 *
 * zMulMatVecNC() and zMulMatVec() multiply a vector \a v1
 * by a matrix \a m from the left side. The result is put
 * into \a v.
 * zMulMatVecDRC() directly multiplies \a v by \a m from
 * the left side.
 *
 * zMulMatTVecNC() and zMulMatTVec() multiply a vector \a v1
 * by transpose of a matrix \a m from the left side. The
 * result is put into \a v.
 * zMulVecMatDRC() directly multiplies \a v by transpose
 * of \a m from the left side.
 *
 * zMulMatMatNC() and zMulMatMat() calculate a multiplication
 * \a m1 \a m2. The result is put into \a m.
 *
 * zMulMatMatTNC() and zMulMatMatT() multiply transpose of
 * \a m2 by \a m1 from the left side. The result is put into
 * \a m.
 *
 * zMulMatTMatNC() and zMulMatTMat() multiply \a m2 by
 * transpose of \a m1. The result is put into \a m.
 * \notes
 * zMul...DRC family requires more time for calculation
 * because temporary memory allocation is done inside.
 * \return
 * These functions return a pointer to the result.
 */
__EXPORT zVec zMulMatVecNC(zMat m, zVec v1, zVec v);
__EXPORT zVec zMulMatTVecNC(zMat m, zVec v1, zVec v);
__EXPORT zMat zMulMatMatNC(zMat m1, zMat m2, zMat m);
__EXPORT zMat zMulMatMatTNC(zMat m1, zMat m2, zMat m);
__EXPORT zMat zMulMatTMatNC(zMat m1, zMat m2, zMat m);

__EXPORT zVec zMulMatVec(zMat m, zVec v1, zVec v);
__EXPORT zVec zMulMatTVec(zMat m, zVec v1, zVec v);
__EXPORT zMat zMulMatMat(zMat m1, zMat m2, zMat m);
__EXPORT zMat zMulMatMatT(zMat m1, zMat m2, zMat m);
__EXPORT zMat zMulMatTMat(zMat m1, zMat m2, zMat m);

__EXPORT zVec zMulMatVecDRC(zMat m, zVec v);
__EXPORT zVec zMulMatTVecDRC(zMat m, zVec v);

/*! \brief quadratic multiplication of matrices.
 *
 * zMatQuadNC() and zMatQuad() calculate a quadratic
 * multiplication of a matrix \a a amplified by a vector
 * \a w.
 * The resultant matrix \a q forms as \a q = \a a diag{\a w} \a a^T.
 *
 * zMatTQuadNC() and zMatTQuad() calculate a quadratic
 * multiplication of the transpose of \a a amplified by
 * \a w.
 * The resultant matrix \a q forms as \a q = \a a^T diag{\a w} \a a.
 *
 * zMatQuad() and zMatTQuad() check if the sizes of \a a, \a w
 * and \a q are consistent, while neither zMatQuadNC() nor
 * zMatTQuadNC() do.
 * \return
 * zMatQuadNC() and zMatTQuadNC()return a pointer \a q.
 * zMatQuad() and zMatTQuad() also return a pointer \a q,
 * if succeeding. Otherwise, the null pointer is returned.
 */
__EXPORT zMat zMatQuadNC(zMat a, zVec w, zMat q);
__EXPORT zMat zMatQuad(zMat a, zVec w, zMat q);
__EXPORT zMat zMatTQuadNC(zMat a, zVec w, zMat q);
__EXPORT zMat zMatTQuad(zMat a, zVec w, zMat q);

/*! \brief input/output of matrix.
 *
 * zMatFRead() reads a 2-dim sequence of double floating-point
 * values from the current position of the file \a fp,
 * and create a new matrix.
 * The format is as follows:
 *  (r, c) {
 *   x11 x12 ... x1c
 *   x21 x22 ... x2c
 *    .   .  .    .
 *    .   .   .   .
 *    .   .    .  .
 *   xr1 xr2 ... xrc
 *  }
 * where \a r and \a c are row and column size of matrix respectively.
 *
 * zMatRead() reads a 2-dim sequence of double values according
 * to the above same format from the standard input.
 *
 * zMatReadFile() reads a matrix from file \a filename or
 * \a filename.zm
 *
 * zMatFWrite() writes the contents of the given matrix
 * \a m to the current position of the file 'fp' in the above
 * format.
 *
 * zMatWrite() writes the contents of the given matrix \a m
 * to the standard output.
 *
 * zMatImg() visualizes \a m using one-charactor collage,
 * grading each component into nine groups represented by
 * '@Oo. ,x*M' in the ascent order.
 * This function is particularly for debug.
 * \return
 * zMatReadFile(), zMatFRead() and zMatRead() return a pointer
 * to the newly created matrix.
 *
 * zMatFWrite() and zMatWrite() return no values.
 */
#define ZMATRIX_SUFFIX "zm"
__EXPORT zMat zMatReadFile(char filename[]);
__EXPORT zMat zMatFRead(FILE *fp);
__EXPORT void zMatFWrite(FILE *fp, zMat m);
#define zMatRead()   zMatFRead( stdin )
#define zMatWrite(m) zMatFWrite( stdout, (m) )
__EXPORT void zMatImg(zMat m);

__END_DECLS

#endif /* __ZM_MAT_H__ */
