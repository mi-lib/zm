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
 * NOTES: each element of matrix(size=r*c) is at (0 - r-1,0 - c-1).
 *//* ******************************************************* */
zArray2Class( zMatStruct, double );
typedef zMatStruct * zMat;

/*! \brief row size of a matrix.
 * \retval the row size of a matrix if \a m is not null.
 * \retval 0 if \a m is the null pointer.
 */
#define zMatRowSizeNC(m)        zArray2RowSize(m)
#define zMatRowSize(m)          ( (m) ? zMatRowSizeNC(m) : 0 )
/*! \brief column size of a matrix.
 * \retval the column size of a matrix if \a m is not null.
 * \retval 0 if \a m is the null pointer.
 */
#define zMatColSizeNC(m)        zArray2ColSize(m)
#define zMatColSize(m)          ( (m) ? zMatColSizeNC(m) : 0 )
/*! \brief set the row size of a matrix. */
#define zMatSetRowSizeNC(m,r)   ( zMatRowSizeNC(m) = (r) )
#define zMatSetRowSize(m,r)     ( (m) ? zMatSetRowSizeNC(m,r) : 0 )
/*! \brief set the column size of a matrix. */
#define zMatSetColSizeNC(m,c)   ( zMatColSizeNC(m) = (c) )
#define zMatSetColSize(m,c)     ( (m) ? zMatSetColSizeNC(m,c) : 0 )

#define zMatSetSizeNC(m,r,c) do{\
  zMatSetRowSizeNC(m,r);\
  zMatSetColSizeNC(m,c);\
} while(0)
#define zMatSetSize(m,r,c) do{\
  zMatSetRowSize(m,r);\
  zMatSetColSize(m,c);\
} while(0)

#define zMatRowSizeEqual(m1,m2)    ( zMatRowSizeNC(m1) == zMatRowSizeNC(m2) )
#define zMatColSizeEqual(m1,m2)    ( zMatColSizeNC(m1) == zMatColSizeNC(m2) )
#define zMatSizeEqual(m1,m2)       ( zMatRowSizeEqual(m1,m2) && zMatColSizeEqual(m1,m2) )
#define zMatColVecSizeEqual(m,v)   ( zMatColSizeNC(m) == zVecSizeNC(v) )
#define zMatRowVecSizeEqual(m,v)   ( zMatRowSizeNC(m) == zVecSizeNC(v) )
#define zMatRowColSizeEqual(m1,m2) ( zMatRowSizeNC(m1) == zMatColSizeNC(m2) )
#define zMatColRowSizeEqual(m1,m2) zMatRowColSizeEqual(m2,m1)
#define zMatIsSqr(m)               zMatRowColSizeEqual(m,m)

/*! \brief pointer to the array buffer of double-precision floating-point values in a matrix. */
#define zMatBufNC(m)      zArray2Buf(m)
#define zMatBuf(m)        ( (m) ? zMatBufNC(m) : NULL )
/*! \brief pointer to the \a r th row array buffer of double-precision floating-point values in a matrix. */
#define zMatRowBufNC(m,r) ( zMatBufNC(m) + (r)*zMatColSizeNC(m) )
#define zMatRowBuf(m,r)   ( (m) ? zMatRowBufNC(m,r) : NULL )

/*! \brief check if the specified row and column of a matrix is valid. */
#define zMatPosIsValid(m,r,c) zArray2PosIsValid(m,r,c)

/*! \brief get an element at the specified row and column of a matrix without checking the size. */
#define zMatElemNC(m,r,c)      zMatRowBufNC(m,r)[c]
/*! \brief get an element at the specified row and column of a matrix. */
#define zMatElem(m,r,c)        ( zMatPosIsValid(m,r,c) ? zMatElemNC(m,r,c) : 0 )
/*! \brief set an alement at the specified row and column of a matrix without checking the size. */
#define zMatSetElemNC(m,r,c,e) ( zMatElemNC(m,r,c) = (e) )
/*! \brief set an alement at the specified row and column of a matrix. */
#define zMatSetElem(m,r,c,e)   ( zMatPosIsValid(m,r,c) ? zMatSetElemNC(m,r,c,e) : 0 )

/*! \brief set elements of a matrix for values in the argument list.
 *
 * zMatSetElemList() sets all elements of a matrix \a m for the
 * values given by the argument list.
 * \return
 * zVecSetElemList() returns a pointer \a m.
 */
__ZM_EXPORT zMat zMatSetElemList(zMat m, ... );

/*! \brief create, destroy and cleanup a matrix.
 *
 * zMatAlloc() allocates memory for a new matrix with the size \a row times \a col.
 * zMatAllocSqr() allocates memory for a square matrix with the size \a size times \a size.
 *
 * zMatCreateList() creates a new matrix from the values given by the argument list.
 * \a row and \a col are the sizes of row and column, respectively.
 *
 * zMatFree() frees the matrix \a m.
 *
 * zMatFreeAtOnce() frees multiple matrices given by the argument list at once. \a n is the number
 * of matrices to be freed.
 *
 * zMatZero() sets all components of a matrix \a m for zeros.
 *
 * zMatTouchup() replaces all components less than \a tol of \a m for zeros.
 * \return
 * zMatAlloc() and zMatAllocSqr() return a pointer to the newly allocated memory.
 *
 * zMatFree() and zMatFreeAtOnce() return no values.
 *
 * zMatZero() and zMatTouchup() return a pointer \a m.
 * \notes
 * Because of a bug in glibc, the following call does not work as expected.
 *   v = zVecCreateList( 3, 1, 2, 3 );
 * It should be written as follows.
 *   v = zVecCreateList( 3, 1.0, 2.0, 3.0 );
 */
__ZM_EXPORT zMat zMatAlloc(int row, int col);
#define zMatAllocSqr(s) zMatAlloc( (s), (s) )
__ZM_EXPORT zMat zMatCreateList(int row, int col, ...);
__ZM_EXPORT void zMatFree(zMat m);
__ZM_EXPORT void zMatFreeAtOnce(int num, ...);
__ZM_EXPORT zMat zMatZero(zMat m);
__ZM_EXPORT zMat zMatTouchup(zMat m, double tol);

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
 * zMatDiagNC() and zMatDiag() make a diagonal matrix, where the elements are given by a vector \a d:
 *  | d_1 0.0 .  .  .     |
 *  | 0.0 d_2             |
 *  | .       .           |
 *  | .          .        |
 *  | .             . 0.0 |
 *  | .           0.0 d_n |
 *
 * zMatRand() makes a random matrix, all elements of which are within the range from \a min to \a max.
 *
 * zMatIdentNC(), zMatDiagNC() and zMatRandNC() do the same operation but without checking the size
 * consistency.
 * \return
 * zMatIdentNC(), zMatDiagNC(), zMatIdent(), zMatDiag() and zMatRand() return a pointer \a m.
 * \notes
 * Since zMatIdentNC() and zMatDiagNC() does not check the size consistency, anything might happen
 * if \a m is not a square matrix. If it is not urgent and you are not hasty, you'd better use
 * zMatIdent() and zMatDiag() for safety.
 */
__ZM_EXPORT zMat zMatIdentNC(zMat m);
__ZM_EXPORT zMat zMatDiagNC(zMat m, zVec d);
__ZM_EXPORT zMat zMatIdent(zMat m);
__ZM_EXPORT zMat zMatDiag(zMat m, zVec d);
__ZM_EXPORT zMat zMatRandUniform(zMat m, double min, double max);
__ZM_EXPORT zMat zMatRand(zMat m, zMat min, zMat max);

/*! \brief copy and clone a matrix.
 *
 * zMatCopy() copies the matrix \a src to the other \a dest.
 * zMatCopyNC() also copies \a src to the other \a dest but without checking the size consistency.
 *
 * zMatCopyArray() copies elements of an array with the size \a row times \a col to \a m.
 *
 * zMatClone() creates a clone of \a m.
 * zMatCloneArray() creates a clone of \a array with the size \a s.
 * \return
 * zMatCopyNC() returns a pointer \a dest.
 *
 * zMatCopy() returns a pointer \a dest, or the null pointer if the size of \a src and \a dest
 * are not consistent.
 *
 * zMatCopyArray() returns a pointer \a m.
 *
 * zMatClone() and zMatCloneArray() returns a pointer to the newly created matrix.
 * \notes
 * Since zMatCopyNC() does not check the size consistency, anything might happen if the sizes of
 * \a src and \a dest are inconsistent. If it is not urgent and you are not hasty, you'd better
 * use zMatCopy() for safety.
 */
__ZM_EXPORT zMat zMatCopyNC(const zMat src, zMat dest);
__ZM_EXPORT zMat zMatCopy(const zMat src, zMat dest);
__ZM_EXPORT zMat zMatCopyArray(const double array[], int r, int c, zMat m);
__ZM_EXPORT zMat zMatClone(const zMat src);
__ZM_EXPORT zMat zMatCloneArray(const double array[], int r, int c);

/*! \brief partially copy a matrix.
 *
 * zMatGet() gets a submatrix of \a src from (\a pr, \a pc) to \a dest, while zMatPut() puts \a src
 * to \a dest as a submatrix at (\a pr, \a pc).
 *
 * It is expected that \a dest for zMatGet() (or \a src for zMatPut()) has larger size than \a pr +
 * the row size of \a src (or \a dest) times \a pc + the column size of \a src (or \a dest).
 * \return
 * zMatGetNC() and zMatPutNC() always return a pointer \a dest without checking the size consistency
 * between \a src and \a dest, while zMatGet() and zMatPut() return the null pointer if the sizes of
 * \a src and \a dest are inconsistent.
 * \sa
 * zRawMatGet, zRawMatPut
 */
__ZM_EXPORT zMat zMatGetNC(const zMat src, int pr, int pc, zMat dest);
__ZM_EXPORT zMat zMatGet(const zMat src, int pr, int pc, zMat dest);
__ZM_EXPORT zMat zMatTGetNC(const zMat src, int pr, int pc, zMat dest);
__ZM_EXPORT zMat zMatTGet(const zMat src, int pr, int pc, zMat dest);
__ZM_EXPORT zMat zMatPutNC(zMat dest, int pr, int pc, const zMat src);
__ZM_EXPORT zMat zMatPut(zMat dest, int pr, int pc, const zMat src);
__ZM_EXPORT zMat zMatTPutNC(zMat dest, int pr, int pc, const zMat src);
__ZM_EXPORT zMat zMatTPut(zMat dest, int pr, int pc, const zMat src);

/*! \brief abstract, put and swap row/column vector of a matrix.
 *
 * zMatRowNC() abstracts the \a row th row vector of a matrix \a m and puts it into a vector \a v
 * without checking the size.
 * zMatRow() abstracts the \a row th row vector of \a m and puts it into \a v.
 *
 * zMatColNC() abstracts the \a col th column vector of \a m and puts it into \a v without checking
 * the size.
 * zMatCol() abstracts the \a col th column vector of \a m and puts it into \a v.
 *
 * zMatPutRowNC() puts the \a row th row vector of \a m for \a v without checking the size.
 * zMatPutRow() puts the \a row th row vector of \a m for \a v.
 *
 * zMatPutColNC() puts the \a col th column vector of \a m for \a v without checking the size.
 * zMatPutCol() puts the \a col th column vector of \a m for \a v.
 *
 * zMatSwapRowNC() swaps \a r1 th row and \a r2 th row of \a m without checking the size.
 * zMatSwapRow() swaps \a r1 th row and \a r2 th row of \a m.
 * zMatSwapColNC() swaps \a c1 th column and \a c2 th column of \a m without checking the size.
 * zMatSwapCol() swaps \a c1 th column and \a c2 th column of \a m.
 * \return
 * zMatRowNC(), zMatColNC(), zMatRow() and zMatCol() return a pointer to the abstracted vector.
 *
 * zMatPutRowNC(), zMatPutColNC(), zMatPutRow(), zMatPutCol(), zMatSwapRow() and zMatSwapCol()
 * return a pointer \a m.
 * \notes
 * If it is not urgent and you are not hasty, you'd better not use NC functions for safety.
 */
__ZM_EXPORT zVec zMatGetRowNC(const zMat m, int row, zVec v);
__ZM_EXPORT zVec zMatGetColNC(const zMat m, int col, zVec v);
__ZM_EXPORT zVec zMatGetRow(const zMat m, int row, zVec v);
__ZM_EXPORT zVec zMatGetCol(const zMat m, int col, zVec v);
__ZM_EXPORT zMat zMatPutRowNC(zMat m, int row, const zVec v);
__ZM_EXPORT zMat zMatPutColNC(zMat m, int col, const zVec v);
__ZM_EXPORT zMat zMatPutRow(zMat m, int row, const zVec v);
__ZM_EXPORT zMat zMatPutCol(zMat m, int col, const zVec v);
__ZM_EXPORT zMat zMatSwapRowNC(zMat m, int r1, int r2);
__ZM_EXPORT zMat zMatSwapColNC(zMat m, int c1, int c2);
__ZM_EXPORT zMat zMatSwapRow(zMat m, int r1, int r2);
__ZM_EXPORT zMat zMatSwapCol(zMat m, int c1, int c2);

/*! \brief shift diagonal values of a matrix.
 *
 * zMatShift() add the specified value \a shift to diagonal elements of a matrix \a m as offsets.
 * \return
 * zMatShift() returns no value.
 */
__ZM_EXPORT void zMatShift(zMat m, double shift);

/*! \brief maximum and minimum of matrix elements.
 *
 * zMatMaxElem() and zMatMinElem() find the maximum and minimum component of all components of a matrix
 * \a m, respectively.
 * zMatAbsMaxElem() and zMatAbsMinElem() find the component of \a m whose absolute value is the maximum
 * and minimum, respectively.
 * For those four functions, the index that gives the maximum/minimum is stored where pointed by \a im,
 * unless it is the null pointer.
 * \return
 * zMatMaxElem(), zMatMinElem(), zMatAbsMaxElem(), and zMatAbsMinElem() return the results.
 */
#define _zMatMaxElem(m,im)    zDataMax( zMatBuf(m), zMatRowSizeNC(m)*zMatColSizeNC(m), im )
#define _zMatMinElem(m,im)    zDataMin( zMatBuf(m), zMatRowSizeNC(m)*zMatColSizeNC(m), im )
#define _zMatAbsMaxElem(m,im) zDataAbsMax( zMatBuf(m), zMatRowSizeNC(m)*zMatColSizeNC(m), im )
#define _zMatAbsMinElem(m,im) zDataAbsMin( zMatBuf(m), zMatRowSizeNC(m)*zMatColSizeNC(m), im )
__ZM_EXPORT double zMatMaxElem(const zMat m, int *im);
__ZM_EXPORT double zMatMinElem(const zMat m, int *im);
__ZM_EXPORT double zMatAbsMaxElem(const zMat m, int *im);
__ZM_EXPORT double zMatAbsMinElem(const zMat m, int *im);

/*! \brief check if two matrices are equal.
 *
 * zMatEqual() checks if the given two matrices \a m1 and \a m2 are equal to each other.
 * \a tol is the tolerance to regard two values as the same.
 *
 * zMatMatch() checks fi the given two matrices \a m1 and \a m2 exactly match with each other.
 * \return
 * zMatEqual() returns the true value if \a m1 and \a m2 are equal, or the false value otherwise.
 * zMatMatch() returns the true value if \a m1 exactly matches with \a m2, or the false value otherwise.
 */
__ZM_EXPORT bool zMatEqual(const zMat m1, const zMat m2, double tol);
__ZM_EXPORT bool zMatMatch(const zMat m1, const zMat m2);

/*! \brief check if a matrix is tiny.
 *
 * zMatIsTol() checks if all elements of a matrix \a m is smaller than the tolerance \a tol.
 *
 * zMatIsTiny() checks if all elements of \a m is smaller than zTOL (defined in zm_misc.h).
 * \return
 * zMatIsTol() and zMatIsTiny() return the result as a boolean value.
 */
__ZM_EXPORT bool zMatIsTol(const zMat m, double tol);
#define zMatIsTiny(m) zMatIsTol( (m), zTOL )

/*! \brief check if a matrix is diagonal.
 *
 * zMatIsDiag() checks if a matrix \a m is diagonal, mamely, the absolute values of all non-diagonal
 * components of \a m are smaller than \a tol.
 * \return
 * zMatIsDiag() returns the true value if \a m is diagonal. Otherwise, it returns the false value.
 */
__ZM_EXPORT bool zMatIsDiag(zMat m, double tol);

/*! \brief check if a matrix is the identity matrix.
 *
 * zMatIsIdent() checks if a matrix \a m is the identity matrix, namely, all diagonal components of \a m
 * are 1 and all non-diagonal components of \a m are 0.
 * \a tol is the tolerance.
 * \return
 * zMatIsIdent() returns the true value if \a m is the identity matrix. Otherwise, it returns the false value.
 */
__ZM_EXPORT bool zMatIsIdent(zMat m, double tol);

/*! \brief check if a matrix is square and symmetric.
 *
 * zMatIsSymmetric() checks if a matrix \a m is square and symmetric.
 * \return
 * zMatIsSymmetric() returns the result as a boolean value.
 */
__ZM_EXPORT bool zMatIsSymmetric(const zMat m);

/*! \brief matrix regression.
 *
 * zMatRowReg() regresses the row size of a matrix \a m, namely, if \a rank is less than the row size
 * of \a m, it regresses \a m in row direction.
 *
 * zMatColReg() regresses the column size of \a m, namely, if \a rank is less than the column size of
 * \a m, it regresses \a m in column direction.
 *
 * Those functions directly modify \a m.
 * \return
 * zMatRowReg() and zMatColReg() return a pointer \a m.
 */
__ZM_EXPORT zMat zMatRowReg(zMat m, int rank);
__ZM_EXPORT zMat zMatColReg(zMat m, int rank);

/*! \brief basic arithmetics for matrix.
 *
 * zMatAddNC() and zMatAdd() add two matrices \a m1 and \a m2. The result is put into \a m.
 *
 * zMatSubNC() and zMatSub() subtract \a m2 from \a m1. The result is put into \a m.
 *
 * zMatRevNC() and zMatRev() reverse \a m1. The result is put into \a m.
 *
 * zMatMulNC() and zMatMul() multiply \a m1 by a scalar value \a k. The result is put into \a m.
 *
 * zMatDivNC() and zMatDiv() divide \a m1 by \a k. The result is put into \a m.
 *
 * zMatCatNC() and zMatCat() concatenate \a m1 by \a m2 multiplied by \a k. The result is put into \a m.
 *
 * zMatAddNCDRC() and zMatAddDRC() directly add \a m2 to \a m1.
 *
 * zMatSubNCDRC() and zMatSubDRC() directly subtract \a m2 from \a m1.
 *
 * zMatRevNCDRC() and zMatRevDRC() directly reverse \a m.
 *
 * zMatMulNCDRC() and zMatMulDRC() directly multiply \a m by \a k.
 *
 * zMatDivNCDRC() and zMatDivDRC() directly divide \a m by \a k.
 *
 * zMatCatNCDRC() and zMatCatDRC() directly concatenate \a m1 by \a m2 multiplied by \a k.
 * \return
 * These functions return a pointer to the result.
 * \notes
 * NC-type functions calculate without checking the size consistency. If it is not urgent and
 * you are not hasty, you'd better not use them.
 */
__ZM_EXPORT zMat zMatAddNC(const zMat m1, const zMat m2, zMat m);
__ZM_EXPORT zMat zMatSubNC(const zMat m1, const zMat m2, zMat m);
__ZM_EXPORT zMat zMatRevNC(const zMat m1, zMat m);
__ZM_EXPORT zMat zMatMulNC(const zMat m1, double k, zMat m);
__ZM_EXPORT zMat zMatDivNC(const zMat m1, double k, zMat m);
__ZM_EXPORT zMat zMatCatNC(const zMat m1, double k, const zMat m2, zMat m);

#define zMatAddNCDRC(m1,m2)   zMatAddNC( (m1), (m2), (m1) )
#define zMatSubNCDRC(m1,m2)   zMatSubNC( (m1), (m2), (m1) )
#define zMatRevNCDRC(m)       zMatRevNC( (m), (m) )
#define zMatMulNCDRC(m,k)     zMatMulNC( (m), (k) , (m) )
#define zMatDivNCDRC(m,k)     zMatDivNC( (m), (k) , (m) )
#define zMatCatNCDRC(m1,k,m2) zMatCatNC( (m1), (k) , (m2), (m1) )

__ZM_EXPORT zMat zMatAdd(const zMat m1, const zMat m2, zMat m);
__ZM_EXPORT zMat zMatSub(const zMat m1, const zMat m2, zMat m);
__ZM_EXPORT zMat zMatRev(const zMat m1, zMat m);
__ZM_EXPORT zMat zMatMul(const zMat m1, double k, zMat m);
__ZM_EXPORT zMat zMatDiv(const zMat m1, double k, zMat m);
__ZM_EXPORT zMat zMatCat(const zMat m1, double k, const zMat m2, zMat m);

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
__ZM_EXPORT double zMatSqrNorm(const zMat m);
#define zMatNorm(m) sqrt( zMatSqrNorm(m) )

__ZM_EXPORT double zMatInfNorm(const zMat m);

/*! \brief transpose a matrix.
 *
 * zMatTNC() and zMatT() transpose a matrix \a m. The result is put into \a tm.
 *
 * zMatTDRC() directly transposes \a m.
 * \return
 * zMatTNC() and zMatT() return a pointer \a tm.
 * zMatTDRC() returns a pointer \a m.
 * \notes
 * zMatTNC() does not check the size consistency. If it is not urgent and you are not hasty, you'd
 * better not use it.
 */
__ZM_EXPORT zMat zMatTNC(const zMat m, zMat tm);
__ZM_EXPORT zMat zMatT(const zMat m, zMat tm);
__ZM_EXPORT zMat zMatTDRC(zMat m);
__ZM_EXPORT zMat zMatTClone(const zMat src);

/*! \brief dyadic product of two vectors.
 *
 * zVecDyadNC() and zVecDyad() calculate the dyadic product of two vectors \a v1 and \a v2, namely,
 * \a v1 \a v2^T. The result is put into \a dyad.
 *
 * zMatAddDyadNC() and zMatAddDyad() add the dyadic product of \a v1 and \a v2 to \a m.
 * zMatSubDyadNC() and zMatSubDyad() subtract the dyad product of \a v1 and \a v2 from \a m.
 * zMatCatDyadNC() and zMatCatDyad() add the dyadic product of \a v1 and \a v2 multiplied by a scalar
 * value \a k to \a m.
 * \return
 * zVecDyadNC() and zVecDyad() return a pointer \a dyad.
 *
 * zMatAddDyad(), zMatSubDyad() and zMatCatDyad() return a pointer \a m.
 * \notes
 * NC-type functions do not check the size consistency between the vectors. If it is not urgent and
 * you are not hasty, you'd better not use them.
 */
__ZM_EXPORT zMat zVecDyadNC(const zVec v1, const zVec v2, zMat dyad);
__ZM_EXPORT zMat zVecDyad(const zVec v1, const zVec v2, zMat dyad);
__ZM_EXPORT zMat zMatAddDyadNC(zMat m, const zVec v1, const zVec v2);
__ZM_EXPORT zMat zMatAddDyad(zMat m, const zVec v1, const zVec v2);
__ZM_EXPORT zMat zMatSubDyadNC(zMat m, const zVec v1, const zVec v2);
__ZM_EXPORT zMat zMatSubDyad(zMat m, const zVec v1, const zVec v2);
__ZM_EXPORT zMat zMatCatDyadNC(zMat m, double k, const zVec v1, const zVec v2);
__ZM_EXPORT zMat zMatCatDyad(zMat m, double k, const zVec v1, const zVec v2);

/*! \brief trace of a matrix.
 *
 * \return
 * zMatTraceNC() and zMatTrace() return the trace value of a matrix \a m, i.e., the sum of diagonal components.
 * zMatTrace() returns 0 if \a m is not square.
 * \notes
 * \a m must be a square matrix.
 * zMatTraceNC() does not check if \a m is square.
 */
__ZM_EXPORT double zMatTraceNC(const zMat m);
__ZM_EXPORT double zMatTrace(const zMat m);

/*! \brief multiplication of a matrix and a vector, or of two matrices.
 *
 * zMulMatVecNC() and zMulMatVec() multiply a vector \a v1 by a matrix \a m from the left side. The
 * result is put into \a v.
 * zMulMatVecDRC() directly multiplies \a v by \a m from the left side.
 *
 * zMulMatTVecNC() and zMulMatTVec() multiply a vector \a v1 by transpose of a matrix \a m from the
 * left side. The result is put into \a v.
 * zMulVecMatDRC() directly multiplies \a v by transpose of \a m from the left side.
 *
 * zMulMatMatNC() and zMulMatMat() calculate a multiplication \a m1 \a m2. The result is put into \a m.
 *
 * zMulMatMatTNC() and zMulMatMatT() multiply transpose of \a m2 by \a m1 from the left side. The
 * result is put into \a m.
 *
 * zMulMatTMatNC() and zMulMatTMat() multiply \a m2 by transpose of \a m1. The result is put into \a m.
 * \notes
 * zMul...DRC family requires more time for calculation because temporary memory allocation is done inside.
 * \return
 * These functions return a pointer to the result.
 */
__ZM_EXPORT zVec zMulMatVecNC(const zMat m, const zVec v1, zVec v);
__ZM_EXPORT zVec zMulMatTVecNC(const zMat m, const zVec v1, zVec v);
__ZM_EXPORT zMat zMulMatMatNC(const zMat m1, const zMat m2, zMat m);
__ZM_EXPORT zMat zMulMatMatTNC(const zMat m1, const zMat m2, zMat m);
__ZM_EXPORT zMat zMulMatTMatNC(const zMat m1, const zMat m2, zMat m);

__ZM_EXPORT zVec zMulMatVec(const zMat m, const zVec v1, zVec v);
__ZM_EXPORT zVec zMulMatTVec(const zMat m, const zVec v1, zVec v);
__ZM_EXPORT zMat zMulMatMat(const zMat m1, const zMat m2, zMat m);
__ZM_EXPORT zMat zMulMatMatT(const zMat m1, const zMat m2, zMat m);
__ZM_EXPORT zMat zMulMatTMat(const zMat m1, const zMat m2, zMat m);

__ZM_EXPORT zVec zMulMatVecDRC(const zMat m, zVec v);
__ZM_EXPORT zVec zMulMatTVecDRC(const zMat m, zVec v);

/*! \brief quadratic multiplication of matrices and a weighting vector.
 *
 * zMatQuadNC() and zMatQuad() calculate a quadratic multiplication of a matrix \a a amplified by a
 * matrix diagonalizing a vector \a w.
 * The result matrix \a q forms as \a q = \a a diag{\a w} \a a^T.
 *
 * zMatTQuadNC() and zMatTQuad() calculate a quadratic multiplication of the transpose of \a a amplified
 * by a matrix diagonalizing \a w.
 * The result matrix \a q forms as \a q = \a a^T diag{\a w} \a a.
 *
 * zMatQuad() and zMatTQuad() check if the sizes of \a a, \a w and \a q are consistent, while neither
 * zMatQuadNC() nor zMatTQuadNC() do.
 * \return
 * zMatQuadNC() and zMatTQuadNC() return a pointer \a q.
 * zMatQuad() and zMatTQuad() also return a pointer \a q, if they succeed.
 * Otherwise, the null pointer is returned.
 */
__ZM_EXPORT zMat zMatQuadNC(const zMat a, const zVec w, zMat q);
__ZM_EXPORT zMat zMatQuad(const zMat a, const zVec w, zMat q);
__ZM_EXPORT zMat zMatTQuadNC(const zMat a, const zVec w, zMat q);
__ZM_EXPORT zMat zMatTQuad(const zMat a, const zVec w, zMat q);

/*! \brief quadratic multiplication of matrices.
 *
 * zMulMatMatMatTNC and zMulMatMatMatT() calculate a quadratic multiplication of a matrix \a a amplified
 * by another matrix \a q. The result matrix \a m forms as \a m = \a a \a q \a a^T.
 *
 * zMulMatTMatMatNC and zMulMatTMatMat() calculate a quadratic multiplication of a matrix \a a amplified
 * by another matrix \a q. The result matrix \a m forms as \a m = \a a^T \a q \a a.
 *
 * zMulMatMatMatT() and zMulMatTMatMat() check if the sizes of \a a, \a q and \a m are consistent, while
 * neither zMulMatMatMatTNC() nor zMulMatTMatMatNC() do.
 * \return
 * zMulMatMatMatTNC() and zMulMatTMatMatNC() return a pointer \a m.
 * zMulMatMatMatT() and zMulMatTMatMat() also return a pointer \a m, if they succeed.
 * Otherwise, the null pointer is returned.
 */
__ZM_EXPORT zMat zMulMatMatMatTNC(const zMat a, const zMat q, zMat m);
__ZM_EXPORT zMat zMulMatMatMatT(const zMat a, const zMat q, zMat m);
__ZM_EXPORT zMat zMulMatTMatMatNC(const zMat a, const zMat q, zMat m);
__ZM_EXPORT zMat zMulMatTMatMat(const zMat a, const zMat q, zMat m);

/*! \brief read a matrix from a ZTK format processor. */
__ZM_EXPORT zMat zMatFromZTK(ZTK *ztk);

/*! \brief scan and print a matrix.
 *
 * zMatFScan() scans a 2-dim sequence of double floating-point values from the current position of
 * a file \a fp, and create a new matrix. The format is as follows:
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
 * zMatScan() scans a 2-dim sequence of double values according to the above same format from the
 * standard input.
 *
 * zMatFPrint() prints a matrix \a m to the current position of a file \a fp in the above format.
 *
 * zMatPrint() prints \a m to the standard output.
 *
 * zMatImg() visualizes \a m using one-charactor collage, grading each component into nine groups
 * represented by '@Oo. ,x*M' in the ascent order. This function is particularly for debug.
 * \return
 * zMatFScan() and zMatScan() return a pointer to the newly created matrix.
 *
 * zMatFPrint() and zMatPrint() return no values.
 */
__ZM_EXPORT zMat zMatFScan(FILE *fp);
__ZM_EXPORT void zMatFPrint(FILE *fp, const zMat m);
#define zMatScan()   zMatFScan( stdin )
#define zMatPrint(m) zMatFPrint( stdout, (m) )
__ZM_EXPORT void zMatImg(const zMat m);

__END_DECLS

#endif /* __ZM_MAT_H__ */
