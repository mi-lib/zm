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
/* CLASS: zMat
 * double precision floating point value matrix class
 * NOTES: each elements of matrix(size=r*c) is at (0 - r-1,0 - c-1).
 * ********************************************************** */

typedef struct{
  int row, col;
  double *elem;
} zMatStruct;
typedef zMatStruct * zMat;

#define zMatRowSizeNC(m)        (m)->row
#define zMatColSizeNC(m)        (m)->col
#define zMatRowSize(m)          ( (m) ? zMatRowSizeNC(m) : 0 )
#define zMatColSize(m)          ( (m) ? zMatColSizeNC(m) : 0 )
#define zMatSetRowSize(m,r)     ( zMatRowSizeNC(m) = (r) )
#define zMatSetColSize(m,c)     ( zMatColSizeNC(m) = (c) )

#define zMatSetSize(m,r,c) do{\
  zMatSetRowSize(m,r);\
  zMatSetColSize(m,c);\
} while(0)
#define zMatRowSizeIsEqual(m1,m2) \
  ( zMatRowSize(m1) == zMatRowSize(m2) )
#define zMatColSizeIsEqual(m1,m2) \
  ( zMatColSize(m1) == zMatColSize(m2) )
#define zMatSizeIsEqual(m1,m2) \
  ( zMatRowSizeIsEqual(m1,m2) && zMatColSizeIsEqual(m1,m2) )
#define zMatColVecSizeIsEqual(m,v) \
  ( zMatColSize(m) == zVecSize(v) )
#define zMatRowVecSizeIsEqual(m,v) \
  ( zMatRowSize(m) == zVecSize(v) )
#define zMatRowColSizeIsEqual(m1,m2) \
  ( zMatRowSize(m1) == zMatColSize(m2) )
#define zMatColRowSizeIsEqual(m1,m2) zMatRowColSizeIsEqual(m2,m1)
#define zMatIsSqr(m) zMatRowColSizeIsEqual(m,m)

/* METHOD:
 * zMatBuf, zMatRowBuf
 * - convert matrix to an array of double-recision floating-point values.
 *
 * 'zMatBuf()' converts the matrix 'm' to the pointer
 * to the array of double precision floating-point values.
 *
 * 'zMatRowBuf()' converts the 'row'th row of 'm' to
 * the pointer to the array.
 * [RETURN VALUE]
 * 'zMatBuf()' and 'zMatRowBuf()' return a pointer
 * to the array converted.
 */
#define zMatBuf(m)      (m)->elem
#define zMatRowBuf(m,r) ( zMatBuf(m) + (r)*zMatColSizeNC(m) )

/* METHOD:
 * zMatElem, zMatSetElem, zMatSetElemList
 * - abstraction and set of matrix element.
 *
 * 'zMatElem()' returns the component at 'r'th row and
 * 'c'th column of a matrix 'm'.
 * 'zMatSetElem()' sets the component at 'r'th row and
 * 'c'th column of 'm' for a scalar value 'value'.
 *
 * 'zMatSetElemList()' sets all the components of 'm'
 * according to the value list given by the
 * rest of arguments.
 * [RETURN VALUE]
 * 'zMatElem()' and 'zMatSetElem()' return the
 * component at 'r'th row and 'c'th column of 'm'.
 *
 * 'zMatSetElemList()' returns a pointer 'm'.
 */
#define zMatElem(m,r,c)      zMatRowBuf(m,r)[c]
#define zMatSetElem(m,r,c,e) ( zMatElem(m,r,c) = (e) )
__EXPORT zMat zMatSetElemList(zMat m, ... );

/* METHOD:
 * zMatAlloc, zMatAllocSqr, zMatCreateList,
 * zMatFree, zMatFreeAO, zMatClear, zMatTouchup
 * - creation, destruction and cleanup of matrix.
 *
 * 'zMatAlloc()' creates a new matrix with 'row' and 'col'
 * for rows and columns, respectively, allocating enough memory.
 * And 'zMatAllocSqr()' creates a square matrix with
 * 'size' for both row and column.
 *
 * 'zMatCreateList()' creates a new matrix from the value
 * list given by arguments. 'row' and 'col' are the number of
 * rows and columns, respectively.
 *
 * 'zMatFree()' frees the given matrix 'm', freeing
 * the memory allocated.
 *
 * 'zMatFreeAO()' frees plural matrices lasted as '...'
 * at once. 'n' is the number of matrices to be freed.
 *
 * 'zMatClear()' clears the matrix 'm', set all of the
 * components for zero.
 *
 * 'zMatTouchup()' replaces all components which are less
 * than zTOL for zeros.
 * [RETURN VALUE]
 * Each of 'zMatAlloc()' and 'zMatAllocSqr()'
 * returns a pointer to newly allocated memory.
 *
 * 'zMatFree()' and 'zMatFreeAO()' return no values.
 *
 * 'zMatClear()' and 'zMatTouchup()' return a pointer to 'm'.
 * [NOTES]
 * Because of a bug in glibc, the following call does not
 * work as expected.
 *   v = zVecCreateList( 3, 1, 2, 3 );
 * It should be written in the following form.
 *   v = zVecCreateList( 3, 1.0, 2.0, 3.0 );
 */
__EXPORT zMat zMatAlloc(int row, int col);
#define zMatAllocSqr(s) zMatAlloc( (s), (s) )
__EXPORT zMat zMatCreateList(int row, int col, ...);
__EXPORT void zMatFree(zMat m);
__EXPORT void zMatFreeAO(int n, ...);
__EXPORT zMat zMatClear(zMat m);
__EXPORT zMat zMatTouchup(zMat m);

/* METHOD:
 * zMatIdentNC, zMatDiagNC, zMatIdent, zMatDiag, zMatRand
 * - identity matrix, diagonal matrix and random matrix.
 *
 * 'zMatIdentNC()' and 'zMatIdent()' makes the given
 * matrix 'm' an identity matrix.
 *  | 1.0 0.0 .  .  .     |
 *  | 0.0 1.0             |
 *  | .       .           |
 *  | .          .        |
 *  | .             . 0.0 |
 *  | .           0.0 1.0 |
 *
 * 'zMatDiagNC()' and 'zMatDiag()' makes 'm' a diagonal
 * matrix. The components are given by a vector 'd'.
 *  | d_1 0.0 .  .  .     |
 *  | 0.0 d_2             |
 *  | .       .           |
 *  | .          .        |
 *  | .             . 0.0 |
 *  | .           0.0 d_n |
 *
 * 'zMatRand()' sets all the components randomly within
 * the range from 'min' to 'max'.
 *
 * 'zMatIdentNC()', 'zMatDiagNC()' and 'zMatRandNC()' do
 * the operation without checking the size consistency.
 * [RETURN VALUE]
 * 'zMatIdentNC()', 'zMatDiagNC()', 'zMatIdent()',
 * 'zMatDiag()' and 'zMatRand()' return a pointer 'm'.
 * [NOTES]
 * Since 'zMatIdentNC()' and 'zMatDiagNC()' does not
 * check the size consistency; when 'm' is not a square
 * matrix, anything might happen.
 * If it is not urgent and you are not hasty, you should
 * use 'zMatIdent()' and 'zMatDiag()' for safety.
 */
__EXPORT zMat zMatIdentNC(zMat m);
__EXPORT zMat zMatDiagNC(zMat m, zVec d);
__EXPORT zMat zMatIdent(zMat m);
__EXPORT zMat zMatDiag(zMat m, zVec d);
__EXPORT zMat zMatRandUniform(zMat m, double min, double max);
__EXPORT zMat zMatRand(zMat m, zMat min, zMat max);

/* METHOD:
 * zMatCopyNC, zMatCopy, zMatCopyArray,
 * zMatClone, zMatCloneArray
 * - copy of matrix.
 *
 * 'zMatCopyNC()' copies the matrix 'src' to the other
 * 'dest' without checking the size consistency between 'src'
 * and 'dest'.
 *
 * 'zMatCopy()' copies the matrix 'src' to the other 'dest'.
 *
 * 'zMatCopyArray()' copies the 'row'x'col' components in
 * 'array' to the matrix 'm'.
 *
 * 'zMatClone()' creates a clone of the matrix 'm'.
 * 'zMatCloneArray()' creates a clone of the 'array'.
 * 's' is the size of 'array'.
 * [RETURN VALUE]
 * 'zMatCopyNC()' returns a pointer to 'dest'.
 *
 * 'zMatCopy()' returns a pointer to 'dest', or the null
 * pointer if the size of 'src' and 'dest' do not coincide.
 *
 * 'zMatCopyArray()' returns a pointer to 'm'.
 *
 * 'zMatClone()' and 'zMatCloneArray()' returns a pointer
 * to the newly created matrix.
 * [NOTES]
 * Since 'zMatCopyNC()' does not check the size
 * consistency, when the size of 'src' and 'dest' are
 * different from each other, anything might happen.
 * If it is not urgent and you are not hasty, you should
 * use 'zMatCopy()' for safety.
 */
__EXPORT zMat zMatCopyNC(zMat src, zMat dest);
__EXPORT zMat zMatCopy(zMat src, zMat dest);
__EXPORT zMat zMatCopyArray(double array[], int r, int c, zMat m);
__EXPORT zMat zMatClone(zMat src);
__EXPORT zMat zMatCloneArray(double array[], int r, int c);

/*! \brief partially copy a matrix.
 *
 * 'zMatGet()' gets a submatrix of 'src' from ('pr', 'pc)
 * to 'dest', while 'zMatPut()' puts 'dest' to 'src'
 * as a submatrix at ('pr', 'pc).
 *
 * It is expected that 'dest' for 'zMatGet()'(or 'src'
 * for 'zMatPut()') has the larger size than 'pr' +
 * the row size of 'src' (or 'dest') times 'pc' + the
 * column size of 'src' (or 'dest').
 * [RETURN VALUE]
 * 'zMatGetNC()' and 'zMatPutNC()' always return a
 * pointer to 'dest' without checking the size validity
 * between 'src' and 'dest', while 'zMatGet()' and
 * 'zMatPut()' return the null pointer if the size
 * of 'src' and 'dest' are inconsistent.
 * [SEE ALSO]
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

/* METHOD:
 * zMatRowNC, zMatColNC, zMatRow, zMatCol
 * zMatSetRowNC, zMatSetColNC,
 * zMatSetRow, zMatSetCol,
 * zMatSwapRowNC, zMatSwapColNC,
 * zMatSwapRow, zMatSwapCol
 * - abstraction, set and swap of row/column vector from matrix.
 *
 * 'zMatRowNC()' and 'zMatRow()' abstracts the 'row'th
 * row vector of the matrix 'm' and put it into a vector 'v'.
 * 'zMatRowNC()' does not do size checking, while
 * 'zMatRow()' does.
 *
 * 'zMatColNC()' and 'zMatCol()' abstracts the 'col'th
 * column vector of the matrix 'm' and put it into a vector 'v'.
 * 'zMatColNC()' does not do size checking, while
 * 'zMatCol()' does.
 *
 * 'zMatSetRowNC()' and 'zMatSetRow()' sets the given
 * vector 'v' for the 'row'th row vector of the matrix 'm'.
 * 'zMatSetRowNC()' does not do size checking, while
 * 'zMatSetRow()' does.
 *
 * 'zMatSetColNC()' and 'zMatSetCol()' sets the given
 * vector 'v' for the 'col'th column vector of the matrix 'm'.
 * 'zMatSetColNC()' does not do size checking, while
 * 'zMatSetCol()' does.
 *
 * 'zMatSwapRowNC()' and 'zMatSwapRow()' swaps
 * 'r1'th row and 'r2'th row of 'm'. 'zMatSwapColNC()'
 * and 'zMatSwapCol()' swaps 'c1'th column and 'c2'th
 * column of 'm'.
 * [RETURN VALUE]
 * Each of 'zMatRowNC()', 'zMatColNC()',
 * 'zMatRow()' and 'zMatCol()' returns a pointer to
 * the abstracted vector.
 *
 * Each of 'zMatSetRowNC()', 'zMatSetColNC()',
 * 'zMatSetRow()', 'zMatSetCol()', 'zMatSwapRow()'
 * and 'zMatSwapCol()' returns a pointer 'm'.
 * [NOTES]
 * The type of NC functions does calculation
 * without checking the size consistency. If it is not
 * urgent and you are not hasty, you should not use
 * them.
 */
__EXPORT zVec zMatGetRowNC(zMat m, int row, zVec v);
__EXPORT zVec zMatGetColNC(zMat m, int col, zVec v);
__EXPORT zVec zMatGetRow(zMat m, int row, zVec v);
__EXPORT zVec zMatGetCol(zMat m, int col, zVec v);
__EXPORT zMat zMatSetRowNC(zMat m, int row, zVec v);
__EXPORT zMat zMatSetColNC(zMat m, int col, zVec v);
__EXPORT zMat zMatSetRow(zMat m, int row, zVec v);
__EXPORT zMat zMatSetCol(zMat m, int col, zVec v);
__EXPORT zMat zMatSwapRowNC(zMat m, int r1, int r2);
__EXPORT zMat zMatSwapColNC(zMat m, int c1, int c2);
__EXPORT zMat zMatSwapRow(zMat m, int r1, int r2);
__EXPORT zMat zMatSwapCol(zMat m, int c1, int c2);

/* METHOD:
 * zMatShift
 * - shift diagonal values of a matrix.
 *
 * 'zMatShift()' shifts diagonal values of a matrix 'm'.
 * [RETURN VALUE]
 * 'zMatShift()' returns no value.
 */
__EXPORT void zMatShift(zMat m, double shift);

/* METHOD:
 * zMatIsEqual - comparison of two matrices.
 *
 * 'zMatIsEqual()' sees if the given two matrices 'm1' and
 * 'm2' are equal to each other.
 * [RETURN VALUE]
 * 'zMatIsEqual()' returns the true value if 'm1' equals
 * to 'm2', or the false value otherwise.
 */
__EXPORT bool zMatIsEqual(zMat m1, zMat m2);

/* METHOD:
 * zMatIsTol, zMatIsTiny - see if matrix is tiny.
 *
 * 'zMatIsTol()' returns the true value if all the
 * components of the vector 'm' are less than 'tol', or the
 * false value otherwise.
 * 'zMatIsTiny()' is the same with 'zMatIsTol()'
 * except it compares each component with zTOL(defined in
 * "zm_misc.h") instead of 'tol'.
 * [RETURN VALUE]
 * 'zMatIsTol()' and 'zMatIsTiny()' return results
 * as a boolean value.
 */
__EXPORT bool zMatIsTol(zMat m, double tol);
#define zMatIsTiny(m) zMatIsTol( (m), zTOL )

/* METHOD:
 * zMatRowReg, zMatColReg - matrix regression.
 *
 * 'zMatRowReg()' regresses the row size of matrix,
 * namely, if 'rank' is less than the row size of
 * a given matrix 'm', it regresses 'm' in row
 * direction.
 *
 * 'zMatColReg()' regresses the column size of matrix,
 * namely, if 'rank' is less than the column size of
 * a given matrix 'm', it regresses 'm' in column
 * direction.
 *
 * They directly modify 'm'.
 * [RETURN VALUE]
 * These functions return a pointer to 'm'.
 */
__EXPORT zMat zMatRowReg(zMat m, int rank);
__EXPORT zMat zMatColReg(zMat m, int rank);

/* METHOD:
 * zMatAddNC, zMatSubNC, zMatRevNC,
 * zMatMulNC, zMatDivNC, zMatCatNC,
 * zMatAddNCDRC, zMatSubNCDRC,
 * zMatRevNCDRC, zMatMulNCDRC,
 * zMatDivNCDRC, zMatCatNCDRC,
 * zMatAdd, zMatSub, zMatRev, zMatMul, zMatDiv, zMatCat,
 * zMatAddDRC, zMatSubDRC, zMatRevDRC,
 * zMatMulDRC, zMatDivDRC, zMatCatDRC
 * - basic arithmetics for matrix.
 *
 * 'zMatAddNC()' and 'zMatAdd()' add the two matrices,
 * 'm1' and 'm2', and put the result into 'm'.
 *
 * 'zMatSubNC()' and 'zMatSub()' subtract 'm2' from
 * 'm1', and put the result into 'm'.
 *
 * 'zMatRevNC()' and 'zMatRev()' reverse 'm1', and
 * put the result into 'm'.
 *
 * 'zMatMulNC()' and 'zMatMul()' multiply 'm1'
 * by a scalar value 'k', and put the result into 'm'.
 *
 * 'zMatDivNC()' and 'zMatDiv()' divide 'm1' by 'k',
 * and put the result into 'm'.
 *
 * 'zMatCatNC()' and 'zMatCat()' concatenate 'm1' by
 * 'm2' multiplied by 'k', and put the result into 'm'.
 *
 * 'zMatAddNCDRC()' and 'zMatAddDRC()' directly add
 * 'm2' to 'm1'.
 *
 * 'zMatSubNCDRC()' and 'zMatSubDRC()' directly
 * subtract 'm2' from 'm1'.
 *
 * 'zMatRevNCDRC()' and 'zMatRevDRC()' directly
 * reverse 'm'.
 *
 * 'zMatMulNCDRC()' and 'zMatMulDRC()' directly
 * multiply 'm' by 'k'.
 *
 * 'zMatDivNCDRC()' and 'zMatDivDRC()' directly
 * divide 'm' by 'k'.
 *
 * 'zMatCatNCDRC()' and 'zMatCatDRC()' directly
 * concatenate 'm1' by 'm2' multiplied by 'k'.
 * [RETURN VALUE]
 * Each of all these functions returns a pointer to the
 * result.
 * [NOTES]
 * NC-typed functions do calculation without checking
 * the size consistency. If it is not urgent and you
 * are not hasty, you should not use them.
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

/* METHOD:
 * zMatSqrNorm, zMatNorm,
 * - calculation of matrix norm.
 *
 * 'zMatSqrNorm()' calculates the squared norm of the
 * matrix 'm'. And 'zMatNorm()' calculates the norm
 * of 'm'.
 * [RETURN VALUE]
 * 'zMatSqrNorm()' and 'zMatNorm()' return the
 * value calculated.
 */
__EXPORT double zMatSqrNorm(zMat m);
#define zMatNorm(m) sqrt( zMatSqrNorm(m) )

__EXPORT double zMatInfNorm(zMat m);

/* METHOD:
 * zMatTNC, zMatT, zMatTDST
 * - transpose of matrix.
 *
 * 'zMatTNC()' and 'zMatT()' gets transpose matrix
 * of 'm' and set it into 'tm'.
 *
 * 'zMatTDST()' destructively modifies 'm' to the
 * transpose of itself.
 * [RETURN VALUE]
 * Each of 'zMatTNC()' and 'zMatT()' returns a
 * pointer to 'tm'.
 * 'zMatTDST()' returns a pointer to 'm'.
 * [NOTES]
 * 'zMatTNC()' does calculation without checking
 * the size consistency. If it is not urgent and
 * you are not hasty, you should not use them.
 */
__EXPORT zMat zMatTNC(zMat m, zMat tm);
__EXPORT zMat zMatT(zMat m, zMat tm);
__EXPORT zMat zMatTDST(zMat m);
__EXPORT zMat zMatTClone(zMat src);

/* METHOD:
 * zVecDyadNC, zVecDyad, zMatAddDyadNC, zMatAddDyad,
 * zMatCatDyadNC, zMatCatDyad
 * - dyad of vector.
 *
 * 'zVecDyadNC()' and 'zVecDyad()' calculate
 * the dyad of the two vector 'v1' and 'v2', namely,
 * 'v1' 'v2'^T. The result will be put into 'dyad'.
 *
 * 'zMatAddDyadNC()' and 'zMatAddDyad()' add dyad of
 * 'v1' and 'v2' to a matrix 'm'.
 * 'zMatSubDyadNC()' and 'zMatSubDyad()' subtract dyad
 * of 'v1' and 'v2' to a matrix 'm'.
 * 'zMatCatDyadNC()' and 'zMatCatDyad()' add the
 * multiplied dyad of 'v1' and 'v2' by 'k' to 'm'.
 *
 * The difference between functions with 'NC'
 * and without it is that the latters check if
 * the sizes of 'v1' and 'v2' coincide with
 * the row and column size of 'm', respectively.
 * [RETURN VALUE]
 * 'zVecDyadNC()' and 'zVecDyad()' return a pointer
 * 'dyad'.
 *
 * 'zMatAddDyad()', 'zMatSubDyad()' and 'zMatCatDyad()'
 * return a pointer 'm'.
 * [NOTES]
 * Since 'NC' type functions do not check the size
 * consistency between 'v1' and 'v2', you should
 * prefer 'zVecDyad()', unless it is urgent and
 * you are hasty.
 */
__EXPORT zMat zVecDyadNC(zVec v1, zVec v2, zMat dyad);
__EXPORT zMat zVecDyad(zVec v1, zVec v2, zMat dyad);
__EXPORT zMat zMatAddDyadNC(zMat m, zVec v1, zVec v2);
__EXPORT zMat zMatAddDyad(zMat m, zVec v1, zVec v2);
__EXPORT zMat zMatSubDyadNC(zMat m, zVec v1, zVec v2);
__EXPORT zMat zMatSubDyad(zMat m, zVec v1, zVec v2);
__EXPORT zMat zMatCatDyadNC(zMat m, double k, zVec v1, zVec v2);
__EXPORT zMat zMatCatDyad(zMat m, double k, zVec v1, zVec v2);

/* METHOD:
 * zMatTrNC, zMatTr
 * - trace of matrix.
 *
 * Both 'zMatTrNC()' and 'zMatTr()' return a trace value
 * of matrix 'm', Tr(m), a summation of diagonal components.
 *
 * 'm' must be a square matrix. 'zMatTrNC()' does not check
 * if 'm' is square.
 * [RETURN VALUE]
 * 'zMatTrNC()' and 'zMatTr()' return a value calculated.
 * When 'm' is not square, 'zMatTr()' returns 0.
 */
__EXPORT double zMatTrNC(zMat m);
__EXPORT double zMatTr(zMat m);

/* METHOD:
 * zMulMatVecNC, zMulVecMatNC,
 * zMulMatTVecNC,
 * zMulMatMatNC, zMulMatMatTNC,
 * zMulMatTMatNC,
 * zMulMatVec, zMulVecMat, zMulMatTVec,
 * zMulMatMat, zMulMatMatT, zMulMatTMat,
 * zMulMatVecDRC, zMulVecMatDRC,
 * zMulMatTVecDRC,
 * zMulMatMatDRC, zMulMatMatTDRC,
 * zMulMatTMatDRC
 * - multiplication of a matrix and a vector, or of two matrices.
 *
 * All these functions are for multiplication of
 * a matrix and a vector or two matrices.
 * zMul...NC family does calculation without checking
 * size consistency of vectors and matrices.
 * zMul...DRC family overrides vector or matrix
 * pointed by the given argument.
 *
 * 'zMulMatVecNC()' and 'zMulMatVec()' multiplies a
 * column vector 'v1' by matrix 'm' and put the result into 'v'.
 * 'zMulMatVecDRC()' directly multiplies 'v' by 'm'.
 *
 * 'zMulVecMatNC()' and 'zMulVecMat()' multiplies a
 * row vector 'v1' by matrix 'm' from rightside and put the
 * result into 'v'.
 * 'zMulVecMatDRC()' directly multiplies 'v' by 'm'.
 * It is obvious that the result done by these functions
 * exactly the same with 'zMulMatTVecNC()',
 * 'zMulMatTVec()', and 'zMulMatTVecDRC()'
 * by definition.
 *
 * 'zMulMatMatNC()' and 'zMulMatMat()' calculates
 * a multiplication 'm1 m2' and put it into 'm'.
 * 'zMulMatMatDRC()' directly multiplies 'm2' by 'm1'
 * from leftside. Note that 'm2' must be a square matrix.
 *
 * 'zMulMatMatTNC()' and 'zMulMatMatT()'
 * multiplies a transpose matrix of 'm2' by 'm1' and put the
 * result into 'm'.
 * 'zMulMatMatTDRC()' directly multiplies a
 * transpose matrix of 'm2' by 'm1' from leftside.
 *
 * 'zMulMatTMatNC()' and 'zMulMatTMat()'
 * multiplies 'm2' by a transpose matrix of 'm1' and put the
 * result into 'm'.
 * 'zMulMatMatTDRC()' directly multiplies 'm2' by
 * a transpose matrix of 'm1' from leftside.
 * [NOTES]
 * zMuli...DRC family requires more time for calculation
 * because temporary buffer allocation would be done as
 * inner operation of them.
 * [RETURN VALUE]
 * Each of these functions returns a pointer to the result.
 */
__EXPORT zVec zMulMatVecNC(zMat m, zVec v1, zVec v);
__EXPORT zVec zMulVecMatNC(zVec v1, zMat m, zVec v);
#define zMulMatTVecNC(m,v1,v) zMulVecMatNC(v1,m,v)
__EXPORT zMat zMulMatMatNC(zMat m1, zMat m2, zMat m);
__EXPORT zMat zMulMatMatTNC(zMat m1, zMat m2, zMat m);
__EXPORT zMat zMulMatTMatNC(zMat m1, zMat m2, zMat m);

__EXPORT zVec zMulMatVec(zMat m, zVec v1, zVec v);
__EXPORT zVec zMulVecMat(zVec v1, zMat m, zVec v);
#define zMulMatTVec(m,v1,v) zMulVecMat(v1,m,v)
__EXPORT zMat zMulMatMat(zMat m1, zMat m2, zMat m);
__EXPORT zMat zMulMatMatT(zMat m1, zMat m2, zMat m);
__EXPORT zMat zMulMatTMat(zMat m1, zMat m2, zMat m);

__EXPORT zVec zMulMatVecDRC(zMat m, zVec v);
__EXPORT zVec zMulVecMatDRC(zVec v, zMat m);
#define zMulMatTVecDRC(m,v) zMulVecMatDRC(v,m)
__EXPORT zMat zMulMatMatDRC(zMat m1, zMat m2);
__EXPORT zMat zMulMatMatTDRC(zMat m1, zMat m2);
__EXPORT zMat zMulMatTMatDRC(zMat m1, zMat m2);

/* METHOD:
 * zMatQuadNC, zMatQuad, zMatTQuadNC, zMatTQuad
 * - quadratic multiplication of matrices.
 * [SYNOPSIS]
 * zMat zMatQuadNC(zMat a, zVec w, zMat q);
 * zMat zMatQuad(zMat a, zVec w, zMat q);
 * zMat zMatTQuadNC(zMat a, zVec w, zMat q);
 * zMat zMatTQuad(zMat a, zVec w, zMat q);
 * [DESCRIPTION]
 * 'zMatQuadNC()' and 'zMatQuad()' calculate a quadratic
 * multiplication of a matrix 'a' amplified by a vector 'w'.
 * The resultant matrix 'q' forms as 'q = a diag{w} a^T'.
 *
 * 'zMatTQuadNC()' and 'zMatTQuad()' calculate a quadratic
 * multiplication of the transpose of 'a' amplified by 'w'.
 * The resultant matrix 'q' forms as 'q = a^T diag{w} a'.
 *
 * 'zMatQuad()' and 'zMatTQuad()' check if the sizes of 'a',
 * 'w' and 'q' are in consistent, while neither 'zMatQuadNC()'
 * nor 'zMatTQuadNC()' do.
 * [RETURN VALUE]
 * 'zMatQuadNC()' and 'zMatTQuadNC()'return a pointer 'q'.
 * 'zMatQuad()' and 'zMatTQuad()' also return a pointer 'q',
 * if succeeding. Otherwise, the null pointer is returned.
 */
__EXPORT zMat zMatQuadNC(zMat a, zVec w, zMat q);
__EXPORT zMat zMatQuad(zMat a, zVec w, zMat q);
__EXPORT zMat zMatTQuadNC(zMat a, zVec w, zMat q);
__EXPORT zMat zMatTQuad(zMat a, zVec w, zMat q);

/* METHOD:
 * zMatReadFile, zMatFRead, zMatRead,
 * zMatFWrite, zMatWrite, zMatImg
 * - input/output of matrix.
 *
 * 'zMatFRead()' reads a 2-dim sequence of double floating-point
 * values from the current position of the file 'fp',
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
 * where 'r' and 'c' are row and column size of matrix respectively.
 *
 * 'zMatRead()' reads a 2-dim sequence of double values according
 * to the above same format simply from the standard input.
 *
 * 'zMatReadFile()' reads a matrix from file 'filename' or
 * 'filename'.zm
 *
 * 'zMatFWrite()' writes the contents of the given matrix
 * 'm' to the current position of the file 'fp' in the above
 * format.
 *
 * 'zMatWrite()' writes the contents of the given matrix
 * 'm' simply to the standard output.
 *
 * 'zMatImg()' visualizes 'm' using one-charactor collage,
 * grading each component into nine groups represented by
 * '@Oo. ,x*M' in the ascent order.
 * This function is particularly for debug.
 * [RETURN VALUE]
 * Each of 'zMatReadFile()', 'zMatFRead()' and
 * 'zMatRead()' returns a pointer to the newly created matrix.
 *
 * 'zMatFWrite()' and 'zMatWrite()' return no values.
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
