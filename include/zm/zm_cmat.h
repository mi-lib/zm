/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_cmat - complex matrix class.
 */

#ifndef __ZM_CMAT_H__
#define __ZM_CMAT_H__

#include <zm/zm_mat.h>
#include <zm/zm_cvec.h>

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zCMat
 * complex matrix class
 * NOTES: each elements of matrix(size=r*c) is at (0 - r-1,0 - c-1).
 * ********************************************************** */

typedef struct{
  int row, col;
  zComplex *elem;
} zCMatStruct;
typedef zCMatStruct * zCMat;

#define _zCMatRowSize(m)         (m)->row
#define _zCMatColSize(m)         (m)->col
#define zCMatRowSize(m)          ( (m) ? _zCMatRowSize(m) : 0 )
#define zCMatColSize(m)          ( (m) ? _zCMatColSize(m) : 0 )
#define zCMatSetRowSize(m,r)     ( _zCMatRowSize(m) = (r) )
#define zCMatSetColSize(m,c)     ( _zCMatColSize(m) = (c) )

#define zCMatSetSize(m,r,c) do{\
  zCMatSetRowSize(m,r);\
  zCMatSetColSize(m,c);\
} while(0)

#define zCMatSizeIsEqual(m1,m2) \
  ( _zCMatRowSize(m1) == _zCMatRowSize(m2) && \
    _zCMatColSize(m1) == _zCMatColSize(m2) )

/* METHOD:
 * zCMatElem, zCMatSetElem
 * - abstraction and set of matrix element.
 *
 * 'zCMatElem()' returns the component at 'r'th row and
 * 'c'th column of a matrix 'm'.
 * 'zCMatSetElem()' sets the component at 'r'th row and
 * 'c'th column of 'm' for a scalar value 'value'.
 * [RETURN VALUE]
 * 'zCMatElem()' and 'zCMatSetElem()' return the
 * component at 'r'th row and 'c'th column of 'm'.
 */
#define zCMatElem(m,r,c)      ( &(m)->elem[(r)*_zCMatColSize(m)+(c)] )
#define zCMatSetElem(m,r,c,e) zComplexCopy( (e), zCMatElem(m,r,c) )

/* METHOD:
 * zCMatAlloc, zCMatAllocSqr, zCMatFree, zCMatClear
 * - creation, destruction and cleanup of matrix.
 *
 * 'zCMatAlloc()' creates a new matrix with 'row' and 'col'
 * for rows and columns, respectively, allocating enough memory.
 * And 'zCMatAllocSqr()' creates a square matrix with
 * 'size' for both row and column.
 * #
 * 'zCMatFree()' frees the memory allocated for \a m.
 * #
 * 'zCMatClear()' clears the matrix 'm', set all of the
 * components for zero.
 * [RETURN VALUE]
 * Each of 'zCMatAlloc()' and 'zCMatAllocSqr()'
 * returns a pointer to newly allocated memory.
 * #
 * 'zCMatFree()' returns no values.
 * #
 * 'zCMatClear()' returns a pointer to 'm'.
 */
__EXPORT zCMat zCMatAlloc(int row, int col);
#define zCMatAllocSqr(s) zCMatAlloc( (s), (s) )
__EXPORT void zCMatFree(zCMat m);
__EXPORT zCMat zCMatClear(zCMat m);

/* METHOD:
 * zCMatCopyNC, zCMatCopy, zCMatClone
 * - copy of matrix.
 *
 * 'zCMatCopyNC()' copies the matrix 'src' to the other
 * 'dest' without checking the size consistency between 'src'
 * and 'dest'.
 * #
 * 'zCMatCopy()' copies the matrix 'src' to the other 'dest'.
 * #
 * 'zCMatClone()' creates a clone of the matrix 'm'.
 * [RETURN VALUE]
 * 'zCMatCopyNC()' returns a pointer to 'dest'.
 * #
 * 'zCMatCopy()' returns a pointer to 'dest', or the null
 * pointer if the size of 'src' and 'dest' do not coincide.
 * #
 * 'zCMatClone()' returns a pointer
 * to the newly created matrix.
 * [NOTES]
 * Since 'zCMatCopyNC()' does not check the size
 * consistency, when the size of 'src' and 'dest' are
 * different from each other, anything might happen.
 * If it is not urgent and you are not hasty, you should
 * use 'zCMatCopy()' for safety.
 */
__EXPORT zCMat zCMatCopyNC(zCMat src, zCMat dest);
__EXPORT zCMat zCMatCopy(zCMat src, zCMat dest);
__EXPORT zCMat zCMatClone(zCMat src);

__EXPORT zCMat zMat2CMat(zMat m, zCMat cm);

/* METHOD:
 * zCMatIsTol, zCMatIsTiny - see if matrix is tiny.
 *
 * 'zCMatIsTol()' returns the true value if all the
 * components of the vector 'm' are less than 'tol', or the
 * false value otherwise.
 * 'zCMatIsTiny()' is the same with 'zCMatIsTol()'
 * except it compares each component with zTOL(defined in
 * "zm_misc.h") instead of 'tol'.
 * [RETURN VALUE]
 * 'zCMatIsTol()' and 'zCMatIsTiny()' return results
 * as a boolean value.
 */
__EXPORT bool zCMatIsTol(zCMat m, double tol);
#define zCMatIsTiny(m) zCMatIsTol( (m), zTOL )

/* METHOD:
 * zCMatAddNC, zCMatSubNC, zCMatRevNC,
 * zCMatMulNC, zCMatDivNC,
 * zCMatAddNCDRC, zCMatSubNCDRC,
 * zCMatRevNCDRC, zCMatMulNCDRC,
 * zCMatDivNCDRC,
 * zCMatAdd, zCMatSub, zCMatRev, zCMatMul, zCMatDiv,
 * zCMatAddDRC, zCMatSubDRC, zCMatRevDRC,
 * zCMatMulDRC, zCMatDivDRC,
 * - basic arithmetics for matrix.
 *
 * 'zCMatAddNC()' and 'zCMatAdd()' add the two matrices,
 * 'm1' and 'm2', and put the result into 'm'.
 * #
 * 'zCMatSubNC()' and 'zCMatSub()' subtract 'm2' from
 * 'm1', and put the result into 'm'.
 * #
 * 'zCMatRevNC()' and 'zCMatRev()' reverse 'm1', and
 * put the result into 'm'.
 * #
 * 'zCMatMulNC()' and 'zCMatMul()' multiply 'm1'
 * by a scalar value 'k', and put the result into 'm'.
 * #
 * 'zCMatDivNC()' and 'zCMatDiv()' divide 'm1' by 'k',
 * and put the result into 'm'.
 * #
 * 'zCMatAddNCDRC()' and 'zCMatAddDRC()' directly add
 * 'm2' to 'm1'.
 * #
 * 'zCMatSubNCDRC()' and 'zCMatSubDRC()' directly
 * subtract 'm2' from 'm1'.
 * #
 * 'zCMatRevNCDRC()' and 'zCMatRevDRC()' directly
 * reverse 'm'.
 * #
 * 'zCMatMulNCDRC()' and 'zCMatMulDRC()' directly
 * multiply 'm' by 'k'.
 * #
 * 'zCMatDivNCDRC()' and 'zCMatDivDRC()' directly
 * divide 'm' by 'k'.
 * [RETURN VALUE]
 * Each of all these functions returns a pointer to the
 * result.
 * [NOTES]
 * NC-typed functions do calculation without checking
 * the size consistency. If it is not urgent and you
 * are not hasty, you should not use them.
 */
__EXPORT zCMat zCMatAddNC(zCMat m1, zCMat m2, zCMat m);
__EXPORT zCMat zCMatSubNC(zCMat m1, zCMat m2, zCMat m);
__EXPORT zCMat zCMatRevNC(zCMat m1, zCMat m);
__EXPORT zCMat zCMatMulNC(zCMat m1, zComplex *z, zCMat m);
__EXPORT zCMat zCMatDivNC(zCMat m1, zComplex *z, zCMat m);

#define zCMatAddNCDRC(m1,m2)   zCMatAddNC( (m1), (m2), (m1) )
#define zCMatSubNCDRC(m1,m2)   zCMatSubNC( (m1), (m2), (m1) )
#define zCMatRevNCDRC(m)       zCMatRevNC( (m), (m) )
#define zCMatMulNCDRC(m,z)     zCMatMulNC( (m), (z) , (m) )
#define zCMatDivNCDRC(m,z)     zCMatDivNC( (m), (z) , (m) )

__EXPORT zCMat zCMatAdd(zCMat m1, zCMat m2, zCMat m);
__EXPORT zCMat zCMatSub(zCMat m1, zCMat m2, zCMat m);
__EXPORT zCMat zCMatRev(zCMat m1, zCMat m);
__EXPORT zCMat zCMatMul(zCMat m1, zComplex *z, zCMat m);
__EXPORT zCMat zCMatDiv(zCMat m1, zComplex *z, zCMat m);

#define zCMatAddDRC(m1,m2)     zCMatAdd( (m1), (m2), (m1) )
#define zCMatSubDRC(m1,m2)     zCMatSub( (m1), (m2), (m1) )
#define zCMatRevDRC(m)         zCMatRev( (m), (m) )
#define zCMatMulDRC(m,z)       zCMatMul( (m), (z), (m) )
#define zCMatDivDRC(m,z)       zCMatDiv( (m), (z), (m) )

/* METHOD:
 * zCMulMatVecNC, zCMulMatVec
 * - multiplication of a matrix and a vector, or of two matrices.
 *
 * 'zCMulMatVecNC()' and 'zCMulMatVec()' multiplies a
 * column vector 'v1' by matrix 'm' and put the result into 'v'.
 * 'zCMulMatVecNC()' does calculation without checking
 * size consistency of vectors and matrices.
 * [RETURN VALUE]
 * These functions return a pointer to the result.
 */
__EXPORT zCVec zCMulMatVecNC(zCMat m, zCVec v1, zCVec v);
__EXPORT zCVec zCMulMatVec(zCMat m, zCVec v1, zCVec v);

/* METHOD:
 * zCMatFWrite, zCMatWrite
 * - output of matrix.
 *
 * 'zCMatFWrite()' writes the contents of the given matrix
 * 'm' to the current position of the file 'fp' in the following
 * format.
 *  (r, c) {
 *   z11 z12 ... z1c
 *   z21 z22 ... z2c
 *    .   .  .    .
 *    .   .   .   .
 *    .   .    .  .
 *   zr1 zr2 ... zrc
 *  }
 * where 'r' and 'c' are row and column size of matrix respectively.
 * #
 * 'zCMatWrite()' writes the contents of the given matrix
 * 'm' simply to the standard output.
 * [RETURN VALUE]
 * 'zCMatFWrite()' and 'zCMatWrite()' return no values.
 */
__EXPORT void zCMatFWrite(FILE *fp, zCMat m);
#define zCMatWrite(m) zCMatFWrite( stdout, (m) )

__END_DECLS

#endif /* __ZM_CMAT_H__ */
