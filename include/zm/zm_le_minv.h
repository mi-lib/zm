/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_minv - linear equation: determinant and inverse matrix.
 */

#ifndef __ZM_LE_MINV_H__
#define __ZM_LE_MINV_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief determinant and inverse of matrix.
 *
 * zMatDet() calculates the determinant value of the matrix \a m.
 * zMatDetDST() destructively modifies \a m when calculating the
 * determinant.
 * \a index is an index vector for calculation, which has to be
 * ordered by zIndexOrder() before calling zMatDetDST().
 *
 * zMatAdj() calculates the adjoint matrix of \a m.
 * \a m must be a square matrix.
 * The result is put into \a adj.
 *
 * zMulInvMatMat() calculates a multiplication of the inverse matrix
 * of \a m1 and \a m2. The result is put into \a m.
 * i.e. \a m = \a m1^-1 \a m2
 *
 * zMulMatInvMat() calculates a multiplication of \a m1 and the
 * inverse matrix of \a m2. The result is put into \a m.
 * i.e. \a m = \a m1 \a m2^-1
 *
 * zMatInv() calculates the inverse matrix of \a m and set it into
 * \a im.
 * \return
 * zMatDet() and zMatDetDST() return the determinant of the given
 * matrix.
 *
 * zMatAdj() returns a pointer \a adj.
 *
 * zMulInvMatMat() returns a pointer \a m, if succeeding.
 * The null pointer is returned if \a m1 is not full rank.
 *
 * zMatInv() returns a pointer \a im, if succeeding.
 * The null pointer is returned if \a m is not full rank.
 * \sa
 * zMatInvHotelling
 */
__EXPORT double zMatDetDST(zMat m, zIndex index);
__EXPORT double zMatDet(zMat m);
__EXPORT zMat zMatAdj(zMat m, zMat adj);
__EXPORT zMat zMulInvMatMat(zMat m1, zMat m2, zMat m);
__EXPORT zMat zMulMatInvMat(zMat m1, zMat m2, zMat m);
__EXPORT zMat zMatInv(zMat m, zMat im);

/*! \brief inverse matrix by Hotelling's method.
 *
 * zMatInvHotelling() computes the inverse matrix of a given matrix
 * \a m by Hotelling's method. The result is put into \a im, if it
 * succeeds.
 * \a tol is the tollerance of iteration.
 * \a iter is the maximum number of iteration. If zero is given for
 * it, Z_MAX_ITER_NUM defined in "zm_misc.h" is used instead.
 * \return
 * zMatInvHotelling() returns the pointer \a im if succeeds.
 * If it fails to allocate for working memories, or \a m is not a
 * square matrix, or the sizes of \a m'and \a im do not match, the
 * null pointer is returned.
 * \sa
 * zMatInv
 */
__EXPORT zMat zMatInvHotelling(zMat m, zMat im, double tol, int iter);

__END_DECLS

#endif /* __ZM_LE_MINV_H__ */
