/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_minv - linear equation: determinant and inverse matrix.
 */

#ifndef __ZM_LE_MINV_H__
#define __ZM_LE_MINV_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* METHOD:
 * zMatDet, zMatDetDST, zMatAdj,
 * zMulInvMatMat, zMulMatInvMat, zMatInv
 * - determinant and inverse of matrix.
 *
 * 'zMatDet()' calculates the determinant value of the
 * matrix 'm'.
 * 'zMatDetDST()' destructively modifies 'm' while
 * calculating the determinant.
 * 'index' is an index vector for calculation.
 * Before calling 'zMatDetDST()', 'index' should be
 * ordered by 'zIndexOrder()'.
 * #
 * 'zMatAdj()' calculates adjoint matrix of 'm'.
 * 'm' must be a square matrix.
 * #
 * 'zMulInvMatMat()' calculates a multiplication of
 * inverse matrix of 'm1' and 'm2'. The result is put into 'm'.
 * i.e. 'm' = 'm1'^-1 'm2'.
 * #
 * 'zMulMatInvMat()' calculates a multiplication of
 * 'm1' and inverse matrix of 'm2'. The result is put into 'm'.
 * i.e. 'm' = 'm1' 'm2'^-1.
 * #
 * 'zMatInv()' calculates the inverse matrix of 'm'
 * and set it into 'im'.
 * [RETURN VALUE]
 * 'zMatDet()' and 'zMatDetDST()' return the
 * determinant calculated.
 * #
 * 'zMatAdj()' returns a pointer 'm'.
 * #
 * 'zMulInvMatMat()' returns a pointer to 'm', if succeeding,
 * or the null vector if 'm1' is not full rank.
 * #
 * 'zMatInv()' returns a pointer to 'im', if succeeding,
 * or the null vector if 'm' is not full rank.
 * [SEE ALSO]
 * zMatInvHotelling
 */
__EXPORT double zMatDetDST(zMat m, zIndex index);
__EXPORT double zMatDet(zMat m);
__EXPORT zMat zMatAdj(zMat m, zMat adj);
__EXPORT zMat zMulInvMatMat(zMat m1, zMat m2, zMat m);
__EXPORT zMat zMulMatInvMat(zMat m1, zMat m2, zMat m);
__EXPORT zMat zMatInv(zMat m, zMat im);

/* METHOD:
 * zMatInvHotelling
 * - inverse matrix by Hotelling's method.
 *
 * 'zMatInvHotelling()' computes the inverse matrix of
 * a given matrix 'm' by Hotelling s method. The result
 * is put into 'im', if it succeeds.
 * 'tol' is the tollerance of iteration.
 * 'iter' is the maximum number of iteration. If zero is
 * given for it, Z_MAX_ITER_NUM defined in 'zm_misc.h'
 * is used instead.
 * [RETURN VALUES]
 * 'zMatInvHotelling()' returns the pointer 'im' if succeeds.
 * If it fails to allocatewhen working memories, or 'm'
 * is a non-square matrix, or the sizes of 'm' and 'im
 * are mismatch, the null pointer is returned.
 * [SEE ALSO]
 * zMatInv
 */
__EXPORT zMat zMatInvHotelling(zMat m, zMat im, double tol, int iter);

__END_DECLS

#endif /* __ZM_LE_MINV_H__ */
