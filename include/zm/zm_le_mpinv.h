/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_mpinv - linear equation: Moore=Penrose inverse matrix.
 */

#ifndef __ZM_LE_MPINV_H__
#define __ZM_LE_MPINV_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* METHOD:
 * zMPInv - Moore=Penrose inverse matrix.
 *
 * 'zMPInv()' calculates Moore=Penrose inverse matrix
 * of an arbitrary matrix 'm'. The result is put into
 * 'mp'. Suppose a matrix B is Moore=Penrose inverse
 * matrix of a matrix A, and the following relations
 * are satisfied.
 *   A B A = A
 *   B A B = B
 *   ( A B )^T = A B
 *   ( B A )^T = B A
 * #
 * The computation algorithm is based on Penrose s
 * iterative method.
 * [RETURN VALUE]
 * 'zMPInv()' returns the rank of the original matrix
 * 'm'.
 */
__EXPORT int zMPInv(zMat m, zMat mp);
__EXPORT int zMPInvPenrose(zMat m, zMat mp);

/*! \brief MP inverse with its null space.
 */
__EXPORT int zMPInvNull(zMat m, zMat mp, zMat mn);

__EXPORT zMat zMulMPInvMatMat(zMat m1, zMat m2, zMat m);

__END_DECLS

#endif /* __ZM_LE_MPINV_H__ */
