/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_mpinv - linear equation: Moore=Penrose inverse matrix.
 */

#ifndef __ZM_LE_MPINV_H__
#define __ZM_LE_MPINV_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief Moore=Penrose inverse matrix.
 *
 * zMatMPInv() calculates the Moore=Penrose inverse matrix of an
 * arbitrary matrix \a m. The result is put into \a mp.
 * Suppose a matrix B is the Moore=Penrose inverse matrix of
 * a matrix A. The following relations are satisfied.
 *   A B A = A
 *   B A B = B
 *   ( A B )^T = A B
 *   ( B A )^T = B A
 *
 * The computation algorithm is based on Penrose s iterative method.
 * \return
 * zMatMPInv() returns the rank of the original matrix \a m.
 */
__ZM_EXPORT int zMatMPInv(zMat m, zMat mp);
__ZM_EXPORT int zMatMPInvPenrose(zMat m, zMat mp);

/*! \brief MP inverse with its null space.
 */
__ZM_EXPORT int zMatMPInvNull(zMat m, zMat mp, zMat mn);

__ZM_EXPORT zMat zMulMPInvMatMat(zMat m1, zMat m2, zMat m);

__END_DECLS

#endif /* __ZM_LE_MPINV_H__ */
