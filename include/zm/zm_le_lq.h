/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_lq - linear equation: LQ/QR decomposition.
 */

#ifndef __ZM_LE_LQ_H__
#define __ZM_LE_LQ_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief relative tolerance for rank detection in LQ decomposition based on Gram=Schmidt's method.
 *
 * A row vector is regarded as linearly dependent on the preceding ones if the norm of its residual
 * after orthogonalization is not greater than this value times the norm of the original row vector.
 */
#define ZM_LQ_DECOMP_GRAMSCHMIDT_TOL ( 1.0e-10 )

/*! \brief relative tolerance for rank detection in LQ decomposition based on Householder's method.
 *
 * A row vector is regarded as linearly dependent on the preceding ones if the norm of its components
 * remaining after the preceding Householder reflections is not greater than this value times the norm
 * of the original row vector.
 */
#define ZM_LQ_DECOMP_HOUSEHOLDER_TOL ( 1.0e-10 )

/*! \brief LQ decomposition of a matrix.
 *
 * zMatDecompLQ() decomposes a matrix \a m into a lower triangular matrix \a l and a row-orthonormal
 * matrix \a q based on modified Gram=Schmidt's orthogonalization method, namely,
 *  \a l \a q = \a m.
 * If \a m is with the size of r x c, \a l and \a q have to be with the sizes of r x s and s x c,
 * respectively, where s is the minimum of r and c.
 * If \a m is degenerated, the sizes of \a l and \a q are automatically adjusted to be column-full-rank
 * and row-full-rank, respectively. Namely, suppose the rank of \a m is d (<= s), the sizes of \a l
 * and \a q will be adjusted to r x d and d x c, respectively.
 * The rank is detected with the relative tolerance ZM_LQ_DECOMP_GRAMSCHMIDT_TOL, so that the result
 * does not depend on the scale of \a m.
 *
 * zMatDecompLQAlloc() automatically allocates matrices \a l and \a q, and then conducts the LQ decomposition.
 * \return
 * zMatDecompLQ() and zMatDecompLQAlloc() return the rank of \a m. If \a m is a zero matrix, they
 * return 0, and the sizes of \a l and \a q are adjusted to r x 0 and 0 x c, respectively.
 * If they fail to allocate internal memory for the computation or find mismatches of sizes of matrices,
 * they return -1.
 */
__ZM_EXPORT int zMatDecompLQ(const zMat m, zMat l, zMat q);
__ZM_EXPORT int zMatDecompLQAlloc(const zMat m, zMat *l, zMat *q);

/*! \brief full-sized LQ decomposition of a matirx.
 *
 * zMatDecompLQFull() decomposes a matrix \a m into a lower triangular matrix \a l and a row-orthonormal
 * matrix \a q based on Householder's method, namely,
 *  \a l \a q = \a m.
 * If \a m is with the size of r x c, \a l and \a q have to be with the sizes of r x c and c x c, respectively.
 * Row vectors of \a m are processed in order, and a row vector linearly dependent on the preceding ones
 * (judged with the relative tolerance ZM_LQ_DECOMP_HOUSEHOLDER_TOL) does not consume a column of \a l.
 * Hence, suppose the rank of \a m is d, the last c-d columns of \a l are zero, and the last c-d rows of
 * \a q span the null space of \a m.
 *
 * zMatDecompLQFullAlloc() automatically allocates matrices \a l and \a q, and then conducts the full-sized
 * LQ decomposition.
 * \return
 * zMatDecompLQFull() and zMatDecompLQFullAlloc() return the rank of the matrix \a m, or -1 if they fail
 * to allocate internal memory or find mismathces of sizes of matrices.
 */
__ZM_EXPORT int zMatDecompLQFull(const zMat m, zMat l, zMat q);
__ZM_EXPORT int zMatDecompLQFullAlloc(const zMat m, zMat *l, zMat *q);

/*! \brief LQ decomposition with the null-space projector.
 *
 * zMatDecompLQNull() decomposes a matrix \a m into a lower triangular matrix \a l and a row-orthonormal
 * matrix \a q based on Householder's method. It also finds the null-space projector matrix \a qnull.
 * Namely, the following equations are satisfied:
 *  \a l \a q = \a m
 *  \a m \a qnull = O.
 * If \a m is with the size of r x c, \a l, \a q, and \a qnull have to be with the sizes r x c, c x c, and
 * c x c, respectively, while their sizes will be adjusted to r x d, d x c, and c x (c-d), respectively,
 * where d is the rank of \a m. If \a m is a zero matrix, \a qnull is the identity matrix.
 * \return
 * zMatDecompLQNull() returns the rank of \a m, or -1 if it finds mismatches of sizes of matrices.
 */
__ZM_EXPORT int zMatDecompLQNull(const zMat m, zMat l, zMat q, zMat qnull);

/*! \brief QR decomposition of a matrix based on Gram=Schmidt's method.
 *
 * zMatDecompQR() decomposes a matrix \a m into a column-orthonormal matrix \a q and an upper triangular
 * matrix \a r based on Gram=Schmidt's orthogonalization method, namely,
 *  \a q \a r = \a m
 * If \a m is with the size of r x c, \a q and \a r are supposed to have more than or equal to r x s
 * and s x c sizes, respectively, where s is the minimum of the row size and the column size of \a m.
 * If \a m is degenerated, the sizes of \a q and \a r are automatically adjusted to be column-full-rank
 * and row-full-rank, respectively.
 * \return
 * zMatDecompQR() returns the rank of \a m.
 */
__ZM_EXPORT int zMatDecompQR(const zMat m, zMat q, zMat r);

__END_DECLS

#endif /* __ZM_LE_LQ_H__ */
