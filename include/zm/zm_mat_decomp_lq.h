/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_mat_decomp_lq - matrix decomposition: LQ/QR decomposition.
 */

#ifndef __ZM_MAT_DECOMP_LQ_H__
#define __ZM_MAT_DECOMP_LQ_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* NOTE: the following functions with the prefix of underscore are for internal use in the library,
 * and not supposed to be used in user programs. Their specifications may be changed without notice. */

/*! \brief exactly scale a matrix by a power of two (internal use).
 *
 * _zMatScalePow2() directly divides a matrix \a m by the power of two with which the absolute maximum
 * component of \a m is in [1,2). It is used to avoid overflow and underflow in matrix decompositions.
 * \return
 * _zMatScalePow2() returns the power of two with which \a m is divided, or 1 if \a m is a zero matrix.
 */
__ZM_EXPORT double _zMatScalePow2(zMat m);

/*! \brief Householder reflection (internal use).
 *
 * _zMatHouseholderUnitNormal() computes the unit normal vector \a nu of the Householder reflection
 * I - 2 \a nu \a nu^T that maps an \a n-dimensional vector \a x to \a s e1, where |\a s|=|\a x|.
 * The sign of \a s is chosen to be opposite to that of the first component of \a x in order to avoid
 * cancellation. \a x is scaled by its absolute maximum component to avoid overflow and underflow.
 * \a x and \a nu must not overlap.
 *
 * _zMatHouseholderL() applies the reflection to rows from the \a r0-th to the (\a r0+\a n-1)-th of a
 * matrix \a m from the left side, where only columns from the \a c0-th to the \a c1-th are modified.
 * \a w is a workspace with the column size of \a m.
 *
 * _zMatHouseholderR() applies the reflection to columns from the \a c0-th to the (\a c0+\a n-1)-th of a
 * matrix \a m from the right side, where only rows from the \a r0-th to the \a r1-th are modified.
 * \return
 * _zMatHouseholderUnitNormal() returns the false value if \a x is the zero vector, for which no
 * reflection is needed. Otherwise, it returns the true value.
 * _zMatHouseholderL() and _zMatHouseholderR() return no value.
 */
__ZM_EXPORT bool _zMatHouseholderUnitNormal(const double *x, int n, double *nu, double *s);
__ZM_EXPORT void _zMatHouseholderL(zMat m, const double *nu, int n, int r0, int c0, int c1, double *w);
__ZM_EXPORT void _zMatHouseholderR(zMat m, const double *nu, int n, int c0, int r0, int r1);

/*! \brief LQ decomposition of a matrix.
 *
 * zMatDecompLQ() decomposes a matrix \a m into a lower triangular matrix \a l and a row-orthonormal
 * matrix \a q with row pivoting, namely,
 *  \a l \a q = \a m.
 * \a index is an index vector for order description. The row vectors of \a m are orthogonalized in the
 * order of the norms of their residuals, and the diagonal components of \a l are found in the rows
 * \a index[0], \a index[1], ... of \a l in this order, while the rows of \a l correspond to those of
 * \a m. Namely, \a l is a lower triangular matrix with rows permuted by \a index. The size of \a index
 * has to be the same with the row size of \a m.
 * If \a m is with the size of r x c, \a l and \a q have to be with the sizes of r x s and s x c,
 * respectively, where s is the minimum of r and c.
 * If r is not greater than c, it is based on modified Gram=Schmidt's orthogonalization method. Otherwise,
 * it is based on Householder's method, where \a l and \a q are used as workspaces.
 * If \a m is degenerated, the sizes of \a l and \a q are automatically adjusted to be column-full-rank
 * and row-full-rank, respectively. Namely, suppose the rank of \a m is d (<= s), the sizes of \a l
 * and \a q will be adjusted to r x d and d x c, respectively.
 * In the course of the orthogonalization, the remaining row vectors are regarded as linearly dependent
 * on the preceding ones if the maximum norm of their residuals after orthogonalization (by Gram=Schmidt's
 * method) or of their components remaining after the preceding reflections (by Householder's method) is
 * not greater than zTOL times the maximum norm of the row vectors of \a m, so that the rank does not
 * depend on the scale of \a m.
 * It does not allocate any internal memory.
 *
 * zMatDecompLQAlloc() automatically allocates matrices \a l and \a q in addition to the index vector
 * \a index, and then conducts the LQ decomposition.
 * \return
 * zMatDecompLQ() and zMatDecompLQAlloc() return the rank of \a m. If \a m is a zero matrix, they
 * return 0, and the sizes of \a l and \a q are adjusted to r x 0 and 0 x c, respectively.
 * If they fail to allocate memory or find mismatches of sizes of matrices and the index vector,
 * they return -1.
 */
__ZM_EXPORT int zMatDecompLQ(const zMat m, zMat l, zMat q, zIndex index);
__ZM_EXPORT int zMatDecompLQAlloc(const zMat m, zMat *l, zMat *q, zIndex *index);

/*! \brief full-sized LQ decomposition of a matirx.
 *
 * zMatDecompLQFull() decomposes a matrix \a m into a lower triangular matrix \a l and a row-orthonormal
 * matrix \a q based on Householder's method with row pivoting, namely,
 *  \a l \a q = \a m.
 * If \a m is with the size of r x c, \a l and \a q have to be with the sizes of r x c and c x c, respectively.
 * \a index is an index vector for order description, which has the same meaning with that for
 * zMatDecompLQ(). The size of \a index has to be r.
 * A row vector linearly dependent on the preceding ones (judged in the same way with zMatDecompLQ())
 * does not consume a column of \a l. Hence, suppose the rank of \a m is d, the last
 * c-d columns of \a l are zero, and the last c-d rows of \a q span the null space of \a m.
 * It does not allocate any internal memory.
 *
 * zMatDecompLQFullAlloc() automatically allocates matrices \a l and \a q in addition to the index vector
 * \a index, and then conducts the full-sized LQ decomposition.
 * \return
 * zMatDecompLQFull() and zMatDecompLQFullAlloc() return the rank of the matrix \a m, or -1 if they fail
 * to allocate memory or find mismathces of sizes of matrices and the index vector.
 */
__ZM_EXPORT int zMatDecompLQFull(const zMat m, zMat l, zMat q, zIndex index);
__ZM_EXPORT int zMatDecompLQFullAlloc(const zMat m, zMat *l, zMat *q, zIndex *index);

/*! \brief LQ decomposition with the null-space projector.
 *
 * zMatDecompLQNull() decomposes a matrix \a m into a lower triangular matrix \a l and a row-orthonormal
 * matrix \a q based on Householder's method with row pivoting. It also finds the null-space projector
 * matrix \a qnull. Namely, the following equations are satisfied:
 *  \a l \a q = \a m
 *  \a m \a qnull = O.
 * If \a m is with the size of r x c, \a l, \a q, and \a qnull have to be with the sizes r x c, c x c, and
 * c x c, respectively, while their sizes will be adjusted to r x d, d x c, and c x (c-d), respectively,
 * where d is the rank of \a m. If \a m is a zero matrix, \a qnull is the identity matrix.
 * \a index is an index vector for order description, which has the same meaning with that for
 * zMatDecompLQ(). The size of \a index has to be r.
 * \return
 * zMatDecompLQNull() returns the rank of \a m, or -1 if it finds mismatches of sizes of matrices and
 * the index vector.
 */
__ZM_EXPORT int zMatDecompLQNull(const zMat m, zMat l, zMat q, zMat qnull, zIndex index);

/*! \brief QR decomposition of a matrix.
 *
 * zMatDecompQR() decomposes a matrix \a m into a column-orthonormal matrix \a q and an upper triangular
 * matrix \a r with column pivoting, namely,
 *  \a q \a r = \a m
 * It is computed as the transpose of the LQ decomposition of the transpose of \a m by zMatDecompLQ().
 * \a index is an index vector for order description. The diagonal components of \a r are found in the
 * columns \a index[0], \a index[1], ... of \a r in this order, while the columns of \a r correspond to
 * those of \a m. Namely, \a r is an upper triangular matrix with columns permuted by \a index. The size
 * of \a index has to be the same with the column size of \a m.
 * If \a m is with the size of r x c, \a q and \a r are supposed to have more than or equal to r x s
 * and s x c sizes, respectively, where s is the minimum of the row size and the column size of \a m.
 * If \a m is degenerated, the sizes of \a q and \a r are automatically adjusted to be column-full-rank
 * and row-full-rank, respectively.
 * \return
 * zMatDecompQR() returns the rank of \a m, or -1 if it fails to allocate internal memory or finds
 * mismatches of sizes of matrices and the index vector.
 */
__ZM_EXPORT int zMatDecompQR(const zMat m, zMat q, zMat r, zIndex index);

__END_DECLS

#endif /* __ZM_MAT_DECOMP_LQ_H__ */
