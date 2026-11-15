/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_mat_decomp_lu - matrix decomposition: LU decomposition.
 */

#ifndef __ZM_MAT_DECOMP_LU_H__
#define __ZM_MAT_DECOMP_LU_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief relative tolerance for rank detection in LU decomposition.
 *
 * A pivot is regarded as zero, namely, the remaining part of the matrix is regarded as linearly
 * dependent on the part already decomposed, if the maximum absolute value of the candidates of the
 * pivot is not greater than this value times the maximum absolute component of the original matrix.
 */
#define ZM_LU_DECOMP_TOL ( 1.0e-12 )

/*! \brief LU decomposition.
 *
 * zMatDecompLU() and zMatDecompLUCP() decompose a matrix \a m into a lower triangular matrix \a l and
 * an upper triangular matrix \a u based on Crout method, namely, \a l \a u = \a m. The components of
 * \a u at the pivots are one.
 * If \a m is with the size of r x c, \a l and \a u have to be with the sizes of r x s and s x c,
 * respectively, where s is the minimum of r and c.
 * If \a m is degenerated, the sizes of \a l and \a u are automatically adjusted to be column-full-rank
 * and row-full-rank, respectively. Namely, suppose the rank of \a m is d (<= s), the sizes of \a l
 * and \a u will be adjusted to r x d and d x c, respectively.
 * The rank is detected with the relative tolerance ZM_LU_DECOMP_TOL, so that the result does not
 * depend on the scale of \a m.
 * They do not allocate any internal memory.
 *
 * The two functions differ in the pivoting strategy.
 *
 * zMatDecompLU() is based on partial pivoting, namely, the columns of \a m are processed in the
 * original order, and the component with the largest absolute value in each column among the remaining
 * rows is chosen as the pivot.
 * \a index is an index vector for order description of rows. The diagonal components of \a l are the
 * pivots, which are found in the rows \a index[0], \a index[1], ... of \a l in this order, while the
 * rows of \a l correspond to those of \a m. The size of \a index has to be the same with the row size
 * of \a m.
 * It is suitable for solving linear equations, but could overestimate the rank of \a m if \a m has a
 * column nearly dependent on the preceding ones. Namely, such a column is chosen as a pivot because of
 * its small residual, which amplifies the rounding errors so that the following dependent columns look
 * independent.
 *
 * zMatDecompLUCP() is based on complete pivoting, namely, the component with the largest absolute value
 * among the remaining rows and columns is chosen as the pivot at each step, so that the rank of \a m is
 * detected reliably even in the above case. Use it when the rank matters.
 * \a index_row and \a index_col are index vectors for order description of rows and columns,
 * respectively. The pivots are found at the (\a index_row[0], 0), (\a index_row[1], 1), ... components
 * of \a l, and at the (0, \a index_col[0]), (1, \a index_col[1]), ... components of \a u in this order.
 * The rows of \a l correspond to those of \a m, and the columns of \a u correspond to those of \a m.
 * Namely, \a l is a lower triangular matrix with rows permuted by \a index_row, and \a u is an upper
 * triangular matrix with columns permuted by \a index_col. The sizes of \a index_row and \a index_col
 * have to be the same with the row and column sizes of \a m, respectively.
 *
 * zMatDecompLUAlloc() and zMatDecompLUCPAlloc() automatically allocate matrices \a l and \a u in
 * addition to the index vectors, and then conduct zMatDecompLU() and zMatDecompLUCP(), respectively.
 * \return
 * These functions return the rank of \a m, which becomes the same with the minimum of the row and
 * column sizes of \a m when \a m is full rank.
 * If they find mismatches of sizes of matrices and index vectors, or zMatDecompLUAlloc() and
 * zMatDecompLUCPAlloc() fail to allocate memory, they return -1.
 * \sa
 * zLESolveL, zLESolveU, zLESolveUCP
 */
__ZM_EXPORT int zMatDecompLU(const zMat m, zMat l, zMat u, zIndex index);
__ZM_EXPORT int zMatDecompLUAlloc(const zMat m, zMat *l, zMat *u, zIndex *index);
__ZM_EXPORT int zMatDecompLUCP(const zMat m, zMat l, zMat u, zIndex index_row, zIndex index_col);
__ZM_EXPORT int zMatDecompLUCPAlloc(const zMat m, zMat *l, zMat *u, zIndex *index_row, zIndex *index_col);

/*! \brief relative tolerance for rank detection in Cholesky decomposition.
 *
 * A pivot is regarded as zero if its absolute value is not greater than this value times the maximum
 * absolute diagonal component of the original matrix. A negative pivot beyond it means that the matrix
 * is not positive semi-definite.
 */
#define ZM_CHOLESKY_DECOMP_TOL ( 1.0e-12 )

/*! \brief Cholesky decomposition.
 *
 * zMatDecompCholesky() decomposes a positive semi-definite symmetric matrix \a m into \a l \a l ^T.
 * \a index is an index vector for order description. The rows of \a l correspond to those of \a m,
 * and the rows \a index[0], \a index[1], ... of \a l are lower triangular in this order.
 * \a l has to have the same size with \a m.
 * If \a m is degenerated, the size of \a l is automatically adjusted to be column-full-rank.
 * The rank is detected with the relative tolerance ZM_CHOLESKY_DECOMP_TOL, so that the result does
 * not depend on the scale of \a m.
 *
 * zMatDecompCholeskyAlloc() automatically allocates matrix \a l and the index vector \a index, and then
 * conducts the Cholesky decomposition. If it fails to allocate them, both \a l and \a index are set
 * for the null pointer.
 * \return
 * zMatDecompCholesky() and zMatDecompCholeskyAlloc() return the rank of \a m.
 * If they fail to allocate internal memory for the computation, find mismatches of sizes of matrices,
 * or find that \a m is not positive semi-definite, they return -1.
 */
__ZM_EXPORT int zMatDecompCholesky(const zMat m, zMat l, zIndex index);
__ZM_EXPORT int zMatDecompCholeskyAlloc(const zMat m, zMat *l, zIndex *index);

__END_DECLS

#endif /* __ZM_MAT_DECOMP_LU_H__ */
