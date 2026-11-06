/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_lu - linear equation: LU decomposition.
 */

#ifndef __ZM_LE_LU_H__
#define __ZM_LE_LU_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief relative tolerance for rank detection in LU decomposition.
 *
 * A column is regarded as linearly dependent on the preceding ones if the maximum absolute value of
 * the candidates of the pivot is not greater than this value times the maximum absolute component of
 * the original matrix.
 */
#define ZM_LU_DECOMP_TOL ( 1.0e-12 )

/*! \brief LU decomposition.
 *
 * zMatDecompLU() decomposes a matrix \a m into a lower triangular matrix \a l and an upper triangular
 * matrix \a u based on Crout method with partial pivoting, namely, \a l \a u = \a m.
 * \a index is an index vector for order description. The diagonal components of \a l are the pivots,
 * which are found in the rows \a index[0], \a index[1], ... of \a l in this order, while the rows of
 * \a l correspond to those of \a m. The components of \a u at the pivots are one.
 * If \a m is with the size of r x c, \a l and \a u have to be with the sizes of r x r and r x c,
 * respectively.
 * If \a m is degenerated, the sizes of \a l and \a u are automatically adjusted to be column-full-rank
 * and row-full-rank, respectively. Namely, suppose the rank of \a m is d, the sizes of \a l and \a u
 * will be adjusted to r x d and d x c, respectively.
 * The rank is detected with the relative tolerance ZM_LU_DECOMP_TOL, so that the result does not
 * depend on the scale of \a m.
 *
 * zMatDecompLUAlloc() automatically allocates matrices \a l and \a u in addition to the index vector
 * \a index, and then conducts the LU decomposition.
 * \return
 * zMatDecompLU() and zMatDecompLUAlloc() return the rank of \a m, which becomes the same with the
 * minimum of the row and column size of \a m when \a m is full rank.
 * If they fail to allocate internal memory for the computation or find mismatches of sizes of matrices,
 * they return -1.
 */
__ZM_EXPORT int zMatDecompLU(const zMat m, zMat l, zMat u, zIndex index);
__ZM_EXPORT int zMatDecompLUAlloc(const zMat m, zMat *l, zMat *u, zIndex *index);

/* Cholesky decomposition.
 *
 * zMatDecompCholesky() decomposes a positive semi-definite symmetric matrix \a m into \a l \a l ^T.
 * \a index is an index vector for order description
 * \a l has to have the same size with \a m.
 * If \a m is degenerated, the size of \a l is automatically adjusted to be column-full-rank.
 *
 * zMatDecompCholeskyDST() also conducts the Cholesky decomposition, during which the given matrix \a m
 * is destroyed.
 *
 * zMatDecompCholeskyAlloc() automatically allocates matrix \a l and the index vector \a index, and then
 * conducts the Cholesky decomposition.
 * \return
 * zMatDecompCholeskyDST(), zMatDecompCholesky(), and zMatDecompCholeskyAlloc() return the rank of \a m.
 * If they fail to allocate internal memory for the computation or find mismatches of sizes of matrices,
 * they return -1.
 */
__ZM_EXPORT int zMatDecompCholeskyDST(zMat m, zMat l, zIndex index);
__ZM_EXPORT int zMatDecompCholesky(const zMat m, zMat l, zIndex index);
__ZM_EXPORT int zMatDecompCholeskyAlloc(const zMat m, zMat *l, zIndex *index);

__END_DECLS

#endif /* __ZM_LE_LU_H__ */
