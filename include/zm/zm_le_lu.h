/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_lu - linear equation: LU decomposition.
 */

#ifndef __ZM_LE_LU_H__
#define __ZM_LE_LU_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief LU decomposition.
 *
 * zMatDecompLU() decomposes a matrix \a m into a lower triangular matrix \a l and an upper triangular
 * matrix \a u based on Crout method.
 * \a index is an index vector for order description.
 * If \a m is with the size of r x c, \a l and \a u are supposed to have more than or equal to r x s
 * and s x c sizes, where s is the minimum of the row size and the column size of \a m.
 * If \a m is degenerated, the sizes of \a l and \a u are automatically adjusted to be column-full-rank
 * and row-full-rank, respectively.
 *
 * zMatDecompLUDST() also conducts the LU decomposition, during which the given matrix \a m is destroyed.
 *
 * zMatDecompLUAlloc() automatically allocates matrices \a l and \a u in addition to the index vector
 * \a index, and then conducts the LU decomposition.
 * \return
 * zMatDecompLUDST(), zMatDecompLU(), and zMatDecompLUAlloc() return the rank of \a m, which becomes
 * the same with the minimum of the row and column size of \a m when \a m is full rank.
 * If they fail to allocate internal memory for the computation or find mismatches of sizes of matrices,
 * they return -1.
 */
__ZM_EXPORT int zMatDecompLUDST(zMat m, zMat l, zMat u, zIndex index);
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
