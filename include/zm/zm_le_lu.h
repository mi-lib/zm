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
 * zMatDecompLU() decomposes the matrix \a m into a lower triangular matrix \a l and an upper triangular
 * matrix \a u based on Crout method.
 * \a index is an index vector for order discription.
 * If \a m is degenerated, the sizes of \a l and \a u are automatically adjusted.
 *
 * zMatDecompLUDST() destroys \a m during the LU decomposition.
 * \return
 * zMatDecompLUDST() and zMatDecompLU() return the rank of \a m, which becomes the same with the minimum
 * of the row and column size of \a m when \a m is full rank.
 */
__ZM_EXPORT int zMatDecompLUDST(zMat m, zMat l, zMat u, zIndex idx);
__ZM_EXPORT int zMatDecompLU(const zMat m, zMat l, zMat u, zIndex idx);
__ZM_EXPORT int zMatDecompLUAndResize(const zMat m, zMat l, zMat u, zIndex idx);
__ZM_EXPORT int zMatDecompLUAlloc(const zMat m, zMat *l, zMat *u, zIndex *idx);

/* Cholesky decomposition.
 *
 * zMatDecompCholesky() decomposes a positive semi-definite symmetric matrix \a m into \a l \a l ^T.
 *
 * zMatDecompCholeskyDST() destroys \a m during the Cholesky decomposition.
 * \return
 * zMatDecompCholeskyDST() and zMatDecompCholesky() return the rank of \a m. If they fail to allocate
 * internal memory for the computation, it returns -1.
 */
__ZM_EXPORT int zMatDecompCholeskyDST(zMat m, zMat l, zIndex index);
__ZM_EXPORT int zMatDecompCholesky(const zMat m, zMat l, zIndex index);
__ZM_EXPORT int zMatDecompCholeskyAndResize(const zMat m, zMat l, zIndex idx);
__ZM_EXPORT int zMatDecompCholeskyAlloc(const zMat m, zMat *l, zIndex *idx);

__END_DECLS

#endif /* __ZM_LE_LU_H__ */
