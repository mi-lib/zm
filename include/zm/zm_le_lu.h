/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_lu - linear equation: LU decomposition.
 */

#ifndef __ZM_LE_LU_H__
#define __ZM_LE_LU_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* LU decomposition.
 *
 * zMatDecompLU() decomposes the matrix \a m into a lower
 * triangular matrix \a l and an upper triangular matrix \a u
 * with Crout method.
 * \a index is an index vector for order discription.
 *
 * zMatDecompLUReg() regresses \a l and \a u after decomposing
 * \a m, if the rank of \a m less than its size.
 *
 * zMatDecompLUDST() destroys \a m during decomposition.
 * \return
 * zMatDecompLUDST(), zMatDecompLU() and zMatDecompLUReg() return the
 * rank of \a m, which becomes the same with the minimum of
 * the row and column size of \a m when \a m is full rank.
 */
__EXPORT int zMatDecompLUDST(zMat m, zMat l, zMat u, zIndex idx);
__EXPORT int zMatDecompLU(zMat m, zMat l, zMat u, zIndex idx);
__EXPORT int zMatDecompLUReg(zMat m, zMat l, zMat u, zIndex idx);
__EXPORT int zMatDecompLUAlloc(zMat m, zMat *l, zMat *u, zIndex *idx);

/* Cholesky decomposition.
 *
 * zMatDecompCholesky() decomposes a symmetric matrix \a m into
 * \a l \a l ^T.
 *
 * zMatDecompCholeskyDST() destroys \a m during decomposition.
 * \return
 * zMatDecompCholeskyDST() and zMatDecompCholesky() return the rank
 * of \a m. In case they fail to allocate internal memory for
 * the computation, -1 is returned.
 */
__EXPORT int zMatDecompCholeskyDST(zMat m, zMat l, zIndex index);
__EXPORT int zMatDecompCholesky(zMat m, zMat l, zIndex index);
__EXPORT int zMatDecompCholeskyReg(zMat m, zMat l, zIndex idx);
__EXPORT int zMatDecompCholeskyAlloc(zMat m, zMat *l, zIndex *idx);

__END_DECLS

#endif /* __ZM_LE_LU_H__ */
