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
 * zLUDecomp() decomposes the matrix \a m into a lower
 * triangular matrix \a l and an upper triangular matrix \a u
 * with Crout method.
 * \a index is an index vector for order discription.
 *
 * zLUDecompReg() regresses \a l and \a u after decomposing
 * \a m, if the rank of \a m less than its size.
 *
 * zLUDecompDST() destroys \a m during decomposition.
 * \return
 * zLUDecompDST(), zLUDecomp() and zLUDecompReg() return the
 * rank of \a m, which becomes the same with the minimum of
 * the row and column size of \a m when \a m is full rank.
 */
__EXPORT int zLUDecompDST(zMat m, zMat l, zMat u, zIndex idx);
__EXPORT int zLUDecomp(zMat m, zMat l, zMat u, zIndex idx);
__EXPORT int zLUDecompReg(zMat m, zMat l, zMat u, zIndex idx);
__EXPORT int zLUDecompAlloc(zMat m, zMat *l, zMat *u, zIndex *idx);

/* Cholesky decomposition.
 *
 * zCholeskyDecomp() decomposes a symmetric matrix \a m into
 * \a l \a l ^T.
 *
 * zCholeskyDecompDST() destroys \a m during decomposition.
 * \return
 * zCholeskyDecompDST() and zCholeskyDecomp() return the rank
 * of \a m. In case they fail to allocate internal memory for
 * the computation, -1 is returned.
 */
__EXPORT int zCholeskyDecompDST(zMat m, zMat l, zIndex index);
__EXPORT int zCholeskyDecomp(zMat m, zMat l, zIndex index);
__EXPORT int zCholeskyDecompReg(zMat m, zMat l, zIndex idx);
__EXPORT int zCholeskyDecompAlloc(zMat m, zMat *l, zIndex *idx);

__END_DECLS

#endif /* __ZM_LE_LU_H__ */
