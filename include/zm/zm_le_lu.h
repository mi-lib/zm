/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_lu - linear equation: LU decomposition.
 */

#ifndef __ZM_LE_LU_H__
#define __ZM_LE_LU_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* METHOD:
 * zLUDecomp, zLUDecompReg, zLUDecompDST
 * - LU decomposition.
 *
 * 'zLUDecomp()' decomposes the matrix 'm' into a
 * lower triangular matrix(L) 'l' and upper triangular
 * matrix(U) 'um' with Crout method.
 * 'index' is an index vector for order discription.
 * #
 * 'zLUDecompReg()' regresses L and U, if the rank of 'm'
 * less than its size.
 * #
 * 'zLUDecompDST()' destroys 'm' during decomposition.
 * [RETURN VALUE]
 * 'zLUDecompDST()', 'zLUDecomp()' and 'zLUDecompReg()'
 * return the rank of 'm', which becomes the same with
 * the minimum of the row and column size of 'm' when
 * 'm' is full rank.
 */
__EXPORT int zLUDecompDST(zMat m, zMat l, zMat u, zIndex idx);
__EXPORT int zLUDecomp(zMat m, zMat l, zMat u, zIndex idx);
__EXPORT int zLUDecompReg(zMat m, zMat l, zMat u, zIndex idx);
__EXPORT int zLUDecompAlloc(zMat m, zMat *l, zMat *u, zIndex *idx);

/* METHOD:
 * zCholeskyDecomp, zCholeskyDecompDST
 * - Cholesky decomposition.
 *
 * zCholeskyDecomp() decomposes a symmetric matrix \a m
 * into \a l \a l ^T.
 *
 * 'zCholeskyDecompDST()' destroys 'm' during decomposition.
 * [RETURN VALUE]
 * 'zCholeskyDecompDST()' and 'zCholeskyDecomp()'
 * return the true value if succeeding to get both 'lmat'
 * and 'umat'. When 'm' is not full rank, they fail to
 * calculate and return the false value.
 * [NOTES]
 * 'zCholeskyDecomp()' is to be treated so as to
 * return the rank of 'm' as much as 'zLUDecomp()'.
 */
__EXPORT int zCholeskyDecompDST(zMat m, zMat l, zIndex index);
__EXPORT int zCholeskyDecomp(zMat m, zMat l, zIndex index);
__EXPORT int zCholeskyDecompReg(zMat m, zMat l, zIndex idx);
__EXPORT int zCholeskyDecompAlloc(zMat m, zMat *l, zIndex *idx);

__END_DECLS

#endif /* __ZM_LE_LU_H__ */
