/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_lq - linear equation: LQ/QR decomposition.
 */

#ifndef __ZM_LE_LQ_H__
#define __ZM_LE_LQ_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* METHOD:
 * zLQDecompDST, zLQDecomp, zLQDecompReg, zQRDecomp
 * - LQ/QR decomposition based on Gram=Schmidt's method.
 *
 * 'zLQDecompDST()' decompose given matrix 'mat' into
 * a lower triangular matrix 'l' and a normalized
 * orthogonal matrix 'q', based on Gram=Schmidt s
 * orthogonalization method, namely,
 *  'l' 'q' = 'mat'.
 * It destroys 'mat' in the course of computation.
 * #
 * 'zLQDecomp()' does LQ decomposition without destroying
 * 'mat'.
 * #
 * 'zQRDecomp()' does QR decomposition, a transpose
 * of 'zLQDecomp()', namely,
 *  'q' 'r' = 'm'
 * [RETURN VALUE]
 * 'zLQDecompDST()' and 'zLQDecomp()' return column-rank.
 * 'zQRDecomp()' returns row-rank.
 * [NOTES]
 * The null pointer is assignable for 'l' (or 'r').
 * When the null pointer is given for them, these
 * functions just compute an orthogonal matrix 'q'.
 * #
 * If one column / row or more are dependent on
 * other columns / rows, they will fail to decompose
 * the given matrix into an orthogonal space.
 */
__EXPORT int zLQDecompDST(zMat m, zMat l, zMat q, zIndex idx);
__EXPORT int zLQDecomp(zMat m, zMat l, zMat q, zIndex idx);
__EXPORT int zLQDecompReg(zMat m, zMat l, zMat q, zIndex idx);

/*! \brief LQ decomposition with an automatic matrix allocation and resize.
 */
__EXPORT int zLQDecompAlloc(zMat m, zMat *l, zMat *q, zIndex *idx);

__EXPORT int zQRDecomp(zMat m, zMat q, zMat r, zIndex idx);

__END_DECLS

#endif /* __ZM_LE_LQ_H__ */
