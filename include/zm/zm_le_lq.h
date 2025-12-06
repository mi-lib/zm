/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_lq - linear equation: LQ/QR decomposition.
 */

#ifndef __ZM_LE_LQ_H__
#define __ZM_LE_LQ_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief tolerance for rank detection in LQ decomposition based on Gram=Schmidt's method. */
#define ZM_LQ_DECOMP_GRAMSCHMIDT_TOL ( 1.0e-10 )

/*! \brief LQ decomposition of a matrix.
 *
 * zMatDecompLQ_GramSchmidt_DST() decomposes a matrix \a m into a lower triangular matrix \a l and
 * a row-orthonormal matrix \a q based on Gram=Schmidt's orthogonalization method, namely,
 *  \a l \a q = \a m.
 * It destroys \a m during the computation.
 * zMatDecompLQ_GramSchmidt() does LQ decomposition based on Gram=Schmidt's orthogonalization method
 * without destroying \a m.
 *
 * zMatDecompLQ_Householder() does LQ decomposition based on Householder's method without destroying \a m.
 * \note
 * zMatDecompLQ_GramSchmidt_DST() and zMatDecompLQ_GramSchmidt() are aliased to zMatDecompLQDST()
 * and zMatDecompLQ(), respectively, as macros.
 * \return
 * zMatDecompLQ_GramSchmidt_DST(), zMatDecompLQ_GramSchmidt(), and zMatDecompLQ_Householder()
 * return the rank of \a m.
 * \notes
 * When the null pointer is given for \a l, these functions compute only the orthonormal matrix \a q.
 */
__ZM_EXPORT int zMatDecompLQ_GramSchmidt_DST(zMat m, zMat l, zMat q);
__ZM_EXPORT int zMatDecompLQ_GramSchmidt(const zMat m, zMat l, zMat q);
__ZM_EXPORT int zMatDecompLQ_Householder(const zMat m, zMat l, zMat q);

/* aliases */
#define zMatDecompLQDST zMatDecompLQ_GramSchmidt_DST
#define zMatDecompLQ    zMatDecompLQ_GramSchmidt

/*! \brief LQ decomposition with the null-space projector.
 *
 * zMatDecompLQNull() decomposes a matrix \a m into a lower triangular matrix \a l and a row-orthonormal
 * matrix \a q based on Householder's method. It also finds the null-space projector matrix \a qnull.
 * Namely, the following equations are satisfied:
 *  \a l \a q = \a m
 *  \a m \a qnull = O.
 * \return
 * zMatDecompLQNull() returns the rank of \a m.
 */
__ZM_EXPORT int zMatDecompLQNull(const zMat m, zMat l, zMat q, zMat qnull);

/*! \brief LQ decomposition based on Gram=Schmidt's method and resizing of a matrix.
 */
__ZM_EXPORT int zMatDecompLQAndResize(const zMat m, zMat l, zMat q);

/*! \brief LQ decomposition with an automatic matrix allocation.
 */
__ZM_EXPORT int zMatDecompLQAlloc(const zMat m, zMat *l, zMat *q);

/*! \brief QR decomposition of a matrix based on Gram=Schmidt's method.
 *
 * zMatDecompQR() decomposes a matrix \a m into a column-orthonormal matrix \a q and an upper triangular
 * matrix \a u based on Gram=Schmidt's orthogonalization method, namely,
 *  \a q \a r = \a m
 * \return
 * zMatDecompQR() returns the rank of \a m.
 * \notes
 * When the null pointer is given for \a r, zMatDecompQR() computes only the orthonormal matrix \a q.
 */
__ZM_EXPORT int zMatDecompQR(const zMat m, zMat q, zMat r);

__END_DECLS

#endif /* __ZM_LE_LQ_H__ */
