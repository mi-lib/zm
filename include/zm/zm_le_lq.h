/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_lq - linear equation: LQ/QR decomposition.
 */

#ifndef __ZM_LE_LQ_H__
#define __ZM_LE_LQ_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief LQ/QR decomposition based on Gram=Schmidt's method.
 *
 * zMatDecompLQDST() decompose given matrix \a m into a lower triangular
 * matrix \a l and a normalized orthogonal matrix \a q, based on Gram=Schmidt's
 * orthogonalization method, namely,
 *  \a l \a q = \a m.
 * It destroys \a m during the computation.
 *
 * zMatDecompLQ() does LQ decomposition without destroying \a m.
 *
 * zQRDecomp() does QR decomposition, a transpose of \a zMatDecompLQ(), namely,
 *  \a q \a r = \a m
 * \return
 * zMatDecompLQDST() and zMatDecompLQ() return column-rank of \a m.
 * zQRDecomp() returns row-rank of \a m.
 * \notes
 * The null pointer is assignable for \a l or \a r.
 * When the null pointer is given for them, these functions compute only the
 * orthogonal matrix \a q.
 *
 * If one column / row or more are dependent on other columns / rows, those
 * functions fail to decompose the given matrix into an orthogonal space.
 */
__EXPORT int zMatDecompLQDST(zMat m, zMat l, zMat q, zIndex idx);
__EXPORT int zMatDecompLQ(zMat m, zMat l, zMat q, zIndex idx);
__EXPORT int zMatDecompLQReg(zMat m, zMat l, zMat q, zIndex idx);

/*! \brief LQ decomposition with an automatic matrix allocation and resize.
 */
__EXPORT int zMatDecompLQAlloc(zMat m, zMat *l, zMat *q, zIndex *idx);

__EXPORT int zMatDecompQR(zMat m, zMat q, zMat r, zIndex idx);

__END_DECLS

#endif /* __ZM_LE_LQ_H__ */
