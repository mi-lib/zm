/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le - linear equation.
 */

#ifndef __ZM_LE_H__
#define __ZM_LE_H__

#include <zm/zm_mat.h>

__BEGIN_DECLS

/*! \brief make a pair of matrix and vector balanced.
 *
 * zBalancingColDST() destructively makes a matrix column-
 * balanced. Namely, each column of \a m is divided by the
 * absolute-maximum component in the same column.
 *
 * zBalancingDST() destructively makes a pair of matrix
 * \a m and vector \a v balanced. Namely, each column of
 * \a m is primarily divided by the absolute-maximum
 * component in the same column, and then each row of
 * \a m and the corresponding component of \a v by the
 * absolute-maximum value in the row.
 *
 * Scalings are skipped when the balancing factor with
 * respect to the working column or row is zero.
 *
 * zBalancing() makes \a morg and \a vorg balanced, and
 * puts the result into \a m and \a v.
 *
 * They work as preprocesses of a linear equation solver.
 * Actually, zLESolveGaussDST() internally calls
 * zBalancingDST().
 *
 * The column-balancing factors are stored in \a s for all
 * the functions if it is not the null pointer. It is used
 * for re-balancing. For the original equation \a morg x = \a vorg,
 * the balanced equation is \a m \a s x = v. Thus, the solution
 * \a y of \a m y = v has to be amplified by \a s in such a
 * way as zVecAmpDRC( y, s ).
 * \return
 * zBalancingColDST() and zBalancingDST() return no values.
 *
 * zBalancing() returns the false value if size-mismatch
 * happens among arguments. Otherwise, it returns the true
 * value.
 * \sa
 * zLESolveGaussDST
 */
__EXPORT void zBalancingColDST(zMat m, zVec s);
__EXPORT void zBalancingDST(zMat m, zVec v, zVec s);
__EXPORT bool zBalancing(zMat morg, zVec vorg, zMat m, zVec v, zVec s);

/*! \brief residual b - a x.
 *
 * \notes Sizes are not checked. It is not recommended to use
 * this function in user programs.
 */
__EXPORT zVec zLEResidual(zMat a, zVec b, zVec x, zVec res);

/*! \brief linear equation solver by Gauss elimination method.
 *
 * zLESolveGauss() solves linear equation \a a x = \a b based
 * on Gaussian elimination.
 * The answer is put into \a ans.
 * \a s is used for column-balancing, if it is not the null
 * pointer.
 *
 * zLESolveGaussDST() destructively modifies \a a and \a b
 * while calculating.
 * \return
 * zLESolveGaussDST() returns a pointer \a ans if succeeds.
 * When \a a is a singular matrix, the equation does not have
 * a unique answer and the null pointer is returned.
 */
__EXPORT zVec zLESolveGaussDST(zMat a, zVec b, zVec ans, zIndex index, zVec s);
__EXPORT zVec zLESolveGauss(zMat a, zVec b, zVec ans);

/*! \brief linear equation solver by LU decomposition method.
 *
 * zLESolveLU() solves linear equation \a l \a u x = \a b based
 * on LU decomposition, where \a l and \a u are supposed to be
 * a lower and an upper triangular matrices, respectively.
 * The answer is put into \a ans.
 * \a index is an index vector for order discription.
 * \a l y = \a b and \a u x = y are solved by zLESolveL() and
 * zLESolveU(), respectively.
 *
 * zLESolveL() needs the pivot index since zLUDecomp() pivots
 * the original coefficient matrix.
 * \return
 * zLESolveLU(), zLESolveL() and zLESolveU() return a pointer
 * \a ans, if succeeding. Otherwise, they return the null
 * pointer if the original matrix is singular.
 * \sa
 * zLUDecomp
 */
__EXPORT zVec zLESolveL(zMat l, zVec b, zVec ans, zIndex idx);
__EXPORT zVec zLESolveU(zMat u, zVec b, zVec ans);
__EXPORT zVec zLESolveLU(zMat l, zMat u, zVec b, zVec ans, zIndex index);

/*! \brief linear equation solver by residual iteration.
 *
 * zLESolveRI() solves linear equation \a a x = \a b based on
 * LU decomposition and improves the accuracy by residual
 * iteration. The answer is put into \a ans.
 * \return \a ans if succeeds.
 * When \a a is a singular matrix, the equation does not have
 * a unique answer and the null pointer is returned.
 */
__EXPORT zVec zLESolveRI(zMat a, zVec b, zVec ans);

/*! \brief linear equation solver by Gauss-Seidel method.
 *
 * zLESolveGS() solves linear equation \a a x = \a b based on
 * Gauss-Seidel's method.
 * The answer (x in the equation) is put into \a ans.
 * It is typically utilized for sparse equation, namely,
 * an equation with a sparse coefficient matrix \a a.
 * \return
 * zLESolveGS() returns a pointer \a ans if succeeds or the
 * iteration does not converge within ZM_MAX_ITER_NUM times
 * (defined in zm_misc.h).
 *
 * When the size mismatch occurs between \a a, \a b and \a ans,
 * or it fails to allocate working memory, the null pointer
 * is returned.
 * \notes
 * Since Gauss-Seidel's method is an iterating computation,
 * it might not converge to a certain vector.
 * Iteration is executed up to Z_MAX_ITER_NUM times.
 * If the iteration does not finish even after trying over
 * Z_MAX_ITER_NUM times, the function gives it up to calculate.
 */
__EXPORT zVec zLESolveGS(zMat a, zVec b, zVec ans);

__END_DECLS

#include <zm/zm_le_pivot.h>   /* pivoting */
#include <zm/zm_le_lu.h>      /* LU decomposition */
#include <zm/zm_le_lq.h>      /* LQ/QR decomposition */
#include <zm/zm_le_minv.h>    /* determinant and inverse matrix */
#include <zm/zm_le_mpinv.h>   /* Moore=Penrose inverse matrix */
#include <zm/zm_le_tridiag.h> /* tridiagonal equation */
#include <zm/zm_le_gen.h>     /* generalized linear equation */
#include <zm/zm_le_lyapnov.h> /* Lyapnov equation */

#endif /* __ZM_LE_H__ */
