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
 * zMatColBalancingDST() destructively makes a matrix column-balanced. Namely, each column of \a m is
 * divided by the absolute-maximum component in the same column.
 *
 * zMatVecRowBalancingDST() destructively makes a pair of matrix \a m and vector \a v row-balanced.
 * Namely, each column of \a m is primarily divided by the absolute-maximum component in the same column,
 * and then, each row of \a m and the corresponding component of \a v by the absolute-maximum value in
 * the row.
 *
 * zMatVecBalancingDST() destructively makes a pair of matrix \a m and vector \a v both column- and row-
 * balanced.
 *
 * zMatVecBalancing() makes \a morg and \a vorg balanced, and puts the result into \a m and \a v.
 *
 * They work as preprocesses of a linear equation solver. Actually, zLESolveGaussDST() internally calls
 * zMatVecBalancingDST().
 *
 * The column-balancing factors are stored in \a s for all the functions if it is not the null pointer.
 * It is used for re-balancing. For the original equation \a morg x = \a vorg, the balanced equation is
 * \a m \a s x = v. Thus, the solution \a y of \a m y = v has to be amplified by \a s in such a way as
 * zVecAmpDRC( y, s ).
 * \return
 * zMatColBalancingDST(), zMatVecRowBalancingDST(), and zMatVecBalancingDST() return no values.
 *
 * zMatVecBalancing() returns the false value if size-mismatch happens among arguments. Otherwise, it
 * returns the true value.
 * \sa
 * zLESolveGaussDST
 */
__ZM_EXPORT void zMatColBalancingDST(zMat m, zVec s);
__ZM_EXPORT void zMatVecRowBalancingDST(zMat m, zVec v);
__ZM_EXPORT void zMatVecBalancingDST(zMat m, zVec v, zVec s);
__ZM_EXPORT bool zMatVecBalancing(const zMat morg, const zVec vorg, zMat m, zVec v, zVec s);

/*! \brief make a pair of matrices balanced.
 *
 * zMatMatBalancingDST() destructively makes a pair of matrices \a m1 and \a m2 both column- and row-
 * balanced.
 *
 * The column-balancing factors are stored in \a s to be used for re-balancing, if it is not the null
 * pointer.
 * \return
 * zMatMatBalancingDST() does not return any value.
 * \sa
 * zMatVecBalancingDST
 */
__ZM_EXPORT void zMatMatBalancingDST(zMat m1, zMat m2, zVec s);

/*! \brief directly make a pair of matrices and a vector balanced.
 *
 * zMatMatVecRowBalancingDST() destructively makes a pair of matrices \a m1 and \a m2 and a vector \a v
 * row-balanced simultaneously.
 *
 * zMatMatVecBalancingDST() destructively makes a pair of matrices \a m1 and \a m2 and a vector \a v
 * both column- and row- balanced simultaneously.
 * The column-balancing factors are stored in \a s to be used for re-balancing, if it is not the null
 * pointer.
 * \return
 * zMatMatVecRowBalancingDST() and zMatMatVecBalancingDST() do not return any value.
 * \sa
 * zMatVecBalancingDST, zMatMatBalancingDST
 */
__ZM_EXPORT void zMatMatVecRowBalancingDST(zMat m1, zMat m2, zVec v);
__ZM_EXPORT void zMatMatVecBalancingDST(zMat m1, zMat m2, zVec v, zVec s);

/*! \brief residual b - a x.
 *
 * \notes Sizes are not checked. It is not recommended to use this function in user programs.
 */
__ZM_EXPORT zVec zLEResidual(const zMat a, const zVec b, const zVec x, zVec res);

/*! \brief linear equation solver by Gauss elimination method.
 *
 * zLESolveGauss() solves linear equation \a a x = \a b based on Gaussian elimination.
 * The answer is put into \a ans.
 * \a s is used for column-balancing, if it is not the null pointer.
 *
 * zLESolveGaussDST() destructively modifies \a a and \a b while calculating.
 * \return
 * zLESolveGaussDST() returns a pointer \a ans if succeeds. When \a a is a singular matrix, the equation
 * does not have a unique answer and the null pointer is returned.
 */
__ZM_EXPORT zVec zLESolveGaussDST(zMat a, zVec b, zVec ans, zIndex index, zVec s);
__ZM_EXPORT zVec zLESolveGauss(const zMat a, const zVec b, zVec ans);

/*! \brief linear equation solver by LU decomposition method.
 *
 * zLESolveLU() solves linear equation \a l \a u x = \a b based on LU decomposition, where \a l and \a u
 * are supposed to be a lower and an upper triangular matrices, respectively. The answer is put into
 * \a ans.
 * \a index is an index vector for order discription.
 * \a l y = \a b and \a u x = y are solved by zLESolveL() and zLESolveU(), respectively.
 *
 * zLESolveL() needs the pivot index since zLUDecomp() pivots the original coefficient matrix.
 * \return
 * zLESolveLU(), zLESolveL() and zLESolveU() return a pointer \a ans, if succeeding. Otherwise, they
 * return the null pointer if the original matrix is singular.
 * \sa
 * zLUDecomp
 */
__ZM_EXPORT zVec zLESolveL(const zMat l, const zVec b, zVec ans, const zIndex idx);
__ZM_EXPORT zVec zLESolveU(const zMat u, const zVec b, zVec ans);
__ZM_EXPORT zVec zLESolveLU(const zMat l, const zMat u, const zVec b, zVec ans, const zIndex index);

/*! \brief linear equation solver by LU decomposition method and residual iteration.
 *
 * zLESolveLURI() solves linear equation \a a x = \a b based on LU decomposition and improves accuracy
 * by residual iteration. The answer is put into \a ans.
 * \return
 * zLESolveLURI() returns a pointer \a ans if succeeds.
 * If \a a is not a regular matrix, meaning that the equation does not have a unique answer, it returns
 * the null pointer.
 */
__ZM_EXPORT zVec zLESolveLURI(const zMat a, const zVec b, zVec ans);

/*! \brief linear equation solver by Gauss-Seidel method.
 *
 * zLESolveGaussSeidel() solves linear equation \a a x = \a b based on Gauss-Seidel's method.
 * The answer (x in the equation) is put into \a ans.
 * It is typically utilized for sparse equation, namely, an equation with a sparse coefficient matrix \a a.
 * \return
 * zLESolveGaussSeidel() returns a pointer \a ans if succeeds or the iteration does not converge within
 * ZM_MAX_ITER_NUM times (defined in zm_misc.h).
 *
 * If the sizes of \a a, \b, and \a ans do not match, or it fails to allocate working memory, it returns
 * the null pointer.
 * \notes
 * Since Gauss-Seidel's method is an iterating computation, it might not converge to a certain vector.
 * Iteration is executed up to Z_MAX_ITER_NUM times. If the iteration does not finish even after trying
 * over Z_MAX_ITER_NUM times, the function gives it up to calculate.
 */
__ZM_EXPORT zVec zLESolveGaussSeidel(const zMat a, const zVec b, zVec ans);

__END_DECLS

#include <zm/zm_le_pivot.h>     /* pivoting */
#include <zm/zm_le_lu.h>        /* LU decomposition */
#include <zm/zm_le_lq.h>        /* LQ/QR decomposition */
#include <zm/zm_le_mat_inv.h>   /* determinant and inverse matrix */
#include <zm/zm_le_mat_mpinv.h> /* Moore-Penrose inverse matrix */
#include <zm/zm_le_tridiag.h>   /* tridiagonal equation */
#include <zm/zm_le_gen.h>       /* generalized linear equation */
#include <zm/zm_le_lyapnov.h>   /* Lyapnov equation */

#endif /* __ZM_LE_H__ */
