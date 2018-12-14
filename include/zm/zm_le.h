/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le - linear equation.
 */

#ifndef __ZM_LE_H__
#define __ZM_LE_H__

#include <zm/zm_mat.h>

__BEGIN_DECLS

/* METHOD:
 * zBalancingColDST, zBalancingDST, zBalancing
 * - make a pair of matrix and vector balanced.
 *
 * 'zBalancingColDST()' destructively makes a matrix
 * column-balanced. Namely, each column of 'm' is divided
 * by the absolute-maximum component in the same column.
 * #
 * 'zBalancingDST()' destructively makes a pair of matrix
 * 'm' and vector 'v' balanced. Namely, each column of 'm'
 * is primarily divided by the absolute-maximum component
 * in the same column, and then each row of 'm' and the
 * corresponding component of 'v' by the absolute-maximum
 * value in the row.
 * #
 * Scalings are skipped when the balancing factor with
 * respect to the working column or row is zero.
 * #
 * 'zBalancing()' makes 'morg' and 'vorg' balanced, and
 * puts the result into 'm' and 'v'.
 * #
 * They work as preprocesses of a linear equation solver.
 * Actually, 'zLESolveGaussDST()' internally calls
 * 'zBalancingDST()'.
 * #
 * The column-balancing factors are stored in 's' for all
 * the functions if it is not the null pointer. It is used
 * for re-balancing. For the original equation 'morg x = vorg',
 * the balanced equation is 'm s x = v'. Thus, the solution
 * 'y' of 'm y = v' has to be amplified by 's' in such a
 * way as 'zVecAmpDRC( y, s )'.
 * [RETURN VALUES]
 * 'zBalancingColDST()' and 'zBalancingDST()' return no values.
 * #
 * 'zBalancing()' returns the false value if size-mismatch
 * happens among arguments. Otherwise, it returns the true
 * value.
 * [SEE ALSO]
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

/* METHOD:
 * zLESolveGaussDST, zLESolveGauss
 * - linear equation solver by Gauss elimination method.
 *
 * 'zLESolveGauss()' solves linear equation 'a x = b' based
 * on Gaussian elimination.
 * The answer is put into 'ans'.
 * 's' is used for column-balancing, if it is not the null
 * pointer.
 * #
 * 'zLESolveGaussDST()' destructively modifies 'a' and 'b'
 * while calculating.
 * [RETURN VALUE]
 * 'zLESolveGaussDST()' returns a pointer 'ans' if succeeds.
 * When 'a' is a singular matrix, the equation does not have
 * a unique answer and the null pointer is returned.
 */
__EXPORT zVec zLESolveGaussDST(zMat a, zVec b, zVec ans, zIndex index, zVec s);
__EXPORT zVec zLESolveGauss(zMat a, zVec b, zVec ans);

/* METHOD:
 * zLESolve_L, zLESolve_U, zLESolve_L_U, zLESolveLU
 * - linear equation solver by LU decomposition method.
 *
 * 'zLESolve_LE()' solves linear equation 'a x = b' based on
 * LU decomposition. The answer is put into 'ans'.
 * #
 * 'zLESolve_L_U()' solves the linear equation 'lmat umat x = b'.
 * 'lmat' and 'umat' are a lower and upper triangular matrix,
 * respectively.
 * 'index' is an index vector for order discription.
 * Each of 'lmat y = b' and 'umat x = y' is solved by
 * 'zLESolve_L()' and 'zLESolve_U()', respectively.
 * while the latter solves for 'x'.
 * #
 * Since 'zLUDecomp()' pivots the original coefficient
 * matrix, 'zLESolve_L()' needs the pivot index.
 * [RETURN VALUE]
 * 'zLESolveLU()' family functions return a pointer 'ans'
 * if succeeds.
 * When 'a' is a singular matrix, the equation does not have
 * a unique answer and these functions return the null pointer.
 */
__EXPORT zVec zLESolve_L(zMat lmat, zVec b, zVec ans, zIndex idx);
__EXPORT zVec zLESolve_U(zMat umat, zVec b, zVec ans);
__EXPORT zVec zLESolve_L_U(zMat lmat, zMat umat, zVec b, zVec ans, zIndex index);
__EXPORT zVec zLESolveLU(zMat a, zVec b, zVec ans);

/*\brief linear equation solver by residual iteration.
 *
 * 'zLESolveRI()' solves linear equation 'a x = b' based on
 * LU decomposition and improves the accuracy by residual
 * iteration. The answer is put into 'ans'.
 * \return \a ans if succeeds.
 * When \a a is a singular matrix, the equation does not have
 * a unique answer and the null pointer is returned.
 */
__EXPORT zVec zLESolveRI(zMat a, zVec b, zVec ans);

/* METHOD:
 * zLESolveGS - linear equation solver by Gauss-Seidel method.
 *
 * 'zLESolveGS()' solves linear equation 'a x = b'
 * based on Gauss-Seidel's method.
 * The answer ('x' in the equation) is put into 'ans'.
 * It is typically utilized for sparse equation, namely,
 * an equation with a sparse coefficient matrix 'a'.
 * [RETURN VALUE]
 * 'zLESolveGS()' returns a pointer 'ans' if succeeds
 * or the iteration does not converge within ZM_MAX_ITER_NUM
 * times (defined in 'zm_misc.h').
 * #
 * When the size mismatch occurs between 'a', 'b' and 'ans',
 * or it fails to allocate working memory, the null pointer
 * is returned.
 * [NOTES]
 * Since Gauss-Seidel's method is an iterating computation,
 * it might not converge to a certain vector.
 * Iteration is executed up to Z_MAX_ITER_NUM times.
 * If the iteration does not finish even after trying
 * over Z_MAX_ITER_NUM times, the function gives it up
 * to calculate.
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
