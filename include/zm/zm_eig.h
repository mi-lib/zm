/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_eig - eigenvalue analysis.
 */

#ifndef __ZM_EIG_H__
#define __ZM_EIG_H__

#include <zm/zm_cmat.h>
#include <zm/zm_le.h>

__BEGIN_DECLS

/* Householder conversion to Hessian matrices */

/*! \brief Householder transformation.
 *
 * zHouseholder() directly transforms a given matrix \a m
 * by Householder's method with a projection vector \a u.
 * Only from \a from th to \a to th components of \a u are
 * valid for the operation. \a v and \a w are for working
 * spaces with the same size with \a u.
 * When \a p is not the null pointer, transformation matrix
 * is put into \a p, namely,
 *  \a p \a h \a p^T = \a m
 *
 * zHouseholderVec() computes a projection vector \a u
 * which reflects the \a col th column of \a m to an
 * only-one-component column with zHouseholder().
 * \a from and \a to specify a valid range of the reflection.
 * \return
 * zHouseholder() returns no value.
 * zHouseholderVec() returns a pointer \a u.
 */
__ZM_EXPORT void zHouseholder(zMat m, zMat p, int from, int to, zVec u, zVec v, zVec w);
__ZM_EXPORT zVec zHouseholderVec(zMat m, int col, int from, int to, zVec u);

/*! \brief Hessenberg matrix using Householder transformation.
 *
 * zHess() transforms a given matrix \a m to an upper
 * Hessenberg matrix using Householder's method.
 * The result is put into \a h.
 * Particularly, \a h will be a tridiagonal matrix when \a m
 * is symmetric.
 *
 * If \a h is the null pointer, \a m is directly transformed
 * into an upper Hessenberg matrix.
 *
 * When \a p is not the null pointer, transformation matrix
 * is put into \a p, namely,
 *  \a p \a h \a p^T = \a m
 * \return
 * zHess() returns a pointer \a h.
 */
__ZM_EXPORT zMat zHess(zMat m, zMat h, zMat p);

/* eigenvalue analysis by double QR method and inverse iteration */

/*! \brief double QR method.
 *
 * zEigDQR() computes all complex eigenvalues of a matrix
 * \a a. The resultant values are put into an array \a eig.
 * \a iter is the maximum number of iteration. When
 * zero is given for \a iter, Z_MAX_ITER_NUM is adopted
 * instead.
 *
 * This code is originated from "HANDBOOK OF NUMERICAL
 * METHODS" (H. Togawa, 1992, ISBN 4-7819-0868-3)
 * \return
 * zEigDQR() returns a boolean. When it fails to
 * allocate the working memory, or \a a is a non-square
 * matrix, the false value is returned. Otherwise, the
 * true value is returned.
 */
bool zEigDQR(zMat a, zComplex eig[], int iter);

/*! \brief calculate the dominant eigenvalue.
 *
 * zEigPower() calculates the dominant eigenvalue of matrix
 * \a a by the power method. The corresponding eivenvector
 * is put into \a evec.
 *
 * Since the power method is an iterative computation, it
 * might not converge to the true value. If the iteration
 * does not be finished even if it repeats more than
 * Z_MAX_ITER_NUM(defined in zm_misc.h), the function gives
 * it up to calculate.
 * \return
 * zEigPower() returns the dominant eigenvalue calculated.
 */
__ZM_EXPORT double zEigPower(zMat a, zVec evec, int iter);
__ZM_EXPORT double zEigPowerInv(zMat a, zVec evec, int iter);

__ZM_EXPORT int zEigSystem(zMat m, zComplex eig[], zCVec eigv[], int iter);

/* diagonalization by bisection method */

/*! \brief compute eigenvalues and eigenvectors based on bisection method
 * (by J. W. Givens 1954).
 */
__ZM_EXPORT zVec zEigSymBisec(zMat m, zVec eig, zMat r);

/* diagonalization by Jacobi's method */

/*! \brief compute eigenvalues of a symmetric matrix by Jacobi's method.
 *
 * zEigSymJacobi() diagonalizes a symmetric matrix \a m by
 * Jacobi's method. The iteration will finish when the whole
 * non-diagonal components become tiny.
 * In practice, the iteration is done up to Z_MAX_ITER_NUM
 * (defined in zm_misc.h).
 * If the iteration does not finish even after Z_MAX_ITER_NUM
 * times, the function gives up the calculation.
 *
 * This implimentation of Jacobi's rotation transformation
 * is based on Wilkinson's formula.
 * \return
 * zEigSymJacobi() returns a pointer \a eig.
 * \notes
 * When \a m is not symmetric, zEigSymJacobi() does not work
 * as expected.
 */
__ZM_EXPORT zVec zEigSymJacobi(zMat m, zVec eig, zMat r);

/* singular value decomposition */

/*! \brief singular value decomposition.
 *
 * zSVD() does the singular value decomposition.
 * It decomposes matrix \a m into the following form:
 *   M = U S V
 * where U and V are orthogonal matrices and S is a diagonal
 * matrix. Note that, different from a standard SVD from, V
 * is defined in non-transposed form.
 * The diagonal components of S are so-called singular values,
 * which are stored in a vector \a sv.
 * The number of non-zero singular values coincides with the
 * rank of \a m.
 *
 * The original sizes of each matrix have to be
 *  \a sv: n x 1,
 *  \a u: n x n
 *  \a v: n x m,
 * respectively. However, the row size of \a v possibly changes
 * after the function call, due to the rank regression.
 * \return
 * zSVD() returns the rank of \a m, namely, the number of
 * non-zero singular values of \a m.
 */
__ZM_EXPORT int zSVD(zMat m, zVec sv, zMat u, zMat v);

/*! \brief maximum singular value of a matrix.
 */
__ZM_EXPORT double zSVMax(zMat m);

/*! \brief minimum singular value of a matrix.
 */
__ZM_EXPORT double zSVMin(zMat m);

/*! \brief condition number of a matrix.
 */
__ZM_EXPORT double zMatCondNum(zMat m);

__END_DECLS

#endif /* __ZM_EIG_H__ */
