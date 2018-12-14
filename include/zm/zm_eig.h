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

/* METHOD:
 * zHouseholder, zHouseholderVec
 * - Householder transformation.
 *
 * 'zHouseholder()' directly transforms a given
 * matrix 'm' by Householder's transformation, using
 * a projection vector 'u'.
 * Only from 'from'th to 'to'th components of 'u'
 * are valid for the operation. 'v' and 'w' are
 * for working spaces with the same size with 'u'.
 * When 'p' is not the null pointer, transformation matrix
 * is put where pointed by 'p', namely,
 *  'p' 'h' 'p'^T = 'm'
 *
 * 'zHouseholderVec()' computes a projection vector
 * 'u' which reflects the 'col'th column of 'm' to
 * an only-one-component column with 'zHouseholder()'.
 * 'from' and 'to' specify a valid range of the reflection.
 * [RETURN VALUE]
 * 'zHouseholder()' returns no value.
 * 'zHouseholderVec()' returns a pointer 'u'.
 */
__EXPORT void zHouseholder(zMat m, zMat p, int from, int to, zVec u, zVec v, zVec w);
__EXPORT zVec zHouseholderVec(zMat m, int col, int from, int to, zVec u);

/* METHOD:
 * zHess
 * - Hessenberg matrix using Householder transformation.
 *
 * 'zHess()' transforms a given matrix 'm' to an upper
 * Hessenberg matrix using Householder's transformation.
 * The result will be put into 'h'.
 * Particularly, 'h' will be a tridiagonal matrix when
 * 'm' is symmetric.
 *
 * If 'h' is the null pointer, 'm' is directly transformed
 * into an upper Hessenberg matrix.
 *
 * When 'p' is not the null pointer, transformation matrix
 * is put where pointed by 'p', namely,
 *  'p' 'h' 'p'^T = 'm'
 * [RETURN VALUE]
 * 'zHess()' returns a pointer 'h'.
 */
__EXPORT zMat zHess(zMat m, zMat h, zMat p);

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

/* METHOD:
 * zEigPower - calculate the dominant eigenvalue.
 *
 * 'zEigPower()' calculates the dominant eigenvalue of matrix 'a'
 * by power method. The corresponding eivenvector is stored where
 * is pointed by 'evec'.
 *
 * Since power method is an iterative computation, it might not
 * converge to the true value. If the iteration does not be finished
 * even if it repeats more than Z_MAX_ITER_NUM(defined in
 * 'zm_misc.h'), the function gives it up to calculate.
 * [RETURN VALUE]
 * 'zEigPower()' returns the dominant eigenvalue calculated.
 */
__EXPORT double zEigPower(zMat a, zVec evec, int iter);
__EXPORT double zEigPowerInv(zMat a, zVec evec, int iter);

__EXPORT int zEigSystem(zMat m, zComplex eig[], zCVec eigv[], int iter);

/* diagonalization by bisection method */

/* zEigSymBisec
 * - compute eigenvalues and eigenvectors based on bisection method
 *   (by J. W. Givens 1954).
 */
__EXPORT zVec zEigSymBisec(zMat m, zVec eig, zMat r);

/* diagonalization by Jacobi's method */

/* METHOD:
 * zEigSymJacobi
 * - compute eigenvalues of a symmetric matrix by Jacobi's method.
 *
 * 'zEigSymJacobi()' diagonalizes a symmetric matrix 'm'
 * by Jacobi's method. The iteration will be finished
 * when the whole non-diagonal components become tiny.
 * In practice, the iteration is done up to Z_MAX_ITER_NUM
 * (defined in 'zm_misc.h').
 * If the iteration does not be finished even after
 * Z_MAX_ITER_NUM times, the function gives up
 * the calculation.
 *
 * This implimentation of Jacobi's rotation transformation
 * is based on Wilkinson's formula.
 * [RETURN VALUE]
 * 'zEigSymJacobi()' returns a pointer 'eig'.
 * [NOTES]
 * When 'm' is not symmetric, 'zEigSymJacobi()' does not
 * work validly.
 */
__EXPORT zVec zEigSymJacobi(zMat m, zVec eig, zMat r);

/* singular value decomposition */

/* METHOD:
 * zSVD
 * - singular value decomposition.
 *
 * 'zSVD()' is an implementation of singular value decomposition.
 * It decomposes matrix 'm' into the following form:
 *   M = U S V
 * where U and V are orthogonal matrices and S is a diagonal matrix.
 * Note that, different from a standard SVD from, V matrix is
 * defined in non-transpose form.
 * The diagonal components of S are so-called singular values,
 * which will be stored in a vector 'sv'.
 * The number of non-zero singular values coincides with the
 * rank of 'm'.
 *
 * The original sizes of each matrix have to be
 *  'sv': n x 1,
 *  'u': n x n
 *  'v': n x m,
 * respectively. However, the row size of 'v' possibly changes
 * after the function call, due to the rank regression.
 * [RETURN VALUE]
 * 'zSVD()' returns the rank of 'm', namely, the number of
 * non-zero singular values of 'm'.
 */
__EXPORT int zSVD(zMat m, zVec sv, zMat u, zMat v);

/*! \brief maximum singular value of a matrix.
 */
__EXPORT double zSVMax(zMat m);

/*! \brief minimum singular value of a matrix.
 */
__EXPORT double zSVMin(zMat m);

/*! \brief condition number of a matrix.
 */
__EXPORT double zMatCondNum(zMat m);

__END_DECLS

#endif /* __ZM_EIG_H__ */
