/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_mat_eig - eigensystem of matrices.
 */

#ifndef __ZM_MAT_EIG_H__
#define __ZM_MAT_EIG_H__

#include <zm/zm_cmat.h>
#include <zm/zm_le.h>

__BEGIN_DECLS

/*! \brief transform a matrix to a Hessenberg matrix.
 *
 * zMatToHessenberg() transforms a matrix \a m to an upper Hessenberg matrix using Householder's method.
 * The result is put into \a h.
 * Particularly, \a h will be a tridiagonal matrix when \a m is symmetric.
 *
 * If \a h is the null pointer, \a m is directly transformed into an upper Hessenberg matrix.
 *
 * When \a p is not the null pointer, the transformation matrix is put into \a p, namely,
 *  \a p \a h \a p^T = \a m
 * \return
 * zMatToHessenberg() returns a pointer \a h.
 */
__ZM_EXPORT zMat zMatToHessenberg(zMat m, zMat h, zMat p);

/*! \brief eigensystem of a matrix by double QR method.
 *
 * zMatEigDQR() computes all complex eigenvalues of a matrix \a a. The eigenvalues are put into an array
 * \a eig. \a iter is the maximum number of iteration. When zero is given for \a iter, Z_MAX_ITER_NUM is
 * adopted instead.
 *
 * This code is originated from "HANDBOOK OF NUMERICAL METHODS" (H. Togawa, 1992, ISBN 4-7819-0868-3).
 * \return
 * zMatEigDQR() returns a boolean value. When it fails to allocate working memory, or \a a is a non-square
 * matrix, the false value is returned. Otherwise, the true value is returned.
 */
bool zMatEigDQR(zMat a, zComplex eig[], int iter);

/*! \brief the dominant eigenvalue of a matrix.
 *
 * zMatEigPower() calculates the dominant eigenvalue of matrix \a a by the power method. The corresponding
 * eivenvector is put into \a evec.
 *
 * The method might not converge to the true value. If the iteration does not be finished even if it
 * repeats more than Z_MAX_ITER_NUM (defined in zm_misc.h), the function gives up to calculate.
 * \return
 * zMatEigPower() returns the computed eigenvalue.
 */
__ZM_EXPORT double zMatEigPower(zMat a, zVec evec, int iter);
__ZM_EXPORT double zMatEigPowerInv(zMat a, zVec evec, int iter);

__ZM_EXPORT int zMatEig(zMat m, zComplex eig[], zCVec eigv[], int iter);

/*! \brief diagonalize a symmetric matrix by bisection method.
 *
 * zMatSymEigBisec() diagonalizes a symmetric matrix \a m based on the bisection method proposed by
 * J. W. Givens in 1954. The result diagonal values are stored in a vector form \a eig, while the
 * transformation matrix is stored in \a r.
 * \return
 * zMatSymEigBisec() returns a pointer \a eig.
 * \notes
 * When \a m is not symmetric, zMatSymEigBisec() does not work as expected.
 */
__ZM_EXPORT zVec zMatSymEigBisec(zMat m, zVec eig, zMat r);

/*! \brief diagonalize a symmetric matrix by Jacobi method.
 *
 * zMatSymEigJacobi() diagonalizes a symmetric matrix \a m by Jacobi method. The result diagonal values
 * are stored in a vector form \a eig, while the transformation matrix is stored in \a r.
 * The iteration will finish when the whole non-diagonal components become smaller than zTOL, or the
 * number of steps exceeds Z_MAX_ITER_NUM (defined in zm_misc.h).
 *
 * This implimentation of Jacobi rotation transformation is based on Wilkinson's formula.
 * \return
 * zMatSymEigJacobi() returns a pointer \a eig.
 * \notes
 * When \a m is not symmetric, zMatSymEigJacobi() does not work as expected.
 */
__ZM_EXPORT zVec zMatSymEigJacobi(zMat m, zVec eig, zMat r);

/* singular value decomposition */

/*! \brief singular value decomposition.
 *
 * zMatSVD() decomposes a matrix \a m into the following form:
 *   \a m = \a u \a s \a v.
 * where \a u and \a v are orthogonal matrices and \a s is a diagonal matrix. Note that, different from
 * a standard SVD from, \a v is defined in non-transposed form.
 * The diagonal components of \a s are so-called singular values, which are stored in a vector \a sv.
 * The number of non-zero singular values coincides with the rank of \a m.
 *
 * The original sizes of each matrix have to be
 *  \a sv: n x 1
 *  \a u: n x n
 *  \a v: n x m,
 * respectively. However, the row size of \a v possibly changes after the function call, due to the rank
 * regression.
 * \return
 * zMatSVD() returns the rank of \a m, namely, the number of non-zero singular values of \a m.
 */
__ZM_EXPORT int zMatSVD(zMat m, zVec sv, zMat u, zMat v);

/*! \brief maximum singular value of a matrix.
 */
__ZM_EXPORT double zMatSingularValueMax(zMat m);

/*! \brief minimum singular value of a matrix.
 */
__ZM_EXPORT double zMatSingularValueMin(zMat m);

/*! \brief condition number of a matrix.
 */
__ZM_EXPORT double zMatCondNum(zMat m);

__END_DECLS

#endif /* __ZM_MAT_EIG_H__ */
