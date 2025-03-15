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
__ZM_EXPORT zMat zMatToHessenberg(const zMat m, zMat h, zMat p);

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
 * \sa
 * zMatEig
 */
__ZM_EXPORT bool zMatEigDQR(const zMat m, zCVec eigval, int iter);

/*! \brief the greatest and least eigenvalue of a matrix.
 *
 * zMatEigPower() calculates the greatest real eigenvalue of a square matrix \a m by the power method.
 * zMatEigPowerInv() calculates the least real eigenvalue of a square matrix \a m by the power method.
 * For the both functions, the corresponding eivenvector is put into \a eigvec.
 *
 * The methods are based on an iterative computation. \a iter is the maximum iteration number. If 0 is
 * assigned, it is replaced by Z_MAX_ITER_NUM (defined in zm_misc.h).
 * \return
 * zMatEigPower() and zMatEigPowerInv() return the computed eigenvalues.
 */
__ZM_EXPORT double zMatEigPower(const zMat a, zVec eigvec, int iter);
__ZM_EXPORT double zMatEigPowerInv(const zMat a, zVec eigvec, int iter);

/*! \brief eigenvalues and eigenvectors of a real square matrix.
 *
 * zMatEig() finds all eigenvalues and corresponding eigenvectors of a real square matrix \a m.
 * The eigenvalues are stored in a complex vector \a eigval, while the corresponding eigenvectors are
 * put in columns of a complex matrix \a eigbase.
 * Namely, the following equation holds:
 *  \a m \a eigbase = \a eigbase diag{ \a eigval }
 * It is based on the double QR method, and available for arbitrary real square matrices.
 * \return
 * zMatEig() returns the false value if
 *  - \a m is not a square matrix.
 *  - the sizes of \a m and \a eigbase are not equal.
 *  - the size of \a eigval does not coincide with the row size of \a m.
 *  - it fails to allocate internal working memory.
 * Otherwise, it returns the true value.
 * \sa
 * zMatEigDQR
 */
__ZM_EXPORT bool zMatEig(const zMat m, zCVec eigval, zCMat eigbase, int iter);

/*! \brief diagonalize a symmetric matrix.
 *
 * zMatSymEigBisec() and zMatSymEigJacobi() diagonalizes a symmetric matrix \a m. The former is based
 * on the bisection method proposed by J. W. Givens in 1954, while the latter on Jacobi's method with
 * Wilkinson's formula.
 * The result diagonal values are stored in a vector form \a eigval, while the transformation matrix
 * is stored in \a eigbase. Namely, the following equation holds:
 *   \a m \a eigbase = \a eigbase diag{ \a eigval }
 * \return
 * zMatSymEigBisec() and zMatSymEigJacobi() return the false value if
 *  - \a m is not a square matrix
 *  - the sizes of \a m and \a eigbase are not equal
 *  - the size of \a eigval does not coincide with the row size of \a m
 *  - it fails to allocate internal working memory
 * Otherwise, they return the true value.
 * \notes
 * When \a m is not symmetric, those functions do not work as expected.
 */
__ZM_EXPORT bool zMatSymEigBisec(const zMat m, zVec eigval, zMat eigbase);
__ZM_EXPORT bool zMatSymEigJacobi(const zMat m, zVec eigval, zMat eigbase);

/* singular value decomposition */

/*! \brief singular value decomposition (SVD).
 *
 * zMatSVD() decomposes a matrix \a m into the following form:
 *   \a m = \a u s \a v.
 * where \a u and \a v are orthonormal matrices and s is a diagonal matrix. Note that, different from
 * the standard SVD from, \a v is defined in non-transposed form.
 * The diagonal components of s are so-called singular values, which are stored in a vector \a sv.
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
__ZM_EXPORT int zMatSVD(const zMat m, zMat u, zVec sv, zMat v);

/*! \brief maximum singular value of a matrix.
 */
__ZM_EXPORT double zMatSingularValueMax(const zMat m);

/*! \brief minimum singular value of a matrix.
 */
__ZM_EXPORT double zMatSingularValueMin(const zMat m);

/*! \brief condition number of a matrix.
 */
__ZM_EXPORT double zMatCondNum(const zMat m);

__END_DECLS

#endif /* __ZM_MAT_EIG_H__ */
