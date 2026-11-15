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
 * zMatToHessenberg() transforms a matrix \a m to an upper Hessenberg matrix using Householder's method:
 *  A. S. Householder, Unitary Triangularization of a Nonsymmetric Matrix, Journal of the ACM,
 *  5(4):339-342, 1958.
 * The result is put into \a h.
 * Particularly, \a h will be a tridiagonal matrix when \a m is symmetric.
 *
 * If \a h is the null pointer, \a m is directly transformed into an upper Hessenberg matrix.
 *
 * When \a p is not the null pointer, the transformation matrix is put into \a p, namely,
 *  \a p \a h \a p^T = \a m
 * Each Householder reflection is skipped only if the column to be reduced is exactly zero, so that
 * the result does not depend on the scale of \a m.
 * \return
 * zMatToHessenberg() returns a pointer \a h, or \a m if \a h is the null pointer. If \a m is not
 * a square matrix, the sizes of matrices mismatch, or it fails to allocate internal working memory,
 * the null pointer is returned.
 */
__ZM_EXPORT zMat zMatToHessenberg(const zMat m, zMat h, zMat p);

/*! \brief eigenvalues and eigenvectors of a real square matrix.
 *
 * zMatEigDQR() finds all eigenvalues and corresponding eigenvectors of a real square matrix \a m
 * based on the implicit double-shift QR method by J. G. F. Francis:
 *  J. G. F. Francis, The QR Transformation: A Unitary Analogue to the LR Transformation -- Part 1,
 *  The Computer Journal, 4(3):265-271, 1961.
 *  J. G. F. Francis, The QR Transformation -- Part 2, The Computer Journal, 4(4):332-345, 1962.
 * The implementation combines the following techniques.
 * A subdiagonal component is regarded as zero relatively to its neighboring diagonal components,
 * the first column of the double-shift polynomial is computed with the differences of the diagonal
 * components from the shifts to avoid cancellation, and the exceptional shifts relative to the last
 * diagonal component are applied at the 10th and 20th iterations, all of which follow a procedure by
 * R. S. Martin, G. Peters and J. H. Wilkinson:
 *  R. S. Martin, G. Peters and J. H. Wilkinson, The QR Algorithm for Real Hessenberg Matrices,
 *  Numerische Mathematik, 14:219-231, 1970.
 * Eigenvectors by back substitution on the real Schur form by G. Peters and J. H. Wilkinson:
 *  G. Peters and J. H. Wilkinson, Eigenvectors of Real and Complex Matrices by LR and QR
 *  Triangularizations, Numerische Mathematik, 16:181-204, 1970.
 * A pivot smaller than the norm of the matrix multiplied by a small epsilon value is replaced with
 * that value, which extends the replacement of a zero pivot in another procedure in the same paper.
 *
 * The eigenvalues are stored in a complex vector \a eigval, while the corresponding eigenvectors are
 * put in columns of a complex matrix \a eigbase. Namely, the following equation holds:
 *  \a m \a eigbase = \a eigbase diag{ \a eigval }
 * If \a eigbase is the null pointer, only the eigenvalues are computed, which saves the computation
 * of the Schur vectors.
 * The eigenvalues are stored in the order of the diagonal blocks of the real Schur form of \a m,
 * and a pair of complex conjugate eigenvalues are stored in adjacent components, where the one with
 * the positive imaginary part comes first. Each eigenvector is normalized so that its norm is one
 * and its absolutely maximum component is a positive real number.
 * The eigenvectors are computed by back substitution on the real Schur form, so that they are not
 * necessarily linearly independent for a defective matrix.
 *
 * \a iter is the maximum number of iterations to find each eigenvalue (or each pair of complex
 * conjugate eigenvalues). When zero is given for \a iter, Z_MAX_ITER_NUM is adopted instead.
 *
 * zMatEig() is defined as an alias of zMatEigDQR().
 * \return
 * zMatEigDQR() and zMatEig() return the false value if
 *  - \a m is not a square matrix.
 *  - the size of \a eigval does not coincide with the row size of \a m.
 *  - \a eigbase is not the null pointer and its size does not coincide with that of \a m.
 *  - they fail to allocate internal working memory.
 *  - the iteration does not converge within \a iter times for an eigenvalue.
 * Otherwise, they return the true value.
 */
__ZM_EXPORT bool zMatEigDQR(const zMat m, zCVec eigval, zCMat eigbase, int iter);
#define zMatEig(mat,eigval,eigbase,iter) zMatEigDQR( mat, eigval, eigbase, iter )

/*! \brief the dominant and least-dominant eigenvalues of a matrix.
 *
 * zMatEigPower() calculates the dominant eigenvalue, namely, the real eigenvalue with the largest absolute
 * value, of a square matrix \a m by the power method.
 * zMatEigPowerInv() calculates the real eigenvalue with the smallest absolute value of a square matrix
 * \a m by the inverse power method, where \a m is LU-decomposed in advance, and the linear equation
 * with the matrix is solved at each step instead of computing the inverse of \a m.
 * For the both functions, the corresponding eigenvector is put into \a eigvec as a unit vector.
 *
 * The methods are based on an iterative computation from a random initial vector, which is generated by
 * a local random number generator with a fixed seed. Hence, the result is reproducible.
 * \a iter is the maximum iteration number. If 0 is assigned, it is replaced by Z_MAX_ITER_NUM (defined
 * in zm_misc.h).
 * The iteration is regarded as converged if the residual |a v - l v| for the eigenvalue l and the unit
 * eigenvector v (for zMatEigPowerInv(), those of the inverse of \a m) is not greater than zTOL times
 * |l|, so that the result does not depend on the scale of \a m. If it does not converge within \a iter
 * steps, a warning message is output and the last estimates are returned. Note that the power method
 * does not converge if the dominant eigenvalue is not unique in the absolute value, e.g., if a pair of
 * eigenvalues with opposite signs or of complex conjugates are dominant. The convergence is slow if the
 * absolute values of the dominant and the second dominant eigenvalues are close.
 * \return
 * zMatEigPower() and zMatEigPowerInv() return the computed eigenvalues. If \a m is not square or the
 * size of \a eigvec is not consistent with \a m, they return 0. zMatEigPowerInv() also returns 0 if \a m
 * is singular, and \a eigvec is set for the zero vector.
 */
__ZM_EXPORT double zMatEigPower(const zMat a, zVec eigvec, int iter);
__ZM_EXPORT double zMatEigPowerInv(const zMat a, zVec eigvec, int iter);

/*! \brief diagonalize a symmetric matrix.
 *
 * zMatSymEigBisec() diagonalizes a symmetric matrix \a m based on the bisection method proposed by
 * J. W. Givens:
 *  W. Givens, Numerical Computation of the Characteristic Values of a Real Symmetric Matrix,
 *  Oak Ridge National Laboratory, ORNL-1574, 1954.
 * The eigenvalues are stored in descending order, where multiple eigenvalues are repeated. The
 * eigenvectors are computed by the inverse iteration, in which those for a cluster of eigenvalues, the
 * differences of which are not greater than 1e-3 times the norm of \a m, are orthogonalized with each
 * other, following a procedure by G. Peters and J. H. Wilkinson:
 *  G. Peters and J. H. Wilkinson, The Calculation of Specified Eigenvectors by Inverse Iteration,
 *  in J. H. Wilkinson and C. Reinsch (eds.), Handbook for Automatic Computation, Vol. II, Linear
 *  Algebra, Contribution II/18, pp. 418-439, Springer, 1971.
 * \a m is scaled by a power of two in advance to avoid overflow and underflow.
 *
 * zMatSymEigJacobi() also diagonalizes a symmetric matrix \a m based on the cyclic Jacobi's method
 * with an update formula by Rutishauser:
 *  C. G. J. Jacobi, Ueber ein leichtes Verfahren, die in der Theorie der Saecularstoerungen
 *  vorkommenden Gleichungen numerisch aufzuloesen, Journal fuer die reine und angewandte Mathematik,
 *  30:51-94, 1846.
 *  H. Rutishauser, The Jacobi Method for Real Symmetric Matrices, Numerische Mathematik, 9:1-10, 1966.
 * It regards an off-diagonal component m_ij as zero if |m_ij| <= eps sqrt(|m_ii| |m_jj|), where eps
 * is the machine epsilon, so that the result does not depend on the scale of \a m.
 * The threshold of convergence is determined based on a criterion by J. Demmel and K. Veselic:
 *  J. Demmel and K. Veselic, Jacobi's Method is More Accurate than QR, SIAM Journal on Matrix
 *  Analysis and Applications, 13(4):1204-1245, 1992.
 *
 * For the both functions, the result diagonal values are stored in a vector form \a eigval, while the
 * transformation matrix is stored in \a eigbase. Namely, the following equation holds:
 *   \a m \a eigbase = \a eigbase diag{ \a eigval }
 * The diagonal values computed by zMatSymEigJacobi() are not sorted.
 * \return
 * zMatSymEigBisec() and zMatSymEigJacobi() return the false value if
 *  - \a m is not a square matrix
 *  - the sizes of \a m and \a eigbase are not equal
 *  - the size of \a eigval does not coincide with the row size of \a m
 *  - it fails to allocate internal working memory
 * zMatSymEigJacobi() also returns the false value if the off-diagonal components do not converge to
 * zero within Z_MAX_ITER_NUM sweeps, where \a eigval and \a eigbase hold the incomplete results.
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
 * where s is a diagonal matrix, the columns of \a u and the rows of \a v are orthonormal. Note that,
 * different from the standard SVD form, \a v is defined in non-transposed form.
 * The diagonal components of s are so-called singular values, which are stored in a vector \a sv in
 * descending order. A singular value not greater than zTOL times the greatest one is regarded as zero.
 * The number of non-zero singular values is the rank of \a m.
 *
 * The original sizes of each matrix have to be
 *  \a sv: n x 1
 *  \a u: n x n
 *  \a v: n x m,
 * respectively, where \a m is an n x m matrix. The row size of \a v is changed to the rank of \a m.
 * Only the first rank columns of \a u are meaningful, and the rest of the components of \a sv are zeros.
 *
 * \a m is first reduced to a triangular matrix by Householder reflections, the rows of which are
 * orthogonalized by the one-sided Jacobi's method. It applies LQ (or QR) decomposition by Householder
 * reflections as a preprocessing by Z. Drmac and K. Veselic, and then, reduces the size of the matrix
 * to be orthogonalized with the one-sided Jacobi's method by M. R. Hestenes:
 *  Z. Drmac and K. Veselic, New Fast and Accurate Jacobi SVD Algorithm I, SIAM Journal on Matrix
 *  Analysis and Applications, 29(4):1322-1342, 2008.
 *  Z. Drmac and K. Veselic, New Fast and Accurate Jacobi SVD Algorithm II, SIAM Journal on Matrix
 *  Analysis and Applications, 29(4):1343-1362, 2008.
 *  M. R. Hestenes, Inversion of Matrices by Biorthogonalization and Related Results, Journal of the
 *  Society for Industrial and Applied Mathematics, 6(1):51-90, 1958.
 * \return
 * zMatSVD() returns the rank of \a m. If the sizes of the matrices and the vector do not match, it
 * fails to allocate internal working memory, or the one-sided Jacobi's method does not converge, it
 * returns -1.
 */
__ZM_EXPORT int zMatSVD(const zMat m, zMat u, zVec sv, zMat v);

/*! \brief maximum and minimum singular values and condition number of a matrix.
 *
 * zMatSingularValueMax() and zMatSingularValueMin() calculate the maximum and minimum singular values of
 * a matrix \a m, respectively. zMatCondNum() calculates the condition number of \a m, namely, the ratio
 * of the maximum singular value to the minimum one.
 * They are computed by the singular value decomposition zMatSVD(). Hence, a singular value not greater
 * than zTOL times the maximum one is regarded as zero.
 * \return
 * zMatSingularValueMax() and zMatSingularValueMin() return the maximum and minimum singular values of
 * \a m, respectively. zMatSingularValueMin() returns 0 if \a m is not full-rank.
 * zMatCondNum() returns the condition number of \a m, or HUGE_VAL if \a m is not full-rank.
 * If they fail to allocate internal working memory or zMatSVD() fails, they return NaN.
 */
__ZM_EXPORT double zMatSingularValueMax(const zMat m);
__ZM_EXPORT double zMatSingularValueMin(const zMat m);
__ZM_EXPORT double zMatCondNum(const zMat m);

__END_DECLS

#endif /* __ZM_MAT_EIG_H__ */
