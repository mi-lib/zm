/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_pex_eq - polynomial expression class: polynomial equation solver.
 */

#ifndef __ZM_PEX_EQ_H__
#define __ZM_PEX_EQ_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief numerical solution of polynomial equation.
 *
 * zPexBHDST(), zPexBH() and zPexDKA() solve a polynomial equation \a pex=0 based on Bairstow-Hitchcock's
 * method and Durand-Kerner-Aberth's method, respectively.
 * All answers are put into an array pointed by \a ans as complex numbers. \a ans must have enough size.
 * \a tol is a tolerance to terminate iteration relative to the scale of the problem. If zero or a negative
 * value is given, zTOL is used instead.
 * \a iter is the maximum number of iteration. When it runs over \a iter times, it is aborted. If zero is
 * given for \a iter, Z_MAX_ITER_NUM defined in zm_misc.h is used instead.
 *
 * zPexDKA() regards the iteration as converged if the residual |\a pex(z)| at every approximate root z is
 * not greater than \a tol times the scale of the rounding error of the value, sum_k |a_k| |z|^k, where
 * a_k is the coefficient of k-dimension term, after which the correction continues while the residuals
 * decrease at most ten times to improve the accuracy of ill-conditioned roots such as close roots.
 *
 * zPexBH() and zPexBHDST() regard the iteration for a quadratic factor x^2 + p x + q as converged if the
 * updates of p and q are not greater than \a tol times the scale of the roots of the factor
 * m = |p| + sqrt(|q|) and its square, respectively. The iteration for each factor starts from several
 * factors with a pair of roots on the unit circle of the scaled variable until it converges.
 *
 * zPexBHDST() is a destructive version of zPexBH(), in which \a pex is destroyed in the computation process.
 * Different from zPexBH(), it does neither refine nor polish the answers, since it does not keep the original
 * polynomial.
 * \note
 * Since Bairstow-Hitchcock's method is an iteration with deflations, zPexBH() may not converge or lose
 * accuracy for a polynomial of a high degree (about fifteen or higher), even though it refines and polishes
 * the answers by further iterations. It may also result in a wrong root for a polynomial of which the
 * magnitudes of the roots are widely spread. zPexDKA() is advantageous to it from these viewpoints.
 * \return
 * zPexBHDST(), zPexBH() and zPexDKA() return a pointer \a ans if they succeed. Otherwise, they return the
 * null pointer.
 */
__ZM_EXPORT zCVec zPexBHDST(zPex pex, zCVec ans, double tol, int iter);
__ZM_EXPORT zCVec zPexBH(zPex pex, zCVec ans, double tol, int iter);
__ZM_EXPORT zCVec zPexDKA(zPex pex, zCVec ans, double tol, int iter);

__END_DECLS

#endif /* __ZM_PEX_EQ_H__ */
