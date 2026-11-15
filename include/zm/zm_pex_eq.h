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
 * zPexBHDST(), zPexBH() and zPexDKA() solve a polynomial equation
 * \a a=0 based on Bairstow-Hitchcock's method and Durand-Kerner-Aberth's
 * method, respectively.
 *
 * All answers are put into an array pointed by \a ans as complex numbers.
 * \a ans must have enough size.
 * The variable of \a a is scaled in advance by a power of two close to an
 * upper bound of the absolute values of the roots, which is estimated from
 * the coefficients, so that the results do not depend on the scale of the
 * roots. Since the coefficients are multiplied by powers of two, no rounding
 * error is introduced by the scaling.
 * \a tol is a tolerance to terminate iteration relative to the scale of the
 * problem. If zero or a negative value is given, zTOL is used instead.
 * zPexDKA() regards the iteration as converged if the residual |a(z)| at
 * every approximate root z is not greater than \a tol times the scale of the
 * rounding error of the value, sum_k |a_k| |z|^k, after which the correction
 * continues while the residuals decrease (at most ten times) to improve the
 * accuracy of ill-conditioned roots such as close roots.
 * zPexBH() and zPexBHDST() regard the iteration for a quadratic factor
 * x^2 + p x + q as converged if the updates of p and q are not greater than
 * \a tol times the scale of the roots of the factor m = |p| + sqrt(|q|) and
 * its square, respectively. The iteration for each factor starts from
 * several factors with a pair of roots on the unit circle (of the scaled
 * variable) until it converges.
 * zPexBH() refines each quadratic factor for the original polynomial before
 * deflating the polynomial by it, and finally polishes all the roots by
 * Newton=Raphson's method for the original polynomial, in order to prevent
 * the accumulation of errors through the deflations.
 * zPexBHDST() is a destructive version of zPexBH(), which destroys \a a
 * in the computation process. It does neither the refinement nor the
 * polishing, since it does not keep the original polynomial.
 * \a iter is the maximum number of iteration. When it runs over \a iter
 * times, it is aborted. If zero is given for \a iter, Z_MAX_ITER_NUM
 * defined in zm_misc.h is used instead.
 * \note
 * Since Bairstow-Hitchcock's method is an iteration with deflations,
 * zPexBH() may not converge or lose accuracy for a polynomial of a high
 * degree (about fifteen or higher), even with the refinement and the
 * polishing. zPexDKA() is advantageous to it from this viewpoint.
 * \return
 * zPexBHDST(), zPexBH() and zPexDKA() return a pointer \a ans if succeeding.
 * Otherwise, the null pointer is returned.
 */
__ZM_EXPORT zCVec zPexBHDST(zPex a, zCVec ans, double tol, int iter);
__ZM_EXPORT zCVec zPexBH(zPex a, zCVec ans, double tol, int iter);
__ZM_EXPORT zCVec zPexDKA(zPex a, zCVec ans, double tol, int iter);

__END_DECLS

#endif /* __ZM_PEX_EQ_H__ */
