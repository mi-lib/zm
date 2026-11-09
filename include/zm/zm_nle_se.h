/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_nle_se - nonlinear equation: single nonlinear equation solver.
 */

#ifndef __ZM_NLE_SE_H__
#define __ZM_NLE_SE_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief solve a nonlinear equation of a single variable.
 *
 * These functions solve an equation \a f(x) = 0 of a single variable x, and store the answer where
 * \a answer points. \a util is a pointer to a utility data for \a f.
 *
 * zNLE1Bisec() is the bisection method.
 * zNLE1RegulaFalsi() is the Regula-Falsi method with the Illinois modification.
 * zNLE1VDB() is the Van Wijngaarden-Dekker-Brent method.
 * They work in a range that brackets a root, namely, the function has opposite signs at the end points.
 *
 * The Illinois modification halves the value of the function at an end point that is retained, if the
 * new estimate replaces the latest one, which accelerates the convergence of the Regula-Falsi method.
 * It is based on the following paper:
 *  M. Dowell and P. Jarratt, A Modified Regula Falsi Method for Computing the Root of an Equation, BIT,
 *  11, pp. 168--174, 1971.
 *
 * The Van Wijngaarden-Dekker-Brent method combines the bisection, the linear interpolation and the
 * inverse quadratic interpolation, so that it converges superlinearly for a well-behaved function, and
 * is guaranteed to converge with at most about the square of the number of evaluations required by
 * the bisection method. It is based on the procedure zero in the following book:
 *  R. P. Brent, Algorithms for Minimization without Derivatives, Chapter 4, Prentice-Hall, 1973.
 * The tolerance of Brent's procedure is replaced by a half of \a tol_step ( |x| + s0 ), where s0 is
 * defined below.
 *
 * zNLE1Secant() is the secant method, which starts from two initial guesses \a x1 and \a x2. It does not
 * require them to bracket a root, and can find a root outside of them or a root where the function does
 * not change its sign, e.g., a double root, though the convergence to a multiple root is linear. It is
 * suitable to refine a root from good initial guesses, while it might diverge or stagnate.
 *
 * Use zNLE1VDB() to find a root in a range reliably.
 *
 * If \a f(\a x1) and \a f(\a x2) have the same sign, zNLE1Bisec(), zNLE1RegulaFalsi(), and zNLE1VDB()
 * (and zNLE1Secant() if \a f(\a x1) = \a f(\a x2)) search a bracketing range in [ \a x1, \a x2 ]
 * by dividing it into 2^k equal subranges for k=1, 2, ... until the subranges become smaller than
 * \a tol_step ( 2 s0 ), where s0 = max{|\a x1|,|\a x2|} (1 if both are zero). Which root is found depends on the procedure, if the range
 * contains more than one root. Note that a root where the function does not change its sign, e.g., a
 * double root, cannot be found by the search.
 *
 * If the function is not a number or infinity at a trial point, the step to the point is shortened by
 * halving it until the function is evaluated validly. The function has to be valid at \a x1 and \a x2.
 * The iteration succeeds if |\a f(x)| becomes smaller than \a tol_err. It fails if the step becomes
 * smaller than \a tol_step ( |x| + s0 ), or the number of iterations exceeds \a iter. If \a iter is
 * zero, Z_MAX_ITER_NUM is used instead.
 * \return
 * These functions return the number of iterations, i.e., the number of evaluations of \a f except
 * those at \a x1 and \a x2, if they succeed. They return 0 if \a x1 or \a x2 is a root.
 * They return -1 if they fail, i.e., the iteration stagnates, the number of iterations exceeds
 * \a iter, no bracketing range is found, the function is invalid at \a x1 or \a x2, or the secant
 * becomes horizontal (zNLE1Secant()). Even if they fail, the latest estimate is stored where
 * \a answer points.
 */
__ZM_EXPORT int zNLE1Bisec(double (* f)(double,void*), double *answer, double x1, double x2, void *util, double tol_step, double tol_err, int iter);
__ZM_EXPORT int zNLE1Secant(double (* f)(double,void*), double *answer, double x1, double x2, void *util, double tol_step, double tol_err, int iter);
__ZM_EXPORT int zNLE1RegulaFalsi(double (* f)(double,void*), double *answer, double x1, double x2, void *util, double tol_step, double tol_err, int iter);
__ZM_EXPORT int zNLE1VDB(double (* f)(double,void*), double *answer, double x1, double x2, void *util, double tol_step, double tol_err, int iter);

__END_DECLS

#endif /* __ZM_NLE_SE_H__ */
