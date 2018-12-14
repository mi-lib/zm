/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_lp - optimization tools: linear programming.
 */

#ifndef __ZM_OPT_LP_H__
#define __ZM_OPT_LP_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* METHOD:
 * zLPIneq2Std, zLPUnb2Std
 * - convert linear programming problem from inequality
 *   constraint to standard equation form.
 *
 * The standard equation form of linear programming problem
 * is as follows.
 *   'c^T x' -> minimum
 *   subject to 'a x = b' and 'x >= 0'
 * #
 * 'zLPIneq2Std()' converts linear programming problem from
 * the following inequality constraint form:
 *   'c^T x' -> minimum
 *   subject to 'a x <= b' and 'x >= 0'
 * to the standard equation form, where the set of 'a', 'b'
 * and 'c' denotes the original inequality form and the set
 * of 'as', 'bs' and 'cs' does the standard form. Namely:
 *   'cs^T x' -> minimum
 *   subject to 'as x = bs' and 'x >= 0'
 * #
 * 'zLPUnb2Std()' converts linear programming problem from
 * the following inequality constraint form with unbounded
 * variables:
 *   'c^T x' -> minimum
 *   subject to 'a x <= b'
 * to the standard equation form with positive variables.
 * The set of 'as', 'bs' and 'cs' of this function are also
 * for the following standard form:
 *   'cs^T x' -> minimum
 *   subject to 'as x = bs' and 'x >= 0'
 * [RETURN VALUE]
 * 'zLPIneq2Std()' and 'zLPUnb2Std()' return the true value
 * if succeeding. Or, the false value is returned when they
 * fail to allocate working memory.
 */
__EXPORT bool zLPIneq2Std(zMat a, zVec c, zVec x, zMat *as, zVec *cs, zVec *xs);
__EXPORT bool zLPUnb2Std(zMat a, zVec c, zVec x, zMat *as, zVec *cs, zVec *xs);

/* METHOD:
 * zLPSolveSimplex
 * - linear programming solver with the simplex method.
 *
 * 'zLPSolveSimplex()' solves a linear programming problem
 * by dual phase simplex method (1947 G. B. Danzig) which forms:
 *   'c^T x' -> minimum
 *   subject to 'a x = b' and 'x >= 0'
 * where 'c' is a coefficient vector of the cost funtion,
 * 'a' and 'b' are a coefficient matrix and constant vector
 * to describe the constraints,
 * 'ans' is a vector which the answer will be put into,
 * and 'index' is an index vector for bases.
 * #
 * The result optimum value is stored where 'cost' points,
 * unless 'cost' is the null pointer.
 * [RETURN VALUE]
 * 'zLPSolveSimplex()' returns the true value if succeeding
 * to get optimum (i.e. minimum) cost, or the the false
 * value if there are vector/matrix size mismatch, bad
 * memory allocation, non-feasible solution and infinite
 * solution found.
 */
__EXPORT bool zLPSolveSimplex(zMat a, zVec b, zVec c, zVec ans, double *cost);

/*! \brief find a feasible base under Ax=b and x>=0 based on simplex method.
 */
__EXPORT bool zLPFeasibleBase(zMat a, zVec b, zVec base);

/* METHOD:
 * zLPSolvePDIP_PC
 * - linear programming solver based on primal-dual
 *   interior-point method.
 *
 * 'zLPSolvePDIP_PC()' solves a linear programming problem
 * which forms:
 *   'c^T x' -> minimum
 *   subject to 'a x = b' and 'x >= 0'
 * where 'c' is a coefficient vector of the cost funtion,
 * 'a' and 'b' are a coefficient matrix and constant vector
 * to describe the constraints, and 'ans' is a vector which
 * the answer will be put into.
 * Note that all constraints are equations. Programmers
 * have to add slack variables in order to cope with
 * inequality constraints.
 * The resultant minimum cost value is put where 'cost' points
 * unless it is not the null pointer.
 * #
 * 'zLPSolvePDIP_PC()' is based on primal-dual interior-point
 * method with Mehrotra s predictor-corrector (1992 S. Mehrotra).
 * [RETURN VALUE]
 * 'zLPSolvePDIP_PC()' returns the true value if it
 * succeeds to obtain the minimum value. Or, the false
 * value is returned when it fails to allocate enough
 * working memory or to solve the problem if a set of
 * constraints is unbounded.
 */
__EXPORT bool zLPSolvePDIP_PC(zMat a, zVec b, zVec c, zVec x, double *cost);

__END_DECLS

#endif /* __ZM_OPT_LP_H__ */
