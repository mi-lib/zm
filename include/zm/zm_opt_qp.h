/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_qp - optimization tools: quadratic programming.
 */

#ifndef __ZM_OPT_QP_H__
#define __ZM_OPT_QP_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief calculation of a quadratic value.
 *
 * zQuadraticValue() calculates a quadratic form value with
 * the matrix \a q, the vector \a c and \a x, namely,
 * 0.5 xT q x + cT x.
 * \return
 * zQuadraticValue() returns the calculated quadratic form value.
 */
__EXPORT double zQuadraticValue(zMat q, zVec c, zVec x);

/*! \brief convex quadratic programming solver.
 *
 * zQPSolve() solves a quadratic programming with equality
 * constraint conditions to find the vector x which
 * minimizes 0.5 xT q x + cT x subject to a x = b
 * where \a q is a positive definite quadratic coefficient matrix.
 * \a c is an offset coefficient vector.
 * \a a and \a b are the coefficient matrix and constant vector
 * to describe equality constraint.
 *
 * zQPSolveLemke() solves a quadratic programming with inequality
 * constraint to find the vector x which minimizes 0.5 xT q x + cT x
 * subject to a x >= b and x >= 0 based on Lemke method. For this
 * function, \a q must be positive definite.
 *
 * zQPSolveLemke() solves a quadratic programming with inequality
 * constraint to find the vector x which minimizes 0.5 xT q x + cT x
 * subject to a x >= b and x >= 0 based on interior point method.
 * For this function, \a q must be positive definite.
 *
 * For these functions, the result is put into \a ans, and
 * if \a cost is not the null pointer, the optimum value
 * is put into it.
 * \return
 * zQPSolve() and zQPSolveLemke() returns the true value if
 * it succeeds to get the solution, or the the false value if
 * there are vector/matrix size mismatch, bad memory allocation,
 * non-feasible solution and infinite solution are found.
 * \sa
 * zLCPSolveLemke
 */
__EXPORT bool zQPSolve(zMat q, zVec c, zMat a, zVec b, zVec ans, double *cost);
__EXPORT bool zQPSolveLemke(zMat q, zVec c, zMat a, zVec b, zVec ans, double *cost);
__EXPORT bool zQPSolveIP(zMat q, zVec c, zMat a, zVec b, zVec ans, double *cost);

/*! \brief quadratic programming solver by active-set method
 *
 * zQPSolveASM() solves a quadratic programming with inequality
 * constraint to find the vector x which minimizes 0.5 xT q x + cT x
 * subject to a x >= b and x >= 0 based on active-set method.
 */
__EXPORT bool zQPSolveASM(zMat q, zVec c, zMat a, zVec b, zVec ans, double *cost, zVec init(zMat,zVec,zVec,void*), void *util);

/*! \brief quadratic programming solver by conjugate gradient method.
 *
 * zCGSolve() is an implementation of quadratic programming
 * without constraint conditions - finds the vector x
 * which minimizes 0.5 xT q x + cT x.
 * \a q is a positive definite quadratic coefficient matrix.
 * \a c is an offset coefficient vector.
 * \a ans is a vector which the answer is put into.
 * \a iter is the maximum iteration number. Iteration would be
 * done up to \a iter times. When \a iter is zero, Z_MAX_ITER_NUM
 * defined in zm_misc.h is chosen instead.
 *
 * Before calling the function, \a ans should be set for a proper
 * initial vector. An unsuitable initial value would lead to an
 * illegal solution.
 * \sa
 * zCGSolve() returns the obtained minimum value.
 */
__EXPORT double zCGSolve(zMat q, zVec c, zVec ans, int iter);

__END_DECLS

#endif /* __ZM_OPT_QP_H__ */
