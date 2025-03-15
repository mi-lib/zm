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
 * zQuadraticValue() calculates a quadratic form value with the matrix \a q, the vector \a c and \a x,
 * namely, 0.5 xT q x + cT x.
 * \return
 * zQuadraticValue() returns the calculated quadratic form value.
 */
__ZM_EXPORT double zQuadraticValue(const zMat q, const zVec c, const zVec x);

/* TODO: remove description about zQPSolve(). */
/*! \brief convex quadratic programming solver.
 *
 * zQPSolve() solves a quadratic programming with equality constraint conditions
 * to find the vector x that minimizes 0.5 x^T \a q x + \a c^T x subject to
 * \a a x = \a b, where \a q is a positive definite quadratic coefficient matrix,
 * \a c is an offset coefficient vector, \a a and \a b are the coefficient matrix
 * and a constant vector to describe the equality constraint.
 *
 * The result is put into \a ans. If \a cost is not the null pointer, the optimum
 * value is stored in it.
 * \note
 * zQPSolve() does not check if \a q is positive-definite. If not, the solution
 * does not make sense.
 * \return
 * zQPSolve() returns the true value if it succeeds to get the solution, or the
 * false value in failure cases where there is a mismatch of vector and matrix,
 * or it fails to allocate memory.
 * \sa
 * zLCPSolveLemke
 */
#if 0
__ZM_EXPORT bool zQPSolve(const zMat q, const zVec c, const zMat a, const zVec b, zVec ans, double *cost);
#endif

/*! \brief convex quadratic programming solver.
 *
 * zQPSolveLemke() and zQPSolveIP() solve a quadratic programming with inequality constraint conditions
 * to find the vector x that minimizes 0.5 x^T \a q x + \a c^T x subject to \a a x >= \a b and x >= 0
 * based on the Lemke method and the interior point method, respectively.
 * \a q must be positive definite.
 *
 * For the both functions, the result is put into \a ans. If \a cost is not the null pointer, the optimum
 * value is stored in it.
 * \note
 * Neither zQPSolveLemke() nor zQPSolveIP() check if \a q is positive-definite. If not, the solution does
 * not make sense.
 * \return
 * zQPSolveLemke() and zQPSolveIP() return the true value if they succeed to get the solution. If a
 * vector/matrix size mismatch, bad memory allocation, non-feasible solution and infinite solution are
 * found, the false value is returned.
 * \sa
 * zLCPSolveLemke, zLCPSolveIP
 */
__ZM_EXPORT bool zQPSolveLemke(const zMat q, const zVec c, const zMat a, const zVec b, zVec ans, double *cost);
__ZM_EXPORT bool zQPSolveIP(const zMat q, const zVec c, const zMat a, const zVec b, zVec ans, double *cost);

/*! \brief convex quadratic programming solver by active set method.
 *
 * zQPSolveASM() solves a quadratic programming with inequality constraint conditions to find the vector
 * x that minimizes 0.5 x^T \a q x + \a c^T x subject to \a a x >= \a b based on the active set method.
 * \a q must be positive definite.
 *
 * The result is put into \a ans. If \a cost is not the null pointer, the optimum value is stored in it.
 * \note
 * zQPSolveASM() does not check if \a q is positive-definite. If not, the solution does not make sense.
 * \return
 * zQPSolveASM() returns the true value if it succeeds to get the solution. If a vector/matrix size
 * mismatch, bad memory allocation, non-feasible solution and infinite solution are found, it returns
 * the false value.
 */
__ZM_EXPORT bool zQPSolveASM(const zMat q, const zVec c, const zMat a, const zVec b, zVec ans, double *cost);

/*! \brief quadratic programming solver by conjugate gradient method.
 *
 * zCGSolve() is an implementation of quadratic programming without constraint conditions - finds the
 * vector x that minimizes 0.5 xT q x + cT x.
 * \a q is a positive definite quadratic coefficient matrix.
 * \a c is an offset coefficient vector.
 * \a ans is a vector which the answer is put into.
 * \a iter is the maximum iteration number. Iteration would be done up to \a iter times. When \a iter is
 * zero, Z_MAX_ITER_NUM defined in zm_misc.h is chosen instead.
 *
 * Before calling the function, \a ans should be set for a proper initial vector. An unsuitable initial
 * value would lead to an illegal solution.
 * \sa
 * zCGSolve() returns the obtained minimum value.
 */
__ZM_EXPORT double zCGSolve(const zMat q, const zVec c, zVec ans, int iter);

__END_DECLS

#endif /* __ZM_OPT_QP_H__ */
