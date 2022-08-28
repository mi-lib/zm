/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_lp - optimization tools: linear programming.
 */

#ifndef __ZM_OPT_LP_H__
#define __ZM_OPT_LP_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief simplex method tableau. */
typedef struct{
  zMat a;    /*!< coefficient matrix */
  zVec b;    /*!< constraint vector */
  zVec c;    /*!< cost vector */
  double d;  /*!< cost */
  zIndex ib; /*!< base index */
  zIndex in; /*!< non-base index */
  zIndex ir; /*!< relevant constraint index */
} zLPTableau;

/*! \brief create initial simplex tableau with slack variables. */
__EXPORT bool zLPTableauCreate(zLPTableau *tab, zMat a, zVec b);
/*! \brief destroy simplex tableau. */
__EXPORT void zLPTableauDestroy(zLPTableau *tab);
/*! \brief simplex method for initialized tableau. */
__EXPORT bool zLPTableauSimplex(zLPTableau *tab);
/*! \brief find initial feasible base for the second stage from tableau. */
__EXPORT bool zLPTableauFindBase(zLPTableau *tab);

/*! \brief convert linear programming problem from inequality
 * constraint to standard equation form.
 *
 * The standard equation form of linear programming problem
 * is as follows.
 *   \a c^T x -> minimum
 *   subject to \a a x = \a b and x >= 0
 *
 * zLPIneq2Std() converts linear programming problem from
 * the following inequality constraint form:
 *   \a c^T x -> minimum
 *   subject to \a a x <= \a b and x >= 0
 * to the standard equation form, where the set of \a a, \a b
 * and \a c denotes the original inequality form and the set
 * of \a as, \a bs and \a cs does the standard form. Namely:
 *   \a cs^T \a x -> minimum
 *   subject to \a as x = \a bs and x >= 0
 *
 * zLPUnb2Std() converts linear programming problem from the
 * following inequality constraint form with unbounded
 * variables:
 *   \a c^T x -> minimum
 *   subject to \a a x <= \a b
 * to the standard equation form with positive variables.
 * The set of \a as, \a bs and \a cs of this function are also
 * for the following standard form:
 *   \a cs^T x -> minimum
 *   subject to \a as x = \a bs and x >= 0
 * \return
 * zLPIneq2Std() and zLPUnb2Std() return the true value
 * if succeeding. Or, the false value is returned when they
 * fail to allocate working memory.
 */
__EXPORT bool zLPIneq2Std(zMat a, zVec c, zVec x, zMat *as, zVec *cs, zVec *xs);
__EXPORT bool zLPUnb2Std(zMat a, zVec c, zVec x, zMat *as, zVec *cs, zVec *xs);

/*! \brief linear programming solver with the simplex method.
 *
 * zLPSolveSimplex() solves a linear programming problem
 * by dual phase simplex method (1947 G. B. Danzig) which forms:
 *   \a c^T x -> minimum
 *   subject to \a a x = \a b and x >= 0
 * where \a c is a coefficient vector of the cost funtion,
 * \a a and \a b are a coefficient matrix and constant vector
 * to describe the constraints,
 * \a ans is a vector which the answer will be put into,
 * and \a index is an index vector for bases.
 *
 * The result optimum value is stored where \a cost points,
 * unless \a cost is the null pointer.
 * \return
 * zLPSolveSimplex() returns the true value if it succeeds
 * to get optimum (i.e. minimum) cost, or the the false
 * value if there are vector/matrix size mismatch, bad
 * memory allocation, non-feasible solution and infinite
 * solution found.
 */
__EXPORT bool zLPSolveSimplex(zMat a, zVec b, zVec c, zVec ans, double *cost);

/*! \brief find a feasible base under Ax=b and x>=0 based on simplex method.
 */
__EXPORT bool zLPFeasibleBase(zMat a, zVec b, zVec base);

/*! \brief linear programming solver based on primal-dual
 * interior-point method.
 *
 * zLPSolvePDIP_PC() solves a linear programming problem
 * which forms:
 *   \a c^T x -> minimum
 *   subject to \a a x = \a b and x >= 0
 * where \a c is a coefficient vector of the cost funtion,
 * \a a and \a b are a coefficient matrix and constant vector
 * to describe the constraints, and \a ans is a vector which
 * the answer will be put into.
 * Note that all constraints are equations. Programmers
 * have to add slack variables in order to cope with
 * inequality constraints.
 * The resultant minimum cost value is put where \a cost points
 * unless it is not the null pointer.
 *
 * zLPSolvePDIP_PC() is based on primal-dual interior-point
 * method with Mehrotra s predictor-corrector (1992 S. Mehrotra).
 * \return
 * zLPSolvePDIP_PC() returns the true value if it succeeds
 * to obtain the minimum value. Or, the false value is
 * returned when it fails to allocate enough working memory
 * or to solve the problem if a set of constraints is unbounded.
 */
__EXPORT bool zLPSolvePDIP_PC(zMat a, zVec b, zVec c, zVec x, double *cost);

__END_DECLS

#endif /* __ZM_OPT_LP_H__ */
