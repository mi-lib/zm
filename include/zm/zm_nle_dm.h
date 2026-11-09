/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_nle_dm - nonlinear equation:
 * simultaneous nonlinear equation solver based on descent method.
 */

#ifndef __ZM_NLE_DM_H__
#define __ZM_NLE_DM_H__

/* NOTE: never include this header file in user programs. */

#include <zm/zm_opt.h>

/* ********************************************************** */
/*! \defgroup descent method for simultaneous nonlinear equations
 * \{ *//* ************************************************** */

/*! \struct zNLE
 * \brief descent method for simultaneous nonlinear equations
 *
 * zNLE is simultaneous nonlinear equation solver class based on descent method.
 * Available methods include
 * 1. steepest descent method,
 * 2. Levenberg-Marquardt's method,
 * 3. variable metric method (quasi Newton's method, secant method),
 * 4. conjugate gradient method (Fletcher-Reeves' method).
 * 5. Newton-Raphson's method, and
 * 6. Broyden's method.
 */
ZDEF_STRUCT( __ZM_CLASS_EXPORT, zNLE ){
  zVec (*f)(const zVec,zVec,void*);
  zMat (*jac)(const zVec,zMat,void*);
  void *util;
  zVec wn; /* relative weights on the components of the step of Newton-Raphson's method and Broyden's method */
  zVec we; /* weights on the components of the residual */
  double nr_damp_coeff; /* coefficient of the damping of Newton-Raphson's method and Broyden's method */
  /*! \cond */
  zMat (*_jac)(zNLE*,zVec,zMat,void*);
  int ( *_solve)(zNLE*,zVec,void*,double,double,int,double*);
  zVec _f;
  zVec _fw;
  zVec _fp;
  zVec _adg, _prg; /* for numerical gradient computation */
  zMat _j;
  zOptDM _opt;
  /*! \endcond */
};

/*! \brief default coefficient of the damping of Newton-Raphson's method and Broyden's method.
 *
 * Newton-Raphson's method of zNLE computes the step d by the SR-inverse of the Jacobian matrix J,
 * namely, ( J^T W J + c diag(wn) ) d = J^T W f, where f is the residual vector, W = diag(we), and
 * c is nr_damp_coeff times the maximum diagonal component of J^T W J. Since c is relative to the
 * magnitude of J^T W J, the step is invariant with respect to the scales of the variable and the
 * residual. The damping is a regularization against a singular Jacobian matrix, and nr_damp_coeff
 * should be small enough not to slow down the convergence, while it should not be too small to keep
 * J^T W J + c diag(wn) regular, in particular, for an underdetermined system, i.e., the number of
 * equations is less than that of the variables.
 * Since c is small, wn works as relative weights on the components of the step. In particular, for
 * an underdetermined system, the step is the one that minimizes d^T diag(wn) d among those satisfying
 * J d = f, namely, a component with a larger weight is moved less. wn is initialized to be 1 for all
 * the components. Note that the convergence might be slow if the ratios of the weights are extreme
 * and components with large weights have to be moved.
 * nr_damp_coeff is initialized by ZM_NLE_NR_DEFAULT_DAMP_COEFF in zNLECreate(), and set by
 * zNLESetNRDampCoeff().
 * Broyden's method of zNLE computes the step in the same way, where J is replaced by its approximation
 * updated by Broyden's rule. The exact Jacobian matrix is computed only at the initial point and when
 * the approximate one fails to give a step that decreases the residual or the step stagnates. Hence,
 * it requires less evaluations of the function than Newton-Raphson's method, if the Jacobian matrix
 * is not given and computed numerically, while the convergence is superlinear instead of quadratic.
 */
#define ZM_NLE_NR_DEFAULT_DAMP_COEFF ( 1.0e-12 )

/*! \brief set the coefficient of the damping of Newton-Raphson's method and Broyden's method. */
#define zNLESetNRDampCoeff(nle,coeff) ( (nle)->nr_damp_coeff = (coeff) )

__BEGIN_DECLS

/*! \brief create a simultaneous nonlinear equations solver.
 *
 * zNLECreate() creates a simultaneous nonlinear equation solver \a nle. \a nv is the size of the
 * answer vector. \a ne is the number of equations. \a f is a pointer to the equation function to
 * be solved, '\a f(x)=0'.
 * \return
 * zNLECreate() returns a pointer \a nle if it succeeds. Otherwise, it returns the null pointer.
 * \sa
 * zNLEAssign
 */
__ZM_EXPORT zNLE *zNLECreate(zNLE *nle, int nv, int ne, zVec (*f)(const zVec,zVec,void*), zMat (*jac)(const zVec,zMat,void*));

/*! \brief destroy a simultaneous nonlinear equation solver.
 *
 * zNLEDestroy() destroys a solver instance \a nle.
 * \sa
 * zNLECreate
 */
__ZM_EXPORT void zNLEDestroy(zNLE *nle);

/*! \brief set the ratio of the step of numerical differentiation.
 *
 * zNLESetEPSMin() sets the ratio of the step of numerical differentiation to the magnitude of the
 * variable to compute the Jacobian matrix, if it is not given to zNLECreate().
 * It is initialized by ZM_OPT_DM_EPS. See ZM_OPT_DM_EPS for the detail.
 */
#define zNLESetEPSMin(nle,e) zOptDMSetEPSMin( &(nle)->_opt, e )

/*! \brief assign simultaneous nonlinear equation solver.
 *
 * zNLEAssign() assignes a numerical method to solve
 * simultaneous nonlinear equations to \a nle.
 * Currently, the following methods are available.
 *  "SD" for steepest descent method,
 *  "LM" for Levenberg-Marquardt's method,
 *  "VM" for variable metric method,
 *  "CG" for conjugate gradient metric method,
 *  "NR" for Newton-Raphson's method, and
 *  "Broyden" for Broyden's method.
 */
__ZM_EXPORT zNLE *zNLEAssignSD(zNLE *nle, const char *stepmethod, double scale);
__ZM_EXPORT zNLE *zNLEAssignLM(zNLE *nle, const char *dampmethod, const char *stepmethod, double scale);
__ZM_EXPORT zNLE *zNLEAssignVM(zNLE *nle, const char *updatemethod, const char *stepmethod, double scale);
__ZM_EXPORT zNLE *zNLEAssignCG(zNLE *nle, const char *stepmethod, double scale);
__ZM_EXPORT zNLE *zNLEAssignNR(zNLE *nle);
__ZM_EXPORT zNLE *zNLEAssignBroyden(zNLE *nle);

/*! \brief solve simultaneous nonlinear equations.
 *
 * zNLESolve() solves simultaneous nonlinear equations.
 * The answer is put into \a var. Note that \a var should
 * be initialized for suitable iteration.
 * \a util is for a programmer's utility to attach any
 * type of data chunk.
 * The iteration succeeds if the maximum absolute component
 * of the residual vector becomes smaller than \a tol_err.
 * It fails if the maximum absolute component of the step
 * becomes smaller than \a tol_step ( 1 + |\a var|inf ),
 * namely, the iteration stagnates at a point where the
 * residual is not small enough, e.g., a local minimum of
 * the norm of the residual.
 * \a iter is the maximum iteration number. Iteration
 * would be done up to \a iter times. When \a iter is
 * zero, Z_MAX_ITER_NUM defined in zm_misc.h is chosen
 * instead.
 * A half of the weighted squared norm of the residual
 * vector at the answer is stored where \a err points, if
 * it is not the null pointer.
 * \return the number of iterations if it succeeds, or -1
 * if it fails, i.e., the iteration stagnates, the number
 * of iterations exceeds \a iter, or something is wrong.
 */
#define zNLESolve(nle,var,util,tol_step,tol_err,iter,err) (nle)->_solve( nle, var, util, tol_step, tol_err, iter, err )

__END_DECLS

/*! \} */

#endif /* __ZM_NLE_DM_H__ */
