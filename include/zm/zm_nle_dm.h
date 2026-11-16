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
 * 1. Levenberg-Marquardt's method,
 * 2. Newton-Raphson's method, and
 * 3. Broyden's method.
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
 * the approximate one fails to give a step that decreases the residual, the step stagnates or the
 * iteration makes little progress (see ZM_NLE_NR_DECREASE_RATIO_MIN). Hence,
 * it requires less evaluations of the function than Newton-Raphson's method, if the Jacobian matrix
 * is not given and computed numerically, while the convergence is superlinear instead of quadratic.
 */
#define ZM_NLE_NR_DEFAULT_DAMP_COEFF ( 1.0e-12 )

/*! \brief criterion of the slow progress of Newton-Raphson's method and Broyden's method.
 *
 * An iteration of Newton-Raphson's method or Broyden's method is regarded to make little progress if
 * the squared norm of the residual decreases less than ZM_NLE_NR_DECREASE_RATIO_MIN times itself, i.e.,
 * |f_{k+1}|^2 > ( 1 - ZM_NLE_NR_DECREASE_RATIO_MIN ) |f_k|^2. If more than ZM_NLE_NR_SLOW_ITER_MAX
 * successive iterations make little progress, the iteration is regarded to fail. For Broyden's method, the
 * Jacobian matrix is recomputed exactly instead if it is approximate, and the iteration fails if the
 * next step still makes little progress. It happens, for example, when the line search has to
 * shorten the step extremely along a curved valley of the norm of the residual, where the step is not
 * small but the iteration would take an enormous number of iterations to converge, if at all.
 * It is similar to the criterion of HYBRD in MINPACK, which judges that the iteration is not making
 * good progress from the improvement in the last ten iterations (J. J. More, B. S. Garbow and
 * K. E. Hillstrom, User Guide for MINPACK-1, Argonne National Laboratory Report ANL-80-74, 1980).
 */
#define ZM_NLE_NR_DECREASE_RATIO_MIN ( 1.0e-3 )
#define ZM_NLE_NR_SLOW_ITER_MAX      10

/*! \brief set the coefficient of the damping of Newton-Raphson's method and Broyden's method. */
#define zNLESetNRDampCoeff(nle,coeff) ( (nle)->nr_damp_coeff = (coeff) )

__BEGIN_DECLS

/*! \brief create a simultaneous nonlinear equations solver.
 *
 * zNLECreate() creates a simultaneous nonlinear equation solver \a nle. \a nv is the size of the
 * answer vector. \a ne is the number of equations. \a f is a pointer to the equation function to
 * be solved, '\a f(x)=0'. \a jac is a pointer to the function to compute the Jacobian matrix of \a f.
 * If it is the null pointer, the Jacobian matrix is computed by the central difference, where the
 * one-sided difference is used instead if the residual is not finite on one side, in the same way with
 * ZM_OPT_DM_EPS.
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
 *  "LM" for Levenberg-Marquardt's method,
 *  "NR" for Newton-Raphson's method, and
 *  "Broyden" for Broyden's method.
 * The steepest descent method, the variable metric method and the conjugate gradient method are
 * not available, since they
 * need the Jacobian matrix to compute the gradient vector of the squared norm of the residual at
 * every evaluation, while they do not utilize it to approximate the Hessian matrix.
 */
__ZM_EXPORT zNLE *zNLEAssignLM(zNLE *nle, const char *dampmethod, const char *stepmethod, double scale);
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
 * becomes smaller than \a tol_step ( |\a var|inf + s0 ),
 * where s0 is the maximum absolute component of the initial
 * guess. If the initial guess is zero, s0 is that of the
 * first step for Newton-Raphson's method and Broyden's
 * method (1 if it is zero as well), or 1 for
 * Levenberg-Marquardt's method,
 * namely, the iteration stagnates at a point where the
 * residual is not small enough, e.g., a local minimum of
 * the norm of the residual.
 * Newton-Raphson's method and Broyden's method also fail
 * if the residual decreases little in successive
 * iterations (see ZM_NLE_NR_DECREASE_RATIO_MIN).
 * \a iter is the maximum iteration number. Iteration
 * would be done up to \a iter times. When \a iter is
 * zero, Z_MAX_ITER_NUM defined in zm_misc.h is chosen
 * instead.
 * A half of the weighted squared norm of the residual
 * vector at the answer is stored where \a err points, if
 * it is not the null pointer.
 * \return the number of iterations if it succeeds, or -1
 * if it fails, i.e., the iteration stagnates or makes
 * little progress, the number of iterations exceeds
 * \a iter, or something is wrong.
 */
#define zNLESolve(nle,var,util,tol_step,tol_err,iter,err) (nle)->_solve( nle, var, util, tol_step, tol_err, iter, err )

__END_DECLS

/*! \} */

#endif /* __ZM_NLE_DM_H__ */
