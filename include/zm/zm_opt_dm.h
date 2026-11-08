/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_dm - optimization tools: descent method.
 */

#ifndef __ZM_OPT_DM_H__
#define __ZM_OPT_DM_H__

/* NOTE: never include this header file in user programs. */

/* ********************************************************** */
/*! \defgroup descent method for unconstrained nonlinear optimization
 * \{ *//* ************************************************** */

/*! \struct zOptDM
 * \brief descent method for unconstrained nonlinear optimization
 *
 * zOptDM is an unconstrained nonlinear optimization solver
 * class based on descent method.
 *
 * Available methods include
 * 1. steepest descent method,
 * 2. Levenberg-Marquardt's method,
 * 3. variable metric method (quasi Newton's method, secant method), and
 * 4. conjugate gradient method (Fletcher-Reeves' method).
 */
ZDEF_STRUCT( __ZM_CLASS_EXPORT, zOptDM ){
  double (* eval)(const zVec,void*);
  zVec (* grad)(const zVec,zVec,void*);
  zMat (* hess)(const zVec,zMat,void*);
  double eps_min; /* ratio of the step of numerical differentiation to the magnitude of the variable */
  double tol_grad; /* ratio of the gradient vector to the initial one for convergence */
  /*! \cond */
  zVec (* _grad)(zOptDM*,const zVec,zVec,void*);
  zMat (* _hess)(zOptDM*,const zVec,zMat,void*);
  zVec (* _vec)(zOptDM*,const zVec,zVec,void*);
  zVec (* _step)(zOptDM*,const zVec,double,void*,double,double*);
  double (* _damp)(zOptDM*,void*);
  zMat (* _update_h)(zOptDM*,int);
  double _scale; /* scale of search range for MT and DS, and scale of stride for fixed method. */
  zVec _x, _d, _g;
  zVec _p, _q, _r;
  zMat _h;
  zMat _hw; /* workspace for the damping coefficient of Levenberg-Marquardt's method */
  zIndex _idx;
  double _b;  /* for Fletcher-Reeves method */
  double _df; /* decrease factor */
  double _cf; /* curvature factor */
  /*! \endcond */
};

/*! \brief default threshold of the gradient vector for convergence of descent method.
 *
 * zOptDMSolve() regards the iteration converged if the maximum absolute component of the gradient
 * vector becomes smaller than \a tol_grad times that at the initial point.
 * \a tol_grad of zOptDM is initialized by ZM_OPT_DM_GRAD_MIN, and is set by zOptDMSetTolGrad().
 * Since a line search judges a decrease of the function value, the gradient vector can be reduced
 * only to the order of the square root of the machine epsilon relatively.
 * If \a tol_grad is zero, the criterion with the gradient vector is disabled except that the gradient
 * vector exactly vanishes, and the iteration terminates only when the step becomes small or the
 * number of iterations exceeds the upper limit. It is useful when the function is not expected to
 * be minimized but to be zero, e.g., the squared norm of the residual of nonlinear equations, where
 * the gradient vector can be small even far from the solution if the Jacobian matrix is nearly
 * singular.
 */
#define ZM_OPT_DM_GRAD_MIN  ( 1.0e-8 )

/*! \brief set the threshold of the gradient vector for convergence relative to the initial one. */
#define zOptDMSetTolGrad(opt,t) ( (opt)->tol_grad = (t) )

/*! \brief default ratio of the step of numerical differentiation to the magnitude of the variable.
 *
 * If the gradient vector or the Hessian matrix is not given to zOptDMCreate(), it is computed by the
 * central difference. The step for the i-th component of the variable x is \a eps_min s_i for the gradient
 * vector and the Hessian matrix computed from the gradient vector, and \a eps_min^(3/4) s_i for the Hessian
 * matrix computed from the function value, where s_i is the magnitude of x_i defined with
 * ZM_OPT_DM_DIFF_SCALE_RATIO.
 * \a eps_min of zOptDM is initialized by ZM_OPT_DM_EPS, which is suitable if the function is computed with
 * the accuracy of the machine epsilon. It should be larger for a function with noise, which is set by
 * zOptDMSetEPSMin().
 * Note that the numerical differentiation is not accurate if the magnitude of the variable does not
 * represent the scale of variation of the function. Then, the gradient vector should be given, or the
 * variable should be scaled.
 */
#define ZM_OPT_DM_EPS ( 1.0e-5 )

/*! \brief ratio of the lower bound of the magnitude of a component of the variable to that of the
 * variable for numerical differentiation.
 *
 * The step of numerical differentiation for the i-th component of the variable x is proportional to
 * the magnitude s_i of x_i. Since |x_i| itself vanishes when x_i crosses zero, it is bounded below
 * relatively to the magnitude of the whole variable |x|inf (the maximum absolute value of the
 * components), namely, s_i = max{ |x_i|, r |x|inf }, where r is ZM_OPT_DM_DIFF_SCALE_RATIO for the
 * gradient vector and the Hessian matrix computed from the gradient vector. For the Hessian matrix
 * computed from the function value, r is the square root of ZM_OPT_DM_DIFF_SCALE_RATIO, since the
 * round-off error of the second-order difference is inversely proportional to the square of the step,
 * while that of the first-order difference is inversely proportional to the step; the maximum
 * amplification of the round-off error by the lower bound is equalized for the both. s_i = 1 if x = 0.
 */
#define ZM_OPT_DM_DIFF_SCALE_RATIO ( 1.0e-2 )

/*! \brief set the ratio of the step of numerical differentiation. */
#define zOptDMSetEPSMin(opt,e) ( (opt)->eps_min = (e) )

/*! \brief default scale of search range for More-Thuente's method and Dennis-Schnabel's method. */
#define ZM_OPT_DM_DEFAULT_SCALE  ( 1.0e3 )
/*! \brief default step length of the line search with fixed step. */
#define ZM_OPT_DM_DEFAULT_STEP     ( 1.0 )

__BEGIN_DECLS

/*! \brief create an unconstrained nonlinear optimization solver based on descent method.
 *
 * zOptDMCreate() creates a solver \a opt of unconstrained nonlinear optimization problems based on
 * descent method. \a dim is the dimension of the variable vector. \a eval is the function to be
 * minimized. \a grad and \a hess are functions to compute the gradient vector and the Hessian matrix
 * of \a eval, respectively. If they are the null pointers, they are computed numerically.
 * The steepest descent method with More-Thuente's line search is assigned by default.
 * \return
 * zOptDMCreate() returns a pointer \a opt if it succeeds. Otherwise, it returns the null pointer.
 */
__ZM_EXPORT zOptDM *zOptDMCreate(zOptDM *opt, int dim, double (* eval)(const zVec,void*), zVec (* grad)(const zVec,zVec,void*), zMat (* hess)(const zVec,zMat,void*));

/*! \brief assign a descent method.
 *
 * zOptDMAssignSD(), zOptDMAssignLM(), zOptDMAssignVM() and zOptDMAssignCG() assign the steepest
 * descent method, Levenberg-Marquardt's method, the variable metric (quasi-Newton) method and the
 * conjugate gradient method (Fletcher-Reeves' method) to \a opt, respectively.
 *
 * \a stepmethod specifies the method to compute the step length as follows:
 *  "MT" for More-Thuente's line search,
 *  "DS" for Dennis-Schnabel's backtracking line search,
 *  "Brent" for Brent's strict line search, and
 *  "fixed" for the fixed step length without line search.
 * If it is the null pointer or another string, More-Thuente's line search is assigned, except that
 * Brent's line search is assigned for zOptDMAssignCG() if it is the null pointer.
 *
 * The step methods are suitable for the following cases:
 *  - "MT" requires the function value and the gradient vector at every trial step.
 *    It is the most robust. Use it if unsure. It guarantees the strong Wolfe's conditions, which keep
 *    the quasi-Hessian matrix of the variable metric method positive definite and the direction of the
 *    conjugate gradient method descent. It can also extrapolate the step length.
 *    It is not suitable if the gradient vector is expensive to compute, e.g., computed numerically.
 *  - "DS" requires only the function value at every trial step.
 *    It is suitable for Newton-like directions, e.g., the variable metric method, for which the unit
 *    step length is natural, when the gradient vector is expensive. It is the cheapest as long as the
 *    full step is accepted.
 *    It is not suitable for the steepest descent method or the conjugate gradient method, for which
 *    the step length depends on the magnitude of the direction vector, since it cannot extrapolate.
 *  - "Brent" requires only the function value at every trial step, but many trials.
 *    It is suitable if each iteration is expensive, e.g., Levenberg-Marquardt's method that computes
 *    the Hessian matrix, since it minimizes the number of iterations. It is also suitable for the
 *    conjugate gradient method on nearly quadratic functions, which converges within the dimension
 *    with the strict line search.
 *    It is not suitable if the function is expensive to evaluate.
 *  - "fixed" requires the function value and the gradient vector once for each iteration.
 *    It is suitable for real-time applications with a known proper step length, where the computation
 *    time for each iteration should be constant.
 *    It is not suitable in general, since the iteration may diverge.
 *
 * \a scale is the upper bound of the components of a step for "MT", "DS" and "Brent", and the step
 * length for "fixed". If it is zero, ZM_OPT_DM_DEFAULT_SCALE is used for "MT", "DS" and "Brent", and
 * ZM_OPT_DM_DEFAULT_STEP is used for "fixed".
 *
 * \a dampmethod for zOptDMAssignLM() specifies the damping coefficient mu added to the diagonal
 * components of the Hessian matrix H, namely,
 *  "grad" for mu = |g|, where g is the gradient vector, and
 *  "dec" for mu = min{ |g|, g^T H^-1 g }, where g^T H^-1 g is the Newton decrement that estimates
 *  twice the deviation of the function from the minimum. It coincides with |F|^2 for the least
 *  square problem of F(x)=0 with the Gauss-Newton approximation of the Hessian matrix if the
 *  Jacobian matrix of F is regular. It requires to solve an additional linear equation at every
 *  iteration, and a workspace is allocated internally when it is assigned at the first time.
 * "grad" is assigned if it is the null pointer or another string.
 *
 * \a updatemethod for zOptDMAssignVM() specifies the update formula of the quasi-Hessian matrix,
 * namely, "DFP" for Davidon-Fletcher-Powell's formula and others for Broyden-Fletcher-Goldfarb-Shanno's.
 * \return
 * These functions return a pointer \a opt. zOptDMAssignLM() returns the null pointer if it fails
 * to allocate the workspace.
 */
__ZM_EXPORT zOptDM *zOptDMAssignSD(zOptDM *opt, const char *stepmethod, double scale);
__ZM_EXPORT zOptDM *zOptDMAssignLM(zOptDM *opt, const char *dampmethod, const char *stepmethod, double scale);
__ZM_EXPORT zOptDM *zOptDMAssignVM(zOptDM *opt, const char *updatemethod, const char *stepmethod, double scale);
__ZM_EXPORT zOptDM *zOptDMAssignCG(zOptDM *opt, const char *stepmethod, double scale);

/*! \brief destroy an unconstrained nonlinear optimization solver based on descent method. */
__ZM_EXPORT void zOptDMDestroy(zOptDM *opt);

/*! \brief solve an unconstrained nonlinear optimization problem by descent method.
 *
 * zOptDMSolve() minimizes the function assigned to \a opt from the initial guess \a var, and stores
 * the solution into \a var. \a util is a pointer to the utility data for the function.
 * The iteration terminates if one of the following conditions holds:
 *  - the maximum absolute component of the gradient vector becomes smaller than \a tol_grad of \a opt
 *    (ZM_OPT_DM_GRAD_MIN by default) times that at the initial point,
 *  - the maximum absolute component of the step becomes smaller than \a tol ( 1 + |\a var|inf ), or
 *  - the number of iterations exceeds \a iter. If \a iter is zero, Z_MAX_ITER_NUM is used instead.
 * The value of the function at the solution is stored where \a eval points, if it is not the null
 * pointer.
 * \return
 * zOptDMSolve() returns the number of iterations if the iteration converges. It returns 0 if the gradient
 * vector at the initial point is zero. Otherwise, i.e., if the number of iterations exceeds \a iter or the
 * line search fails, it returns -1.
 */
__ZM_EXPORT int zOptDMSolve(zOptDM *opt, zVec var, void *util, double tol, int iter, double *eval);

/*! \} */

__END_DECLS

#endif /* __ZM_OPT_DM_H__ */
