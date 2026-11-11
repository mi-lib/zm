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
  double lm_damp_coeff; /* coefficient of the damping of Levenberg-Marquardt's method */
  /*! \cond */
  zVec (* _grad)(zOptDM*,const zVec,zVec,void*);
  zMat (* _hess)(zOptDM*,const zVec,zMat,void*);
  zVec (* _vec)(zOptDM*,const zVec,zVec,void*);
  zVec (* _step)(zOptDM*,const zVec,double,void*,double,double*);
  double (* _damp)(zOptDM*,void*);
  zMat (* _update_h)(zOptDM*,int);
  double _scale; /* scale of search range relative to _xscale for MT, DS and Brent, and step length for fixed method. */
  double _xscale; /* scale of the variable for the judgement of the step, namely, |x0|inf (1 if x0 = 0) */
  double _g0;   /* norm of the initial gradient vector for Levenberg-Marquardt's method */
  double _dec0; /* initial Newton decrement for Levenberg-Marquardt's method */
  bool _to_start; /* true at the initial point or at a restart, where the information of the previous step is not available */
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
 *
 * The function (or the gradient vector for the Hessian matrix computed from it) is evaluated on both
 * sides of the variable along each component. If both values are finite, the central difference is used.
 * If only one of them is finite, the one-sided difference on the finite side is used instead. If neither
 * is finite, the result is not finite, and zOptDMSolve() fails. For the Hessian matrix computed from the
 * function value, the off-diagonal components also use the values at the points displaced along two
 * components. If some of them are not finite, the one-sided difference on a side where they are finite
 * is used instead.
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

/*! \brief coefficient of Armijo's condition (sufficient decrease) for line searches.
 *
 * A step length a along a descent direction is accepted if phi(a) <= phi(0) + c a phi'(0), where phi
 * is the function along the direction and c is ZM_OPT_DM_ARMIJO_COEFF.
 */
#define ZM_OPT_DM_ARMIJO_COEFF ( 1.0e-4 )

/*! \brief default coefficient of the damping of Levenberg-Marquardt's method.
 *
 * Levenberg-Marquardt's method of zOptDM adds the damping coefficient mu to the diagonal components of
 * the Hessian matrix H, where mu = c lambda max_i |H_ii|. c is lm_damp_coeff of zOptDM, and lambda is a
 * dimensionless value that decreases as the iteration converges, which is specified by the damping
 * method (see zOptDMAssignLM()). Since mu is proportional to the magnitude of H, the step is invariant
 * with respect to the scales of the variable and the function. If c is too large, the iteration behaves
 * like the steepest descent method at the beginning, and might go to another local minimum.
 * lm_damp_coeff is initialized by ZM_OPT_DM_LM_DEFAULT_DAMP_COEFF in zOptDMCreate(), and set by
 * zOptDMSetLMDampCoeff().
 */
#define ZM_OPT_DM_LM_DEFAULT_DAMP_COEFF ( 1.0e-2 )

/*! \brief set the coefficient of the damping of Levenberg-Marquardt's method. */
#define zOptDMSetLMDampCoeff(opt,c) ( (opt)->lm_damp_coeff = (c) )

/*! \brief default scale of search range for More-Thuente's method, Dennis-Schnabel's method and Brent's
 * method relative to the magnitude of the initial guess. */
#define ZM_OPT_DM_DEFAULT_SCALE  ( 1.0e3 )
/*! \brief default step length of the line search with fixed step. */
#define ZM_OPT_DM_DEFAULT_STEP     ( 1.0 )

__BEGIN_DECLS

/*! \brief step of numerical differentiation.
 *
 * zOptDMDiffStep() computes the step of numerical differentiation for the \a i th component x_i of
 * the variable \a var, which is \a eps s_i for s_i = max{ |x_i|, \a ratio |x|inf } (s_i = 1 if
 * x = 0). \a xmax is |x|inf, namely, the maximum absolute value of the components of \a var.
 * See ZM_OPT_DM_EPS and ZM_OPT_DM_DIFF_SCALE_RATIO for the choice of \a eps and \a ratio.
 * \return
 * zOptDMDiffStep() returns ( x_i + \a eps s_i ) - x_i, which is the actual step represented in the
 * floating-point arithmetic. A difference quotient should be computed with it in order to cancel the
 * round-off error in representing the perturbed variable.
 */
__ZM_EXPORT double zOptDMDiffStep(const zVec var, int i, double xmax, double eps, double ratio);

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
 * \a scale is the upper bound of the components of a step relative to the maximum absolute component
 * of the initial guess (1 if it is zero) for "MT", "DS" and "Brent", and the step length for "fixed". If it is zero, ZM_OPT_DM_DEFAULT_SCALE is used for "MT", "DS" and "Brent", and
 * ZM_OPT_DM_DEFAULT_STEP is used for "fixed".
 *
 * \a dampmethod for zOptDMAssignLM() specifies the damping coefficient mu = c lambda max_i |H_ii| added
 * to the diagonal components of the Hessian matrix H (see ZM_OPT_DM_LM_DEFAULT_DAMP_COEFF), namely,
 *  "grad" for lambda = |g| / |g0|, where g is the gradient vector and g0 is that at the initial
 *  point, and
 *  "dec" for lambda = min{ |g| / |g0|, d / d0 }, where d = g^T H^-1 g is the Newton decrement that
 *  estimates twice the deviation of the function from the minimum, and d0 is that at the initial
 *  point. d coincides with |F|^2 for the least square problem of F(x)=0 with the Gauss-Newton
 *  approximation of the Hessian matrix if the Jacobian matrix of F is regular. It requires to solve
 *  an additional linear equation at every iteration, and a workspace is allocated internally when it
 *  is assigned at the first time.
 * If H vanishes, mu = |g| is used instead.
 * "grad" is assigned if it is the null pointer or another string.
 *
 * \a updatemethod for zOptDMAssignVM() specifies the update formula of the quasi-Hessian matrix,
 * namely, "DFP" for Davidon-Fletcher-Powell's formula and others for Broyden-Fletcher-Goldfarb-Shanno's.
 * For the latter, the initial identity matrix is scaled by p^T q / q^T q with the first step p and the
 * change of the gradient vector q before the first update, so that the step is invariant with respect
 * to the scales of the variable and the function. It is based on the following paper:
 *  D. F. Shanno and K. H. Phua, Matrix Conditioning and Nonlinear Optimization, Mathematical
 *  Programming, 14(1), pp. 149--160, 1978.
 * See also Eq.(6.20) of J. Nocedal and S. J. Wright, Numerical Optimization, 2nd ed., Springer, 2006.
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
 *  - the maximum absolute component of the step becomes smaller than \a tol ( |\a var|inf + s0 ), where
 *    s0 is the maximum absolute component of the initial guess (1 if it is zero), or
 *  - the number of iterations exceeds \a iter. If \a iter is zero, Z_MAX_ITER_NUM is used instead.
 * The value of the function at the solution is stored where \a eval points, if it is not the null
 * pointer.
 * If the initial guess is at the boundary of the domain of the function, it might fail to compute the
 * numerical gradient. If the minimum is at the boundary of the domain, it might be regarded as a failure
 * even if the solution reaches it.
 * \return
 * zOptDMSolve() returns the number of iterations if the iteration converges. It returns 0 if the gradient
 * vector at the initial point is zero. Otherwise, i.e., if the number of iterations exceeds \a iter, the
 * line search fails, the function value or the gradient vector at the initial point is not finite, or
 * the descent direction is not finite, it returns -1.
 */
__ZM_EXPORT int zOptDMSolve(zOptDM *opt, zVec var, void *util, double tol, int iter, double *eval);

/*! \} */

__END_DECLS

#endif /* __ZM_OPT_DM_H__ */
