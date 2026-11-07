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
  /*!< cond */
  zVec (* _grad)(zOptDM*,const zVec,zVec,void*);
  zMat (* _hess)(zOptDM*,const zVec,zMat,void*);
  zVec (* _vec)(zOptDM*,const zVec,zVec,void*);
  zVec (* _step)(zOptDM*,const zVec,void*,double,double*);
  zMat (* _update_h)(zOptDM*,int);
  double _scale;
  zVec _x, _d, _g;
  zVec _p, _q, _r;
  zMat _h;
  zIndex _idx;
  double _b;  /* for Fletcher-Reeves */
  double _df; /* decrease factor */
  double _cf; /* curvature factor */
};

/*! \brief threshold of the gradient vector for convergence of descent method.
 *
 * zOptDMSolve() regards the iteration converged if the maximum absolute component of the gradient
 * vector becomes smaller than ZM_OPT_DM_GRAD_MIN times that at the initial point.
 * Since a line search judges a decrease of the function value, the gradient vector can be reduced
 * only to the order of the square root of the machine epsilon relatively.
 */
#define ZM_OPT_DM_GRAD_MIN ( 1.0e-8 )

__BEGIN_DECLS

__ZM_EXPORT zOptDM *zOptDMCreate(zOptDM *opt, int dim, double scale, double (* eval)(const zVec,void*), zVec (* grad)(const zVec,zVec,void*), zMat (* hess)(const zVec,zMat,void*));
__ZM_EXPORT zOptDM *zOptDMAssignSD(zOptDM *opt, const char *stepmethod);
__ZM_EXPORT zOptDM *zOptDMAssignLM(zOptDM *opt, const char *stepmethod);
__ZM_EXPORT zOptDM *zOptDMAssignVM(zOptDM *opt, const char *stepmethod, const char *updatemethod);
__ZM_EXPORT zOptDM *zOptDMAssignCG(zOptDM *opt);

__ZM_EXPORT void zOptDMDestroy(zOptDM *opt);
/*! \brief solve an unconstrained nonlinear optimization problem by descent method.
 *
 * zOptDMSolve() minimizes the function assigned to \a opt from the initial guess \a var, and stores
 * the solution into \a var. \a util is a pointer to the utility data for the function.
 * The iteration terminates if one of the following conditions holds:
 *  - the maximum absolute component of the gradient vector becomes smaller than ZM_OPT_DM_GRAD_MIN
 *    times that at the initial point,
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
