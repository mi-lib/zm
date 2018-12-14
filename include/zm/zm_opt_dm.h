/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_dm - optimization tools: descent method.
 */

#ifndef __ZM_OPT_DM_H__
#define __ZM_OPT_DM_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* ********************************************************** */
/*! \defgroup descent method for unconstrained nonlinear optimization
 * \{ *//* ************************************************** */

/* ********************************************************** */
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
 *//* ******************************************************* */
typedef struct _zOptDM{
  double (*eval)(zVec,void*);
  zVec (*grad)(zVec,zVec,void*);
  zMat (*hess)(zVec,zMat,void*);
  /*!< cond */
  zVec (*_grad)(struct _zOptDM*,zVec,zVec,void*);
  zMat (*_hess)(struct _zOptDM*,zVec,zMat,void*);
  zVec (*_vec)(struct _zOptDM*,zVec,zVec,void*);
  zVec (*_step)(struct _zOptDM*,zVec,void*,double,double*);
  zMat (*_update_h)(struct _zOptDM*, int count);
  double _scale;
  zVec _x, _d, _g;
  zVec _p, _q, _r;
  zMat _h;
  zIndex _idx;
  double _b;  /* for Fletcher-Reeves */
  double _df; /* decrease factor */
  double _cf; /* curvature factor */
} zOptDM;

__EXPORT zOptDM *zOptDMCreate(zOptDM *opt, int dim, double scale, double (*eval)(zVec,void*), zVec (*grad)(zVec,zVec,void*), zMat (*hess)(zVec,zMat,void*));
__EXPORT zOptDM *zOptDMAssignSD(zOptDM *opt, const char *stepmethod);
__EXPORT zOptDM *zOptDMAssignLM(zOptDM *opt, const char *stepmethod);
__EXPORT zOptDM *zOptDMAssignVM(zOptDM *opt, const char *stepmethod, const char *updatemethod);
__EXPORT zOptDM *zOptDMAssignCG(zOptDM *opt);

__EXPORT void zOptDMDestroy(zOptDM *opt);
__EXPORT int zOptDMSolve(zOptDM *opt, zVec var, void *util, double tol, int iter, double *eval);

/*! \} */

__END_DECLS

#endif /* __ZM_OPT_DM_H__ */
