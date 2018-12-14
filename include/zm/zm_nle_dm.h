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

__BEGIN_DECLS

/* ********************************************************** */
/*! \defgroup descent method for simultaneous nonlinear equations
 * \{ *//* ************************************************** */

/* ********************************************************** */
/*! \struct zNLE
 * \brief descent method for simultaneous nonlinear equations
 *
 * zNLE is simultatneous nonlinear equation solver class based
 * on descent method.
 *
 * Available methods include
 * 1. steepest descent method,
 * 2. Levenberg-Marquardt's method,
 * 3. variable metric method (quasi Newton's method, secant method),
 * 4. conjugate gradient method (Fletcher-Reeves' method).
 * 5. Newton-Raphson's method, and
 * 6. Broyden's method.
 *//* ******************************************************* */
typedef struct _zNLE{
  zVec (*f)(zVec,zVec,void*);
  zMat (*jac)(zVec,zMat,void*);
  void *util;
  zVec wn;
  zVec we;
  /*! >cond */
  zMat (*_jac)(struct _zNLE*,zVec,zMat,void*);
  int ( *_solve)(struct _zNLE*,zVec,void*,double,int,double*);
  zVec _f;
  zVec _fw;
  zVec _fp;
  zVec _adg, _prg; /* for numerical gradient computation */
  zMat _j;
  zOptDM _opt;
} zNLE;

/*! \brief create simultaneous nonlinear equations solver.
 *
 * zNLECreate() creates a simultaneous nonlinear equation
 * solver \a nle. \a nv is the size of the answer vector.
 * \a ne is the number of equations.
 * \a f is a pointer to the equation function to be solved,
 * '\a f(x)=0'.
 * \return a pointer \a nle if succeed. Otherwise, the
 * null pointer.
 * \sa zNLEAssign
 */
__EXPORT zNLE *zNLECreate(zNLE *nle, int nv, int ne, double scale, zVec (*f)(zVec,zVec,void*), zMat (*jac)(zVec,zMat,void*));

/*! \brief destroy a simultaneous nonlinear equation solver.
 *
 * zNLEDestroy() destroys a solver instance \a nle.
 * \sa zNLECreate
 */
__EXPORT void zNLEDestroy(zNLE *nle);

/*! \brief assign simultaneous nonlinear equation solver.
 *
 * zNLEAssign() assignes a numerical method to solve
 * simultaneous nonlinear equations to \a nle.
 * Currently, the following methods are available.
 *  'SD' for steepest descent method,
 *  'LM' for Levenberg-Marquardt's method,
 *  'VM' for variable metric method,
 *  'CG' for conjugate gradient metric method,
 *  'NR' for Newton-Raphson's method, and
 *  'Broyden' for Broyden's method.
 */
__EXPORT zNLE *zNLEAssignSD(zNLE *nle, const char *stepmethod);
__EXPORT zNLE *zNLEAssignLM(zNLE *nle, const char *stepmethod);
__EXPORT zNLE *zNLEAssignVM(zNLE *nle, const char *stepmethod, const char *updatemethod);
__EXPORT zNLE *zNLEAssignCG(zNLE *nle);
__EXPORT zNLE *zNLEAssignNR(zNLE *nle);
__EXPORT zNLE *zNLEAssignBroyden(zNLE *nle);

/*! \brief solve simultaneous nonlinear equations.
 *
 * zNLESolve() solves simultaneous nonlinear equations.
 * The answer is put into \a var. Note that \a var should
 * be initialized for suitable iteration.
 * \a util is for a programmer's utility to attach any
 * type of data chunk.
 * \a tol is the error tolerance for the termination of
 * iteration.
 * \a iter is the maximum iteration number. Iteration
 * would be done up to \a iter times. When \a iter is
 * zero, Z_MAX_ITER_NUM defined in zm_misc.h is chosen
 * instead.
 * \return the number of iteration, or -1 if something is wrong.
 */
#define zNLESolve(nle,var,util,tol,iter,err) (nle)->_solve( nle, var, util, tol, iter, err )

/*! \} */

__END_DECLS

#endif /* __ZM_NLE_DM_H__ */
