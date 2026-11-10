/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode - ordinary differential equation quadrature.
 */

#ifndef __ZM_ODE_H__
#define __ZM_ODE_H__

#include <zm/zm_vec.h>

/*! \struct zODE
 *! \brief ordinary differential equation (ODE) quadrature class
 */
ZDEF_STRUCT( __ZM_CLASS_EXPORT, zODE ){
  zVec (* f)(double,zVec,void*,zVec); /* differential function */
  zODE *(* create)(zODE*,int,int,zVec (*)(double,zVec,void*,zVec));
  void (* destroy)(zODE*);
  zVec (* update)(zODE*,double,zVec,double,void*);
  zVec (* cat)(zVec,double,zVec,zVec,void*);
  zVec (* sub)(zVec,zVec,zVec,void*);
  void (* init_history)(zODE*,double,zVec,void*); /* only for multistep methods */
  zVec (* create_history)(zODE*,double*,zVec,double,void*); /* only for multistep methods */
  int order; /* order of accuracy */
  void *_ws; /* workspace for utility */
  zVec _x1, _x2; /* workspace for deferred correction */
  double _h; /* magnitude of the time step predicted for adaptive step size control, or zero if not yet */
};

__BEGIN_DECLS

/*! \brief initialize an ODE solver. */
#define zODEInit(ode) do{ \
  (ode)->f = NULL; \
  (ode)->create = NULL; \
  (ode)->destroy = NULL; \
  (ode)->update = NULL; \
  (ode)->cat = NULL; \
  (ode)->sub = NULL; \
  (ode)->init_history = NULL; \
  (ode)->create_history = NULL; \
  (ode)->order = 0; \
  (ode)->_ws = NULL; \
  (ode)->_x1 = (ode)->_x2 = NULL; \
  (ode)->_h = 0; \
} while(0)

/*! \brief concatenate deviation in finite time step to the current variable vector.
 * This is for internal implementations of the library. No need to call this function in user's codes.
 */
__ZM_EXPORT zVec _zODECatDefault(zVec x, double dt, zVec dx, zVec xn, void *util);
__ZM_EXPORT zVec _zODESubDefault(zVec x1, zVec x2, zVec dx, void *util);

/*! \brief assign concatenation and subtraction functions to an ODE solver. */
__ZM_EXPORT void zODEAssignFunc(zODE *ode, zVec (* catf)(zVec,double,zVec,zVec,void*), zVec (* subf)(zVec,zVec,zVec,void*));

/*! \brief assign an ordinary differential equation quadrature.
 *
 * zODEAssign() assignes a numerical solver for an ordinary differential equation to \a ode.
 * \a type is an identifier of the method. Currently, the following methods are available.
 *  "Euler" for Euler method
 *  "Heun" for Heun method
 *  "RK4" for classical Runge-Kutta (Explicit Runge-Kutta,4:4) method
 *  "RKG" for Runge-Kutta-Gill method
 *  "RKF45" for Runge-Kutta-Fehlberg (Embedded Runge-Kutta,4:4-5) method
 *  "CK45" for Cash-Karp (Embedded Runge-Kutta,4:4-5) method
 *  "DP45" for Dormand-Prince (Embedded Runge-Kutta,4:4-5) method
 *  "Adams" for predictor-corrector method with Adams=Bashforth=Moulton formulae.
 *  "BEuler" for backward Euler method
 *  "TR" for trapezoidal formula method
 *  "Gauss" for Gauss method (implicit Runge-Kutta,2:4 method)
 *  "Radau" for Radau method (implicit Runge-Kutta,2:4 method)
 *  "Gear" for Gear's method
 * \sa
 * zODECreate, zODEDestroy, zODEUpdate
 *
 * zODEAssign() also sets the order of accuracy of the method for \a ode, which is used by the deferred
 * correction (see zODEUpdateDC()). For "Adams" and "Gear", it is set by zODECreate() since it depends
 * on the step number.
 */
#define zODEAssign(ode,type,catf,subf) do{\
  (ode)->create = zODECreate##type;\
  (ode)->order = zODEOrder##type;\
  (ode)->destroy = zODEDestroy##type;\
  (ode)->update = zODEUpdate##type;\
  (ode)->init_history = zODEInitHistory##type;\
  (ode)->create_history = zODECreateHistory##type;\
  zODEAssignFunc( ode, catf, subf );\
} while(0)

/*! \brief create, destroy and update an ordinary differential equation quadrature.
 *
 * zODECreate() creates a numerical solver for an ordinary differential equation \a ode.
 * \a dim is the size of vector to be integrated.
 * \a num is a step number of differential vector history to be involved in the computation for "Adams"
 * and "Gear". For implicit methods "BEuler", "TR", "Gauss" and "Radau", it is instead the maximum number
 * of iterations to solve the nonlinear equation in each time step, where Z_MAX_ITER_NUM defined in
 * zm_misc.h is chosen if it is zero. It is ignored by the other methods.
 *
 * \a f is a pointer to a function which returns the differential vector at the time \a t and the state \a x
 * by calling \a f(\a t, \a x,\a util, \a dx), where \a dx is a vector to store the resultant.
 * \a util is for programmer s utility to attach any type of data chunk.
 *
 * zODEDestroy() destroys an instance of solver \a ode.
 *
 * zODEUpdate() updates a given vector \a x destructively by \a ode.
 * \a t is the time and \a dt is the quantized integration time step, which has to be nonzero.
 * \a util is for programmer s utility.
 * \return
 * zODECreate() returns a pointer \a ode if it succeeds, or the null pointer otherwise.
 *
 * zODEDestroy() returns no value.
 *
 * zODEUpdate() returns a pointer \a x, or the null pointer without updating \a x if \a dt is zero.
 * \sa
 * zODEAssign
 */
#define zODECreate(ode,dim,num,f)  ( (ode)->create ? (ode)->create( ode, dim, num, f ) : NULL )
#define zODEDestroy(ode)           do{ if( (ode)->destroy ) (ode)->destroy( ode ); } while(0)
#define zODEUpdate(ode,t,x,dt,u)   ( (dt) == 0 ? NULL : (ode)->update( ode, t, x, dt, u ) )

/*! \brief initialize and create the history of a multistep method.
 *
 * The multistep methods "Adams" and "Gear" refer to the history of the past states (for "Gear") or the
 * past differential vectors (for "Adams"), which is filled with zero by zODECreate(). Unless it is
 * initialized properly, the first time steps are not accurate.
 *
 * zODEInitHistory() initializes the history of \a ode as if the state had stayed at \a x before the
 * time \a t (for "Gear"), or as if the differential vector had stayed at that at \a x and \a t (for
 * "Adams"). It is the simplest initialization requiring little computation, while the order of accuracy
 * degrades to one for "Gear" and two for "Adams".
 *
 * zODECreateHistory() creates the history of \a ode by integrating the state \a x forward from the time
 * \a t with a one-step method (starter), where the time step \a dt is the same with that of the
 * following updates. \a x and \a t are updated destructively to the state and the time just after the
 * history is created, namely, they advance by k-1 time steps for a method with k steps. The starter is
 * chosen from the method and the step number so that the order of accuracy is not lower than k-1, which
 * does not degrade the order of accuracy of the method:
 *  "Adams": "RK4" (k<=5) and "DP45" (k>=6),
 *  "Gear":  "Radau" (k<=4), "Gauss" (k=5) and "Gauss" with the deferred correction (k=6),
 * where A-stable implicit methods are chosen for "Gear" not to spoil the stability for stiff equations.
 * Note that "DP45" for "Adams" with more than six steps is not of sufficient order.
 *
 * Both of them do nothing for one-step methods. \a util is for programmer s utility.
 * \return
 * zODEInitHistory() returns no value.
 *
 * zODECreateHistory() returns a pointer \a x, or the null pointer if \a dt is zero or it fails to
 * create the starter. \a x is not updated by a one-step method.
 */
#define zODEInitHistory(ode,t,x,u) do{ \
  if( (ode)->init_history ) (ode)->init_history( ode, t, x, u ); \
} while(0)
#define zODECreateHistory(ode,t,x,dt,u) ( (dt) == 0 ? NULL : (ode)->create_history ? (ode)->create_history( ode, t, x, dt, u ) : (x) )

/*! \brief create the history of a multistep method by a starter.
 *
 * _zODECreateHistory() integrates \a x from the time *\a t by k-1 time steps \a dt with \a starter,
 * where k is the order of accuracy of \a ode, and puts each state before the time step into the history
 * of \a ode by calling \a push(\a ode, t, x, \a util). \a starter has to be assigned in advance, and
 * is created with the deferred correction if \a is_dc is the true value. It is destroyed at the end.
 * \notes
 * This is for internal implementations of the library. No need to call this function in user's codes.
 * \return
 * _zODECreateHistory() returns a pointer \a x, or the null pointer if it fails to create \a starter.
 */
__ZM_EXPORT zVec _zODECreateHistory(zODE *ode, zODE *starter, bool is_dc, void (* push)(zODE*,double,zVec,void*), double *t, zVec x, double dt, void *util);

/*! \brief tolerance for adaptive step size control of ODE quadratures.
 *
 * In ODE quadratures with an adaptive step size control (deferred correction and embedded Runge-Kutta
 * methods), a trial time step is rejected unless every component of two estimations x1 and x2 of the
 * state vector, which start from x0, satisfies
 *  r |x1_i - x2_i| <= ZM_ODE_TOL ( 1 + max{ |x0_i|, |x2_i| } ),
 * namely, the tolerance is absolute for small components and relative for large components.
 * r |x1 - x2| is the estimated error of x2, where r = 1 for the embedded Runge-Kutta methods, and
 * r = 1 / ( 2^p - 1 ) for the deferred correction with a method of order p.
 */
#define ZM_ODE_TOL ( 1.0e-6 )

/*! \brief check if two estimations of the state vector are close enough.
 *
 * zODEStepIsTol() checks if two estimations \a x1 and \a x2 of the state vector starting from \a x0 are close
 * enough in the sense of ZM_ODE_TOL, which is used to decide if a trial time step is accepted.
 * _zODEStepIsTol() is the same with zODEStepIsTol() except that the difference between \a x1 and \a x2 is
 * multiplied by \a r.
 *
 * zODEStepErr() returns the ratio of the estimated error r |x1 - x2| to the tolerance, namely,
 *  max_i r |x1_i - x2_i| / ( ZM_ODE_TOL ( 1 + max{ |x0_i|, |x2_i| } ) ),
 * which is not greater than one if \a x1 and \a x2 are close enough.
 * \return
 * zODEStepIsTol() and _zODEStepIsTol() return the true value if \a x1 and \a x2 are close enough, or the
 * false value otherwise.
 *
 * zODEStepErr() returns the ratio of the estimated error to the tolerance.
 */
__ZM_EXPORT double zODEStepErr(const zVec x0, const zVec x1, const zVec x2, double r);
#define _zODEStepIsTol(x0,x1,x2,r) ( zODEStepErr( x0, x1, x2, r ) <= 1 )
#define zODEStepIsTol(x0,x1,x2)    _zODEStepIsTol( x0, x1, x2, 1.0 )

/*! \brief minimum time step for adaptive step size control of ODE quadratures.
 *
 * In ODE quadratures with an adaptive step size control (deferred correction and embedded Runge-Kutta
 * methods), a trial time step is accepted if it is not greater than this value, and the predicted time
 * step is not smaller than this value.
 */
#define ZM_ODE_DT_TOL ( 1.0e-8 )

/*! \brief parameters to predict the next time step for adaptive step size control of ODE quadratures.
 *
 * In ODE quadratures with an adaptive step size control (deferred correction and embedded Runge-Kutta
 * methods), the next time step is predicted from the current time step h and the ratio err of the estimated error to the tolerance
 * (see zODEStepErr()) as
 *  h min{ ZM_ODE_STEP_RATIO_MAX, max{ ZM_ODE_STEP_RATIO_MIN, ZM_ODE_STEP_SAFETY err^(-1/(q+1)) } },
 * where q is the order of accuracy of the error estimation. ZM_ODE_STEP_RATIO_MAX is replaced by one just
 * after a trial step is rejected.
 */
#define ZM_ODE_STEP_SAFETY    ( 0.9 )
#define ZM_ODE_STEP_RATIO_MIN ( 0.2 )
#define ZM_ODE_STEP_RATIO_MAX ( 5.0 )

/*! \brief update the state vector with an adaptive step size control.
 *
 * _zODEUpdateAdaptive() updates a given vector \a x destructively from the time \a t to \a t + \a dt, where
 * the interval \a dt is divided into time steps, each of which is predicted from the last step and stored
 * in \a ode for the next update.
 * \a trial is a function to try a time step h from the time t and the state x, which returns the ratio of
 * the estimated error to the tolerance (see zODEStepErr()) by calling \a trial(\a ode, t, x, h, \a util,
 * \a data). It must not modify x.
 * \a accept is a function to put the result of the last trial into x by calling \a accept(\a ode, x,
 * \a util, \a data).
 * \a order is the order of accuracy of the error estimation.
 * \a util is for programmer s utility, and \a data is for the method.
 * \notes
 * This is for internal implementations of the library. No need to call this function in user's codes.
 * \return
 * _zODEUpdateAdaptive() returns a pointer \a x. If \a x becomes not a number, the rest of the interval is
 * abandoned.
 */
__ZM_EXPORT zVec _zODEUpdateAdaptive(zODE *ode, double t, zVec x, double dt, int order, double (* trial)(zODE*,double,zVec,double,void*,void*), void (* accept)(zODE*,zVec,void*,void*), void *util, void *data);

__END_DECLS

#include <zm/zm_ode_dc.h>    /* deferred correction */

#include <zm/zm_ode_euler.h> /* "Euler": Euler method */
#include <zm/zm_ode_heun.h>  /* "Heun": Heun method */
#include <zm/zm_ode_rk4.h>   /* "RK4": classical Runge-Kutta method */
#include <zm/zm_ode_rkg.h>   /* "RKG": Runge-Kutta-Gill method */
#include <zm/zm_ode_erk.h>   /* "RKF45": Runge-Kutta-Fehlberg method
                                "CK45": Cash-Karp method
                                "DP45": Dormand-Prince method */
#include <zm/zm_ode_adams.h> /* "Adams": predictor-corrector method */
#include <zm/zm_ode_beuler.h>/* "BEuler": backward Euler method
                                "TR": trapezoidal formula method */
#include <zm/zm_ode_bk4.h>   /* "Gauss": Gauss method
                                "Radau": Radau method */
#include <zm/zm_ode_gear.h>  /* "Gear": Gear method */

#include <zm/zm_ode2.h> /* second-order differential equation solver */

#endif /* __ZM_ODE_H__ */
