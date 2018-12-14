/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode - ordinary differential equation quadrature.
 */

#ifndef __ZM_ODE_H__
#define __ZM_ODE_H__

#include <zm/zm_vec.h>

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zODE
 * differential equation solver class
 * ********************************************************** */

typedef struct _zODE{
  zVec (* f)(double,zVec,void*,zVec); /* differential function */
  struct _zODE *(* init)(struct _zODE*,int,int,zVec (*)(double,zVec,void*,zVec));
  void (* destroy)(struct _zODE*);
  zVec (* update)(struct _zODE*,double,zVec,double,void*);
  zVec (* cat)(zVec,double,zVec,zVec,void*);
  zVec (* sub)(zVec,zVec,zVec,void*);
  void *_ws; /* workspace for utility */
  zVec _x1, _x2; /* workspace for deferred correction */
} zODE;

/*! \brief concatenate deviation in finite time step to the current variable vector.
 * no need to call this function in user's codes.
 */
__EXPORT zVec _zODECatDefault(zVec x, double dt, zVec dx, zVec xn, void *util);
__EXPORT zVec _zODESubDefault(zVec x1, zVec x2, zVec dx, void *util);

__EXPORT void zODEAssignFunc(zODE *ode, zVec (* catf)(zVec,double,zVec,zVec,void*), zVec (* subf)(zVec,zVec,zVec,void*));

/*! \brief assign differential equation solver.
 *
 * zODEAssign() assignes a numerical differential
 * equation solver to \a ode. \a type is an identifier
 * of the method.
 * Currently, the following methods are available.
 *  'Euler' for Euler method
 *  'Heun' for Heun method
 *  'RK4' for classical Runge-Kutta (Explicit Runge-Kutta,4:4) method
 *  'RKG' for Runge-Kutta-Gill method
 *  'RKF45' for Runge-Kutta-Fehlberg (Embedded Runge-Kutta,4:4-5) method
 *  'Adams' for predictor-corrector method with Adams=Bashforth=Moulton formulae.
 *  'BEuler' for backward Euler method
 *  'TR' for trapezoidal formula method
 *  'BK4' for Butcher-Kuntzmann method (implicit Runge-Kutta,2:4 method)
 *  'Gear' for Gear's method
 * \sa
 * zODEInit, zODEDestroy, zODEUpdate
 */
#define zODEAssign(ode,type,catf,subf) do{\
  (ode)->init = zODEInit_##type;\
  (ode)->destroy = zODEDestroy_##type;\
  (ode)->update = zODEUpdate_##type;\
  zODEAssignFunc( ode, catf, subf );\
} while(0)

/*! \brief initialize, destroy and update differential equation solver.
 *
 * zODEInit() initializes a numerical differential equation solver
 * \a ode. \a dim is the size of vector to be integrated. \a step
 * is a step number of differential vector history to be involved
 * in the computation, only required by 'Adams' and 'Gear'.
 * \a f is a pointer to a function which returns the differential
 * vector at the time \a t and the state \a x by calling \a f(\a t, \a x,\a util, \a dx),
 * where \a dx is a vector to store the resultant. \a util is for
 * programmer s utility to attach any type of data chunk.
 *
 * zODEDestroy() destroys an instance of solver \a ode.
 *
 * zODEUpdate() updates a given vector \a x destructively by \a ode.
 * \a t is the time and \a dt is the quantized integration time step.
 * \a util is for programmer s utility.
 * \return
 * zODEInit() returns a pointer \a ode if succeed to create, or the
 * null pointer otherwise.
 *
 * zODEDestroy() returns no value.
 *
 * zODEUpdate() returns a pointer \a x.
 * \sa
 * zODEAssign
 */
#define zODEInit(ode,dim,step,f) (ode)->init( ode, dim, step, f )
#define zODEDestroy(ode)         (ode)->destroy( ode )
#define zODEUpdate(ode,t,x,dt,u) (ode)->update( ode, t, x, dt, u )

__END_DECLS

#include <zm/zm_ode_dc.h>    /* deferred correction */

#include <zm/zm_ode_euler.h> /* <Euler>: Euler method */
#include <zm/zm_ode_heun.h>  /* <Heun>: Heun method */
#include <zm/zm_ode_rk4.h>   /* <RK4>: classical Runge-Kutta method */
#include <zm/zm_ode_rkg.h>   /* <RKG>: Runge-Kutta-Gill method */
#include <zm/zm_ode_erk.h>   /* <RKF45>: Runge-Kutta-Fehlberg method
                                <CK45>: Cash-Karp method
                                <DP45>: Dormand-Prince method */
#include <zm/zm_ode_adams.h> /* <Adams>: predictor-corrector method */
#include <zm/zm_ode_beuler.h>/* <BEuler>: backward Euler method
                                <TR>: trapezoidal formula method */
#include <zm/zm_ode_bk4.h>   /* <BK4>: Butcher-Kuntzmann method */
#include <zm/zm_ode_gear.h>  /* <Gear>: Gear method */

#include <zm/zm_ode2.h> /* second-order differential equation solver */

#endif /* __ZM_ODE_H__ */
