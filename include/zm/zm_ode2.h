/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode2 - ordinary differential equation quadrature:
 * second-order ordinary differential equation solver.
 */

#ifndef __ZM_ODE2_H__
#define __ZM_ODE2_H__

/* NOTE: never include this header file in user programs. */

/*! \struct zODE2
 *! \brief second-order ordinary differential equation quadrature class
 */
ZDEF_STRUCT( __ZM_CLASS_EXPORT, zODE2 ){
  zVec (* f)(double,zVec,zVec,void*,zVec); /* differential function */
  zODE2 *(* create)(zODE2*,int,int,zVec (*)(double,zVec,zVec,void*,zVec));
  void (* destroy)(zODE2*);
  zVec (* update)(zODE2*,double,zVec,zVec,double,void*);
  void (* init_history)(zODE2*,double,zVec,zVec,void*); /* only for methods with history */
  zVec (* create_history)(zODE2*,double*,zVec,zVec,double,void*); /* only for methods with history */
  zVec (* cat_dis)(zVec,double,zVec,zVec,void*);
  zVec (* cat_vel)(zVec,double,zVec,zVec,void*);
  zVec (* sub_dis)(zVec,zVec,zVec,void*);
  zVec (* sub_vel)(zVec,zVec,zVec,void*);
  void *util;  /* workspace for utility */
  zODE _ode;   /* need for regular solution */
  zVec _x, _v, _a; /* working memory */
};

__BEGIN_DECLS

/*! \brief initialize a second-order ODE solver. */
#define zODE2Init(ode) do{ \
  (ode)->f = NULL; \
  (ode)->create = NULL; \
  (ode)->destroy = NULL; \
  (ode)->update = NULL; \
  (ode)->init_history = NULL; \
  (ode)->create_history = NULL; \
  (ode)->cat_dis = NULL; \
  (ode)->cat_vel = NULL; \
  (ode)->sub_dis = NULL; \
  (ode)->sub_vel = NULL; \
  (ode)->util = NULL; \
  zODEInit( &(ode)->_ode ); \
  (ode)->_x = (ode)->_v = (ode)->_a = NULL; \
} while(0)

/*! \brief assign concatenation and subtraction functions to a second-order ODE solver. */
__ZM_EXPORT void zODE2AssignFunc(zODE2 *ode, zVec (* catf1)(zVec,double,zVec,zVec,void*), zVec (* catf2)(zVec,double,zVec,zVec,void*), zVec (* subf1)(zVec,zVec,zVec,void*), zVec (* subf2)(zVec,zVec,zVec,void*));

/*! \brief assign a second-order ordinary differential equation quadrature.
 *
 * zODE2Assign() assigns a numerical solver for a second-order ordinary differential equation to \a ode.
 * \a type1 is an identifier of the method. See zODEAssign().
 * \a type2 is an identifier of the particular method for second-order differential equation, where
 *  "Regular", "Symplectic" and "Leapfrog" are available.
 * \sa
 * zODEAssign, zODE2Create, zODE2Destroy, zODE2Update, zODE2InitHistory, zODE2CreateHistory
 */
#define zODE2Assign( d, type, catf1, catf2, subf1, subf2 ) do{\
  (d)->create = zODE2Create##type;\
  (d)->destroy = zODE2Destroy##type;\
  (d)->update = zODE2Update##type;\
  (d)->init_history = zODE2InitHistory##type;\
  (d)->create_history = zODE2CreateHistory##type;\
  zODE2AssignFunc( d, catf1, catf2, subf1, subf2 );\
} while(0)

/*! \brief create, destroy and update a second-order ordinary differential equation quadrature.
 *
 * zODE2Create() creates a numerical solver for a second-order ordinary differential equation \a ode.
 * \a dim is the size of vector to be integrated.
 * \a num is passed to zODECreate() for the first-order ODE solver assigned by zODE2AssignRegular() in
 * "Regular", namely, it is a step number of differential vector history for "Adams" and "Gear", or the
 * maximum number of iterations to solve the nonlinear equation for "BEuler", "TR", "Gauss" and "Radau".
 * See zODECreate(). It is ignored by "Symplectic" and "Leapfrog".
 * \a f is a pointer to a function which returns the second-order differential vector at the time \a t and
 * the state \a x and \a dx by calling \a f(\a t, \a x, \a dx, \a p, \a ddx), where \a ddx is a vector to
 * store the resultant.
 * \a p is for programmer s utility to attach any type of data chunk.
 *
 * zODE2Destroy() destroys an instance \a ode.
 *
 * zODE2Update() updates a given vector \a x and \a dx destructively
 * by \a ode. \a t is the time and \a dt is the quantized integration
 * time step, which has to be nonzero.
 * \return
 * zODE2Create() returns a pointer \a ode if it it succeeds, or the null pointer otherwise.
 *
 * zODE2Destroy() returns no value.
 *
 * zODE2Update() returns a pointer \a x, or the null pointer without updating \a x and \a dx if \a dt is
 * zero.
 * \sa
 * zODE2Assign, zODECreate, zODEDestroy, zODEUpdate
 */
#define zODE2Create(ode,dim,num,f)     ( (ode)->create ? (ode)->create( ode, dim, num, f ) : NULL )
#define zODE2Destroy(ode)              do{ if( (ode)->destroy ) (ode)->destroy( ode ); } while(0)
#define zODE2Update(ode,t,x,v,dt,util) ( (dt) == 0 ? NULL : (ode)->update( ode, t, x, v, dt, util ) )

/*! \brief initialize and create the history of a second-order ordinary differential equation quadrature.
 *
 * zODE2InitHistory() initializes the history internally stored in a second-order ordinary differential
 * equation solver \a ode with the state \a x and \a v at the time \a t, supposing that the system was at
 * rest at the state before \a t. \a util is the same with that for zODE2Update().
 * For "Regular", it calls zODEInitHistory() of the first-order ODE solver assigned by zODE2AssignRegular(),
 * which degrades the order of accuracy of a multistep method only in the first k-1 steps, where k is the
 * number of steps of the method. See zODEInitHistory().
 * For "Leapfrog", it sets the velocity at a half step before \a t for \a v, which is accurate only in the
 * first order.
 *
 * zODE2CreateHistory() creates the history of \a ode with the state \a x and \a v at the time *\a t
 * with a time step \a dt.
 * For "Regular", it calls zODECreateHistory() of the first-order ODE solver assigned by
 * zODE2AssignRegular(), which integrates \a x and \a v from *\a t by k-1 time steps. Namely, *\a t,
 * \a x and \a v are updated to those after k-1 steps. See zODECreateHistory().
 * For "Leapfrog", it computes the velocity at a half step before *\a t from the acceleration at *\a t,
 * which keeps the second-order accuracy. *\a t, \a x and \a v are not modified.
 *
 * They do nothing for "Symplectic", which has no history, and for "Regular" with one-step methods.
 * \return
 * zODE2InitHistory() returns no value.
 *
 * zODE2CreateHistory() returns a pointer \a x, or the null pointer if \a dt is zero or it fails to
 * create a starter.
 * \sa
 * zODEInitHistory, zODECreateHistory
 */
#define zODE2InitHistory(ode,t,x,v,util)        do{ if( (ode)->init_history ) (ode)->init_history( ode, t, x, v, util ); } while(0)
#define zODE2CreateHistory(ode,t,x,v,dt,util)   ( (dt) == 0 ? NULL : (ode)->create_history ? (ode)->create_history( ode, t, x, v, dt, util ) : (x) )

/* regular update */
__ZM_EXPORT zVec _zODE2Cat(zVec y, double dt, zVec dy, zVec yn, void *ode);
__ZM_EXPORT zVec _zODE2Sub(zVec y1, zVec y2, zVec dy, void *ode);

__ZM_EXPORT zODE2 *zODE2CreateRegular(zODE2 *ode, int dim, int num, zVec (* f)(double,zVec,zVec,void*,zVec));
__ZM_EXPORT void zODE2DestroyRegular(zODE2 *ode);
__ZM_EXPORT zVec zODE2UpdateRegular(zODE2* ode, double t, zVec x, zVec v, double dt, void *util);
__ZM_EXPORT void zODE2InitHistoryRegular(zODE2 *ode, double t, zVec x, zVec v, void *util);
__ZM_EXPORT zVec zODE2CreateHistoryRegular(zODE2 *ode, double *t, zVec x, zVec v, double dt, void *util);
#define zODE2AssignRegular( d, type ) zODEAssign( &(d)->_ode, type, _zODE2Cat, _zODE2Sub )

/* simplest symplectic update */
__ZM_EXPORT zODE2 *zODE2CreateSymplectic(zODE2 *ode, int dim, int num, zVec (* f)(double,zVec,zVec,void*,zVec));
__ZM_EXPORT void zODE2DestroySymplectic(zODE2 *ode);
__ZM_EXPORT zVec zODE2UpdateSymplectic(zODE2* ode, double t, zVec x, zVec v, double dt, void *util);
#define zODE2InitHistorySymplectic   NULL /* no history */
#define zODE2CreateHistorySymplectic NULL /* no history */

/* leapfrog update */
__ZM_EXPORT zODE2 *zODE2CreateLeapfrog(zODE2 *ode, int dim, int num, zVec (* f)(double,zVec,zVec,void*,zVec));
__ZM_EXPORT void zODE2DestroyLeapfrog(zODE2 *ode);
__ZM_EXPORT zVec zODE2UpdateLeapfrog(zODE2* ode, double t, zVec x, zVec v, double dt, void *util);
__ZM_EXPORT void zODE2InitHistoryLeapfrog(zODE2 *ode, double t, zVec x, zVec v, void *util);
__ZM_EXPORT zVec zODE2CreateHistoryLeapfrog(zODE2 *ode, double *t, zVec x, zVec v, double dt, void *util);

__END_DECLS

#endif /* __ZM_ODE2_H__ */
