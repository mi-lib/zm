/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode2 - ordinary differential equation quadrature:
 * second-order ordinary differential equation solver.
 */

#ifndef __ZM_ODE2_H__
#define __ZM_ODE2_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \struct zODE2
 *! \brief second-order ordinary differential equation quadrature class
 */
ZDEF_STRUCT( __ZM_CLASS_EXPORT, zODE2 ){
  zVec (* f)(double,zVec,zVec,void*,zVec); /* differential function */
  zODE2 *(* create)(zODE2*,int,int,zVec (*)(double,zVec,zVec,void*,zVec));
  void (* destroy)(zODE2*);
  void (* update)(zODE2*,double,zVec,zVec,double,void*);
  zVec (* cat_dis)(zVec,double,zVec,zVec,void*);
  zVec (* cat_vel)(zVec,double,zVec,zVec,void*);
  zVec (* sub_dis)(zVec,zVec,zVec,void*);
  zVec (* sub_vel)(zVec,zVec,zVec,void*);
  void *util;  /* workspace for utility */
  zODE _ode;   /* need for regular solution */
  zVec _x, _v, _a; /* working memory */
};

/*! \brief initialize a second-order ODE solver. */
#define zOD2EInit(ode) do{ \
  (ode)->f = NULL; \
  (ode)->create = NULL; \
  (ode)->destroy = NULL; \
  (ode)->update = NULL; \
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
 *  'Regular', 'Sympl' and 'Leapfrog' are available.
 * \sa
 * zODEAssign, zODE2Create, zODE2Destroy, zODE2Update
 */
#define zODE2Assign( d, type, catf1, catf2, subf1, subf2 ) do{\
  (d)->create = zODE2Create##type;\
  (d)->destroy = zODE2Destroy##type;\
  (d)->update = zODE2Update##type;\
  zODE2AssignFunc( d, catf1, catf2, subf1, subf2 );\
} while(0)

/*! \brief create, destroy and update a second-order ordinary differential equation quadrature.
 *
 * zODE2Create() creates a numerical solver for a second-order ordinary differential equation \a ode.
 * \a dim is the size of vector to be integrated.
 * \a step is a step number of differential vector history to be involved in the computation, only required
 * by 'Adams' and 'Gear'.
 * \a f is a pointer to a function which returns the second-order differential vector at the time \a t and
 * the state \a x and \a dx by calling \a f(\a t, \a x, \a dx, \a p, \a ddx), where \a ddx is a vector to
 * store the resultant.
 * \a p is for programmer s utility to attach any type of data chunk.
 *
 * zODE2Destroy() destroys an instance \a ode.
 *
 * zODE2Update() updates a given vector \a x and \a dx destructively
 * by \a ode. \a t is the time and \a dt is the quantized integration
 * time step.
 * \return
 * zODE2Create() returns a pointer \a ode if it it succeeds, or the null pointer otherwise.
 *
 * zODE2Destroy() and zODE2Update() return no value.
 * \sa
 * zODE2Assign, zODECreate, zODEDestroy, zODEUpdate
 */
#define zODE2Create(ode,dim,step,f)    (ode)->create( ode, dim, step, f )
#define zODE2Destroy(ode)              (ode)->destroy( ode )
#define zODE2Update(ode,t,x,v,dt,util) (ode)->update( ode, t, x, v, dt, util )

/* regular update */
__ZM_EXPORT zVec _zODE2Cat(zVec y, double dt, zVec dy, zVec yn, void *ode);
__ZM_EXPORT zVec _zODE2Sub(zVec y1, zVec y2, zVec dy, void *ode);

__ZM_EXPORT zODE2 *zODE2CreateRegular(zODE2 *ode, int dim, int step, zVec (* f)(double,zVec,zVec,void*,zVec));
__ZM_EXPORT void zODE2DestroyRegular(zODE2 *ode);
__ZM_EXPORT void zODE2UpdateRegular(zODE2* ode, double t, zVec x, zVec v, double dt, void *util);
#define zODE2AssignRegular( d, type ) zODEAssign( &(d)->_ode, type, _zODE2Cat, _zODE2Sub )

/* simplest symplectic update */
__ZM_EXPORT zODE2 *zODE2CreateSymplectic(zODE2 *ode, int dim, int step, zVec (* f)(double,zVec,zVec,void*,zVec));
__ZM_EXPORT void zODE2DestroySymplectic(zODE2 *ode);
__ZM_EXPORT void zODE2UpdateSymplectic(zODE2* ode, double t, zVec x, zVec v, double dt, void *util);

/* leapfrog update */
__ZM_EXPORT zODE2 *zODE2CreateLeapfrog(zODE2 *ode, int dim, int step, zVec (* f)(double,zVec,zVec,void*,zVec));
__ZM_EXPORT void zODE2InitHistoryLeapfrog(zODE2 *ode, zVec x, zVec v, double dt);
__ZM_EXPORT void zODE2DestroyLeapfrog(zODE2 *ode);
__ZM_EXPORT void zODE2UpdateLeapfrog(zODE2* ode, double t, zVec x, zVec v, double dt, void *util);

__END_DECLS

#endif /* __ZM_ODE2_H__ */
