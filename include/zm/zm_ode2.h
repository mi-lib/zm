/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode2 - ordinary differential equation quadrature:
 * second-order differential equation solver.
 */

#ifndef __ZM_ODE2_H__
#define __ZM_ODE2_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zODE2
 * second-order differential equation solver class
 * ********************************************************** */

typedef struct _zODE2{
  zVec (* f)(double,zVec,zVec,void*,zVec); /* differential function */
  struct _zODE2 *(* init)(struct _zODE2*,int,int,zVec (*)(double,zVec,zVec,void*,zVec));
  void (* destroy)(struct _zODE2*);
  void (* update)(struct _zODE2*,double,zVec,zVec,double,void*);
  zVec (* cat_dis)(zVec,double,zVec,zVec,void*);
  zVec (* cat_vel)(zVec,double,zVec,zVec,void*);
  zVec (* sub_dis)(zVec,zVec,zVec,void*);
  zVec (* sub_vel)(zVec,zVec,zVec,void*);
  void *util;  /* workspace for utility */

  zODE _ode;   /* need for regular solution */
  zVec _x, _v, _a; /* working memory */
} zODE2;

__EXPORT void zODE2AssignFunc(zODE2 *ode, zVec (* catf1)(zVec,double,zVec,zVec,void*), zVec (* catf2)(zVec,double,zVec,zVec,void*), zVec (* subf1)(zVec,zVec,zVec,void*), zVec (* subf2)(zVec,zVec,zVec,void*));

/*! \brief assign differential equation solver.
 *
 * zODE2Assign() assigns a numerical differential
 * equation solver to \a ode. \a type1 is an identifier
 * of the method. See zODEAssign().
 * \a type2 is an identifier of the particular method
 * for second-order differential equation.
 *  'Regular', 'Sympl' and 'Leapfrog' are available.
 * \sa
 * zODEAssign, zODE2Init, zODE2Destroy, zODE2Update
 */
#define zODE2Assign( d, type, catf1, catf2, subf1, subf2 ) do{\
  (d)->init = zODE2Init##type;\
  (d)->destroy = zODE2Destroy##type;\
  (d)->update = zODE2Update##type;\
  zODE2AssignFunc( d, catf1, catf2, subf1, subf2 );\
} while(0)

/*! \brief initialize, destroy and update second-order differential
 * equation solver.
 *
 * zODE2Init() initializes a numerical differential equation
 * solver \a ode. \a dim is the size of vector to be integrated.
 * \a step is a step number of differential vector history to be
 * involved in the computation, only required by 'Adams' and 'Gear'.
 * \a f is a pointer to a function which returns the second-order
 * differential vector at the time \a t and the state \a x and
 * \a dx by calling \a f(\a t, \a x, \a dx, \a p, \a ddx), where
 * \a ddx is a vector to store the resultant.
 * \a p is for programmer s utility to attach any type of data chunk.
 *
 * zODE2Destroy() destroys an instance \a ode.
 *
 * zODE2Update() updates a given vector \a x and \a dx destructively
 * by \a ode. \a t is the time and \a dt is the quantized integration
 * time step.
 * \return
 * zODE2Init() returns a pointer \a ode if it succeeds to create,
 * or the null pointer otherwise.
 *
 * zODE2Destroy() and zODE2Update() return no value.
 * \sa
 * zODE2Assign, zODEInit, zODEDestroy, zODEUpdate
 */
#define zODE2Init(ode,dim,step,f)      (ode)->init( ode, dim, step, f )
#define zODE2Destroy(ode)              (ode)->destroy( ode )
#define zODE2Update(ode,t,x,v,dt,util) (ode)->update( ode, t, x, v, dt, util )

/* regular update */
__EXPORT zVec _zODE2Cat(zVec y, double dt, zVec dy, zVec yn, void *ode);
__EXPORT zVec _zODE2Sub(zVec y1, zVec y2, zVec dy, void *ode);

__EXPORT zODE2 *zODE2InitRegular(zODE2 *ode, int dim, int step, zVec (* f)(double,zVec,zVec,void*,zVec));
__EXPORT void zODE2DestroyRegular(zODE2 *ode);
__EXPORT void zODE2UpdateRegular(zODE2* ode, double t, zVec x, zVec v, double dt, void *util);
#define zODE2AssignRegular( d, type ) zODEAssign( &(d)->_ode, type, _zODE2Cat, _zODE2Sub )

/* simplest symplectic update */
__EXPORT zODE2 *zODE2InitSympl(zODE2 *ode, int dim, int step, zVec (* f)(double,zVec,zVec,void*,zVec));
__EXPORT void zODE2DestroySympl(zODE2 *ode);
__EXPORT void zODE2UpdateSympl(zODE2* ode, double t, zVec x, zVec v, double dt, void *util);

/* leapfrog update */
__EXPORT zODE2 *zODE2InitLeapfrog(zODE2 *ode, int dim, int step, zVec (* f)(double,zVec,zVec,void*,zVec));
__EXPORT void zODE2InitHistLeapfrog(zODE2 *ode, zVec x, zVec v, double dt);
__EXPORT void zODE2DestroyLeapfrog(zODE2 *ode);
__EXPORT void zODE2UpdateLeapfrog(zODE2* ode, double t, zVec x, zVec v, double dt, void *util);

__END_DECLS

#endif /* __ZM_ODE2_H__ */
