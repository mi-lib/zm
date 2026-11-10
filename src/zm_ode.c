/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode - ordinary differential equation quadrature.
 */

#include <zm/zm_ode.h>

/* concatenate deviation in finite time step to the current variable vector. */
zVec _zODECatDefault(zVec x, double dt, zVec dx, zVec xn, void *dummy)
{
  return zVecCatNC( x, dt, dx, xn );
}

/* subtract a variable vector from another. */
zVec _zODESubDefault(zVec x1, zVec x2, zVec dx, void *dummy)
{
  return zVecSubNC( x1, x2, dx );
}

/* assign concatenation and subtraction functions to an ODE solver. */
void zODEAssignFunc(zODE *ode, zVec (* catf)(zVec,double,zVec,zVec,void*), zVec (* subf)(zVec,zVec,zVec,void*))
{
  ode->cat = catf ? catf : _zODECatDefault;
  ode->sub = subf ? subf : _zODESubDefault;
}

/* ratio of the estimated error between two estimations of the state vector to the tolerance. */
double zODEStepErr(const zVec x0, const zVec x1, const zVec x2, double r)
{
  int i;
  double e, err = 0;

  for( i=0; i<zVecSizeNC(x0); i++ ){
    e = r * fabs( zVecElemNC(x1,i) - zVecElemNC(x2,i) ) / ( ZM_ODE_TOL * ( 1 + zMax( fabs( zVecElemNC(x0,i) ), fabs( zVecElemNC(x2,i) ) ) ) );
    if( !( e <= err ) ) err = e; /* NaN is propagated */
  }
  return err;
}

/* ratio of the next time step to the current one predicted from the ratio of the estimated error to the
 * tolerance err, where order is the order of accuracy of the error estimation. */
static double _zODEStepRatio(double err, int order, double ratio_max)
{
  if( err == 0 ) return ratio_max;
  if( !zIsFinite( err ) ) return ZM_ODE_STEP_RATIO_MIN;
  return zLimit( ZM_ODE_STEP_SAFETY * pow( err, -1.0/(order+1) ), ZM_ODE_STEP_RATIO_MIN, ratio_max );
}

/* update the state vector with an adaptive step size control. the interval dt is divided into time
 * steps, each of which is predicted from the last step. */
zVec _zODEUpdateAdaptive(zODE *ode, double t, zVec x, double dt, int order, double (* trial)(zODE*,double,zVec,double,void*,void*), void (* accept)(zODE*,zVec,void*,void*), void *util, void *data)
{
  double rest, h, hs, sgn, err;
  bool is_last, is_rejected = false;

  sgn = dt > 0 ? 1 : -1;
  /* ode->_h is read only here, since trial might update it in a nested adaptive step size control. */
  h = ode->_h > 0 ? ode->_h : fabs( dt );
  for( rest=fabs(dt); rest>0; ){
    /* the last step is shortened to reach the end of the interval exactly */
    hs = ( is_last = ( h >= rest * ( 1 - zTOL ) ) ) ? rest : h;
    err = trial( ode, t, x, sgn*hs, util, data );
    if( err <= 1 || hs <= ZM_ODE_DT_TOL ){ /* accept (the step is not shortened any more if it is small enough) */
      accept( ode, x, util, data );
      t += sgn*hs;
      rest = is_last ? 0 : rest - hs;
      if( !zIsFinite( err ) ){
        /* the state is not a number any more, so that the rest of the interval is abandoned, and the
         * step is predicted from scratch in the next update. */
        ode->_h = 0;
        return x;
      }
      hs *= _zODEStepRatio( err, order, is_rejected ? 1 : ZM_ODE_STEP_RATIO_MAX );
      /* the shortened last step does not shorten the prediction */
      h = is_last ? zMax( h, hs ) : hs;
      is_rejected = false;
    } else{ /* reject */
      h = hs * _zODEStepRatio( err, order, 1 );
      is_rejected = true;
    }
    /* the step is not shorter than ZM_ODE_DT_TOL, so that the interval is completed in finite steps. */
    h = zMax( h, ZM_ODE_DT_TOL );
  }
  ode->_h = h;
  return x;
}

/* create the history of a multistep method by a starter one-step method. */
zVec _zODECreateHistory(zODE *ode, zODE *starter, bool is_dc, void (* push)(zODE*,double,zVec,void*), double *t, zVec x, double dt, void *util)
{
  int i;

  if( !( is_dc ? zODECreateDC( starter, zVecSizeNC(x), 0, ode->f ) : zODECreate( starter, zVecSizeNC(x), 0, ode->f ) ) )
    return NULL;
  for( i=1; i<ode->order; i++, *t+=dt ){
    push( ode, *t, x, util );
    if( is_dc ) zODEUpdateDC( starter, *t, x, dt, util ); else zODEUpdate( starter, *t, x, dt, util );
  }
  if( is_dc ) zODEDestroyDC( starter ); else zODEDestroy( starter );
  return x;
}
