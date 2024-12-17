/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * interpolation: trapezoidal velocity profiler.
 */

#include <zm/zm_ip.h>

/* create a trapezoidal velocity profiler. */
zTRVelProf *zTRVelProfCreate(zTRVelProf *trvelprof, double v0, double vT, double vmax, double acc, double distance)
{
  double d1, d3, dt3, vm;

  if( distance < 0 ){
    ZRUNERROR( ZM_ERR_IP_TRVELPROF_NEGATIVEDISTANCE, distance );
    return NULL;
  }
  if( vmax <= 0 ){
    ZRUNERROR( ZM_ERR_IP_TRVELPROF_NONPOSITIVEMAXVEL, vmax );
    return NULL;
  }
  if( v0 > vmax ){
    ZRUNERROR( ZM_ERR_IP_TRVELPROF_INVALID_INITIALVEL, v0, vmax );
    return NULL;
  }
  if( vT > vmax ){
    ZRUNERROR( ZM_ERR_IP_TRVELPROF_INVALID_TERMVEL, vT, vmax );
    return NULL;
  }
  if( acc <= 0 ){
    ZRUNERROR( ZM_ERR_IP_TRVELPROF_NONPOSITIVEACC, acc );
    return NULL;
  }
  trvelprof->v0 = v0;
  trvelprof->vT = vT;
  trvelprof->vmax = vmax;
  trvelprof->acc = acc;
  trvelprof->distance = distance;

  trvelprof->_t1 = ( vmax - v0 ) / acc;
  dt3 = ( vmax - vT ) / acc; /* delta time */
  d1 = ( v0 + 0.5 * acc * trvelprof->_t1 ) * trvelprof->_t1;
  d3 = ( vT + 0.5 * acc * dt3 ) * dt3;
  if( d1 + d3 < distance ){
    trvelprof->_t2 = trvelprof->_t1 + ( distance - d1 - d3 ) / vmax;
    trvelprof->_t = trvelprof->_t2 + dt3;
  } else{
    vm = sqrt( 0.5 * ( v0*v0 + vT*vT+ 2*acc*distance ) );
    trvelprof->_t1 = trvelprof->_t2 = ( vm - v0 ) / acc;
    trvelprof->_t = trvelprof->_t1 + ( vm - vT ) / acc;
  }
  return trvelprof;
}

/* travel distance of a trapezoidal velocity profiler. */
double zTRVelProfDist(zTRVelProf *trvelprof, double t)
{
  double d, dt;

  if( t < trvelprof->_t1 ) return ( trvelprof->v0 + 0.5 * trvelprof->acc * t ) * t;
  d = ( trvelprof->v0 + 0.5 * trvelprof->acc * trvelprof->_t1 ) * trvelprof->_t1;
  if( t < trvelprof->_t2 ) return d + trvelprof->vmax * ( t - trvelprof->_t1 );
  dt = trvelprof->_t - t;
  return trvelprof->distance - ( trvelprof->vT + 0.5 * trvelprof->acc * dt ) * dt;
}

/* velocity of a trapezoidal velocity profiler. */
double zTRVelProfVel(zTRVelProf *trvelprof, double t)
{
  if( t < trvelprof->_t1 ) return trvelprof->v0 + trvelprof->acc * t;
  if( t < trvelprof->_t2 ) return trvelprof->vmax;
  return trvelprof->vT + trvelprof->acc * ( trvelprof->_t - t );
}

/* acceleration of a trapezoidal velocity profiler. */
double zTRVelProfAcc(zTRVelProf *trvelprof, double t)
{
  if( t < trvelprof->_t1 ) return trvelprof->acc;
  if( t < trvelprof->_t2 ) return 0;
  return -trvelprof->acc;
}
