/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_ip_trvelprof.h
 * \brief interpolation: trapezoidal velocity profiler.
 * \author Zhidao
 */

#ifndef __ZM_IP_TRVELPROF_H__
#define __ZM_IP_TRVELPROF_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \struct trapezoidal velocity profiler. */
ZDEF_STRUCT( __ZM_CLASS_EXPORT, zTRVelProf ){
  double v0;       /*!< \brief initial velocity */
  double vT;       /*!< \brief terminal velocity */
  double vmax;     /*!< \brief maximum velocity */
  double acc_max;  /*!< \brief constant acceleration */
  /*! \cond */
  double _distance; /* total travel distance */
  double _t;        /* total time */
  double _t1;       /*   0 - _t1: acceleration duration */
  double _t2;       /* _t2 -  _t: deceleration duration */
  /*! \endcond */
#ifdef __cplusplus
  zTRVelProf(double v0, double vT, double vmax, double acc, double dist);
  zTRVelProf();
  double term();
  zTRVelProf *create(double v0, double vT, double vmax, double acc, double dist);
  double distance(double t);
  double velocity(double t);
  double acceleration(double t);
#endif /* __cplusplus */
};

#define zTRVelProfTerm(trvelprof) (trvelprof)->_t

/* create a trapezoidal velocity profiler. */
__ZM_EXPORT zTRVelProf *zTRVelProfCreate(zTRVelProf *trvelprof, double v0, double vT, double vmax, double acc, double dist);
/* travel distance of a trapezoidal velocity profiler. */
__ZM_EXPORT double zTRVelProfDist(const zTRVelProf *trvelprof, double t);
/* velocity of a trapezoidal velocity profiler. */
__ZM_EXPORT double zTRVelProfVel(const zTRVelProf *trvelprof, double t);
/* acceleration of a trapezoidal velocity profiler. */
__ZM_EXPORT double zTRVelProfAcc(const zTRVelProf *trvelprof, double t);

__END_DECLS

#ifdef __cplusplus
inline zTRVelProf::zTRVelProf(double v0, double vT, double vmax, double acc, double dist){ zTRVelProfCreate( this, v0, vT, vmax, acc, dist ); }
inline zTRVelProf::zTRVelProf(){ zTRVelProfCreate( this, 0, 0, 0, 1, 0 ); }
inline double zTRVelProf::term(){ return zTRVelProfTerm( this ); }
inline zTRVelProf *zTRVelProf::create(double v0, double vT, double vmax, double acc, double dist){ return zTRVelProfCreate( this, v0, vT, vmax, acc, dist ); }
inline double zTRVelProf::distance(double t){ return zTRVelProfDist( this, t ); }
inline double zTRVelProf::velocity(double t){ return zTRVelProfVel( this, t ); }
inline double zTRVelProf::acceleration(double t){ return zTRVelProfAcc( this, t ); }
#endif /* __cplusplus */

#endif /* __ZM_IP_TRVELPROF_H__ */
