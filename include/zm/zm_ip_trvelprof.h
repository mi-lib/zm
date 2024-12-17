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
  double acc;      /*!< \brief constant acceleration */
  double distance; /*!< \brief total travel distance */
  /*! \cond */
  double _t;  /* total time */
  double _t1; /*   0 - _t1: acceleration duration */
  double _t2; /* _t2 -  _t: deceleration duration */
  /*! \endcond */
#ifdef __cplusplus
  zTRVelProf(double v0, double vT, double vmax, double acc, double distance);
  zTRVelProf();
  double term();
  zTRVelProf *create(double v0, double vT, double vmax, double acc, double distance);
  double dist(double t);
  double vel(double t);
  double acc(double t);
#endif /* __cplusplus */
};

#define zTRVelProfTerm(trvelprof) (trvelprof)->_t

/* create a trapezoidal velocity profiler. */
zTRVelProf *zTRVelProfCreate(zTRVelProf *trvelprof, double v0, double vT, double vmax, double acc, double distance);
/* travel distance of a trapezoidal velocity profiler. */
double zTRVelProfDist(zTRVelProf *trvelprof, double t);
/* velocity of a trapezoidal velocity profiler. */
double zTRVelProfVel(zTRVelProf *trvelprof, double t);
/* acceleration of a trapezoidal velocity profiler. */
double zTRVelProfAcc(zTRVelProf *trvelprof, double t);

__END_DECLS

#ifdef __cplusplus
inline zTRVelProf::zTRVelProf(double v0, double vT, double vmax, double acc, double distance){ zTRVelProfCreate( this, v0, vT, vmax, acc, distance ); }
inline zTRVelProf::zTRVelProf(){ zTRVelProfCreate( this, 0, 0, 0, 1, 0 ); }
inline zTRVelProf::double term(){ return zTRVelProfTerm( this ); }
inline zTRVelProf *zTRVelProf::create(double v0, double vT, double vmax, double acc, double distance){ return zTRVelProfCreate( this, v0, vT, vmax, acc, distance ); }
inline double zTRVelProf::dist(double t){ return zTRVelProfDist( this, t ); }
inline double zTRVelProf::vel(double t){ return zTRVelProfVel( this, t ); }
inline double zTRVelProf::acc(double t){ return zTRVelProfAcc( this, t ); }
#endif /* __cplusplus */

#endif /* __ZM_IP_TRVELPROF_H__ */
