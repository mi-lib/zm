/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_ip.h
 * \brief interpolation.
 * \author Zhidao
 */

#ifndef __ZM_IP_H__
#define __ZM_IP_H__

#include <zm/zm_seq.h>

__BEGIN_DECLS

ZDEF_STRUCT( __ZM_CLASS_EXPORT, zIPKnotCell ){
  double t;     /*!< time stamp */
  zSeqCell *cp; /*!< a pointer to the corresponding sequence cell */
};

zArrayClass( zIPKnotArray, zIPKnotCell );

ZDEF_STRUCT( __ZM_CLASS_EXPORT, zIPData ){
  const zSeq *seq;         /*!< sequence to be interpolated */
  /*! \cond */
  zIPKnotArray knot; /* an array of time stamps and sequence cells */
  zVecArray va;      /* workspace */
  /*! \endcond */
};

#define zIPDataSize(dat)         zListSize((dat)->seq)
#define zIPDataKnot(dat,i)       zArrayElemNC(&(dat)->knot,i)
#define zIPDataTime(dat,i)       ( zIPDataKnot(dat,i)->t )
#define zIPDataDeltaTime(dat,i)  ( zIPDataKnot(dat,i)->cp->data.dt )
#define zIPDataSecVec(dat,i)     ( zIPDataKnot(dat,i)->cp->data.v )
#define zIPDataSecVal(dat,i,j)   zVecElemNC( zIPDataSecVec(dat,i), j )

/* allocate and free the internal workspace of an interpolator. */
__ZM_EXPORT bool zIPDataAlloc(zIPData *dat, const zSeq *seq);
__ZM_EXPORT void zIPDataFree(zIPData *dat);

/* get segment of a sequence to be interpolated.
 *
 * zIPDataSeg() finds the segment which includes the given point corresponding to \a t in an
 * interpolation data \a dat.
 * \return
 * zIPDataSeg() returns an integer value \a i when t_i <= t < t_i+1.
 */
__ZM_EXPORT int zIPDataSeg(const zIPData *dat, double t);

ZDEF_STRUCT( __ZM_CLASS_EXPORT, zIPCom ){
  zVec (*vec)(const zIPData*,double,zVec);
  zVec (*vel)(const zIPData*,double,zVec);
  zVec (*acc)(const zIPData*,double,zVec);
  zVec (*sec_vel)(const zIPData*,int,zVec);
  zVec (*sec_acc)(const zIPData*,int,zVec);
};

ZDEF_STRUCT( __ZM_CLASS_EXPORT, zIP ){
  zIPData dat;
  zIPCom *com;
};

#define zIPSize(ip)        zIPDataSize( &(ip)->dat )
#define zIPKnot(ip,i)      zIPDataKnot( &(ip)->dat, i )
#define zIPTime(ip,i)      zIPDataTime( &(ip)->dat, i )
#define zIPDeltaTime(ip,i) zIPDataDeltaTime( &(ip)->dat, i )
#define zIPSecVec(ip,i)    zIPDataSecVec( &(ip)->dat, i )
#define zIPSecVal(ip,i,j)  zIPDataSecVal( &(ip)->dat, i, j )

/*! \brief interpolation values.
 *
 * \return
 * zIPVal() returns the value at time \a t by an interpolator \a ip.
 *
 * zIPVel() and zIPAcc() return velocity and acceleration, i.e., the first and second order derivatives,
 * respectively, of the interpolator \a ip at time \a t.
 *
 * zIPSecVal() returns the \a i th section value of the interpolator \a ip.
 *
 * zIPSecVel() and zIPSecAcc() return velocity and acceleration at the \a i th section of \a ip.
 * \notes
 * An interpolation method and a series of points to be interpolated have to be assigned to \a ip
 * in advance. If not, anything might happen.
 */
#define zIPVec(ip,t,v)    (ip)->com->vec( &(ip)->dat, (t), (v) )
#define zIPVel(ip,t,v)    (ip)->com->vel( &(ip)->dat, (t), (v) )
#define zIPAcc(ip,t,v)    (ip)->com->acc( &(ip)->dat, (t), (v) )
#define zIPSecVel(ip,i,v) (ip)->com->sec_vel( &(ip)->dat, (i), (v) )
#define zIPSecAcc(ip,i,v) (ip)->com->sec_acc( &(ip)->dat, (i), (v) )

/* destroy an interpolator. */
#define zIPDestroy(ip) zIPDataFree( &(ip)->dat )

__END_DECLS

#include <zm/zm_ip_lin.h>       /* linear interpolation */
#include <zm/zm_ip_lagrange.h>  /* Lagrange's interpolation */
#include <zm/zm_ip_chebyshev.h> /* Chebyshev's interpolation */
#include <zm/zm_ip_spline.h>    /* spline interpolation and Piecewise Cubic Hermite Interporating Polynomial */
#include <zm/zm_ip_akima.h>     /* Akima's interpolation */

#include <zm/zm_ip_pex.h>       /* polynomial curve */
#include <zm/zm_ip_clothoid.h>  /* clothoid curve interpolation */
#include <zm/zm_ip_trvelprof.h> /* trapezoidal velocity profiler */

#endif /* __ZM_IP_H__ */
