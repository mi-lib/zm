/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip - interpolation.
 */

#ifndef __ZM_IP_H__
#define __ZM_IP_H__

#include <zm/zm_seq.h>

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zIP
 * interpolation curve class
 * ********************************************************** */

typedef struct{
  double t;         /*!< time stamp */
  zSeqListCell *cp; /*!< a pointer to the corresponding sequence cell */
} zIPKnotCell;

zArrayClass( zIPKnotArray, zIPKnotCell );

typedef struct{
  zSeq *seq;         /*!< sequence to be interpolated */
  /*! \cond */
  zIPKnotArray knot; /* an array of time stamps and sequence cells */
  zVecArray va;      /* workspace */
  /*! \endcond */
} zIPData;

#define zIPSize(dat)       zListNum((dat)->seq)
#define zIPKnot(dat,i)     zArrayElem(&(dat)->knot,i)
#define zIPTime(dat,i)     ( zIPKnot(dat,i)->t )
#define zIPDelta(dat,i)    ( zIPKnot(dat,i)->cp->data.dt )
#define zIPSecVec(dat,i)   ( zIPKnot(dat,i)->cp->data.v )
#define zIPSecVal(dat,i,j) zVecElem( zIPSecVec(dat,i), j )

/* zIPDataAlloc, zIPDataFree
 * - allocate and free the internal workspace of an interpolator.
 */
__EXPORT bool zIPDataAlloc(zIPData *dat, zSeq *seq);
__EXPORT void zIPDataFree(zIPData *dat);

/* zIPSeg - get segment of a sequence to be interpolated.
 *
 * zIPSeg() finds the segment which the given point corresponding
 * to \a t belongs to in an interpolation data \a dat.
 * [RETURN VALUE]
 * zIPSeg() returns an integer value \a i when t_i <= t < t_i+1.
 */
__EXPORT int zIPSeg(zIPData *dat, double t);

typedef struct{
  zVec (*vec)(zIPData*,double,zVec);
  zVec (*vel)(zIPData*,double,zVec);
  zVec (*acc)(zIPData*,double,zVec);
  zVec (*sec_vel)(zIPData*,int,zVec);
  zVec (*sec_acc)(zIPData*,int,zVec);
} zIPCom;

typedef struct{
  zIPData dat;
  zIPCom *com;
} zIP;

/* METHOD:
 * zIPVal, zIPVel, zIPAcc, zIPSecVal, zIPSecVel, zIPSecAcc
 * - interpolation values.
 *
 * 'zIPVal()' returns a value interpolated at 'x' by an
 * interpolator 'ip' and a set of control points 'dat'.
 * #
 * 'zIPVel()' and 'zIPAcc()' return velocity and
 * acceleration, namely, the first and second order
 * derivative, respectively, of the interpolation curve
 * at 'x'.
 * #
 * 'zIPSecVal()' returns the 'i' th section value
 * of the interpolation data 'dat' with the interpolator
 * 'ip'.
 * #
 * 'zIPSecVel()' and 'zIPSecAcc()' return velocity
 * and acceleration at the 'i' th section of 'dat'
 * with 'ip'.
 * [RETURNS]
 * All these functions utilize each method assigned
 * in advance. When no interpolator is assigned, it
 * causes segmentation faults.
 */
#define zIPVec(ip,t,v)    (ip)->com->vec( &(ip)->dat, (t), (v) )
#define zIPVel(ip,t,v)    (ip)->com->vel( &(ip)->dat, (t), (v) )
#define zIPAcc(ip,t,v)    (ip)->com->acc( &(ip)->dat, (t), (v) )
#define zIPSecVel(ip,i,v) (ip)->com->sec_vel( &(ip)->dat, (i), (v) )
#define zIPSecAcc(ip,i,v) (ip)->com->sec_acc( &(ip)->dat, (i), (v) )

/* METHOD:
 * zIPDestroy - destroy an interpolator.
 */
#define zIPDestroy(ip) zIPDataFree( &(ip)->dat )

__END_DECLS

#include <zm/zm_ip_lin.h>      /* linear interpolation */
#include <zm/zm_ip_lagrange.h> /* Lagrange's interpolation */
#include <zm/zm_ip_spline.h>   /* spline interpolation */
#include <zm/zm_ip_akima.h>    /* Akima's interpolation */

#include <zm/zm_ip_pex.h>      /* polynomial curve */
#include <zm/zm_ip_ipio.h>     /* interpolation-in-order */

#endif /* __ZM_IP_H__ */
