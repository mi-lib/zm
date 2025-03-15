/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_lin - interpolation: linear interpolation.
 */

#include <zm/zm_ip.h>

/* value on linear interpolation */
static zVec _zIPVecLinear(const zIPData *dat, double t, zVec v)
{
  int i;

  i = zIPSeg( dat, t );
  return zVecCatNC( zIPSecVec(dat,i), t - zIPTime(dat,i), *zArrayElem(&dat->va,i), v );
}

/* velocity at section on linear interpolation */
static zVec _zIPSecVelLinear(const zIPData *dat, int i, zVec v)
{
  return i < 0 || i + 1 >= zIPSize(dat) ?
    zVecZero(v) : zVecCopyNC( *zArrayElem(&dat->va,i), v );
}

/* acceleration at section on linear interpolation */
static zVec _zIPSecAccLinear(const zIPData *dat, int i, zVec v)
{
  int j;

  if( i <= 0 || i + 1 >= zIPSize(dat) ) return zVecZero(v);
  for( j=0; j<zVecSizeNC(v); j++ )
    zVecSetElemNC( v, j,
      zVecArrayElem(&dat->va,i,j) > zVecArrayElem(&dat->va,i-1,j) ? HUGE_VAL :
      ( zVecArrayElem(&dat->va,i,j) < zVecArrayElem(&dat->va,i-1,j) ? -HUGE_VAL : 0 ) );
  return v;
}

/* velocity on linear interpolation */
static zVec _zIPVelLinear(const zIPData *dat, double t, zVec v)
{
  return _zIPSecVelLinear( dat, zIPSeg(dat,t), v );
}

/* acceleration on linear interpolation */
static zVec _zIPAccLinear(const zIPData *dat, double t, zVec v)
{
  int i;

  return zIsTiny( zIPTime( dat, ( i = zIPSeg(dat,t) ) ) - t ) ?
    _zIPSecAccLinear( dat, i, v ) : zVecZero( v );
}

/* methods */
static zIPCom _zm_ip_com_linear = {
  _zIPVecLinear,
  _zIPVelLinear,
  _zIPAccLinear,
  _zIPSecVelLinear,
  _zIPSecAccLinear,
};

/* create linear interpolator */
bool zIPCreateLinear(zIP *ip, const zSeq *seq)
{
  int i;

  if( !zIPDataAlloc( &ip->dat, seq ) ) return false;
  for( i=1; i<zIPSize(&ip->dat); i++ ){
    zVecSubNC( zIPSecVec(&ip->dat,i), zIPSecVec(&ip->dat,i-1), *zArrayElem(&ip->dat.va,i-1) );
    zVecDivDRC( *zArrayElem(&ip->dat.va,i-1), zIPDelta(&ip->dat,i) );
  }
  ip->com = &_zm_ip_com_linear;
  return true;
}
