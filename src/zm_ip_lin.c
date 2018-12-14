/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_lin - interpolation: linear interpolation.
 */

#include <zm/zm_ip.h>

static zVec zIPVecLinear(zIPData *dat, double t, zVec v);
static zVec zIPVelLinear(zIPData *dat, double t, zVec v);
static zVec zIPAccLinear(zIPData *dat, double t, zVec v);
static zVec zIPSecVelLinear(zIPData *dat, int i, zVec v);
static zVec zIPSecAccLinear(zIPData *dat, int i, zVec v);

/* zIPValLinear
 * - value on linear interpolation.
 */
zVec zIPVecLinear(zIPData *dat, double t, zVec v)
{
  register int i;

  i = zIPSeg( dat, t );
  return zVecCatNC( zIPSecVec(dat,i), t - zIPTime(dat,i), *zArrayElem(&dat->va,i), v );
}

/* zIPVelLinear
 * - velocity on linear interpolation.
 */
zVec zIPVelLinear(zIPData *dat, double t, zVec v)
{
  return zIPSecVelLinear( dat, zIPSeg(dat,t), v );
}

/* zIPAccLinear
 * - acceleration on linear interpolation.
 */
zVec zIPAccLinear(zIPData *dat, double t, zVec v)
{
  register int i;

  return zIsTiny( zIPTime( dat, ( i = zIPSeg(dat,t) ) ) - t ) ?
    zIPSecAccLinear( dat, i, v ) : zVecClear( v );
}

/* zIPSecVelLinear
 * - velocity at section on linear interpolation.
 */
zVec zIPSecVelLinear(zIPData *dat, int i, zVec v)
{
  return i < 0 || i >= zIPSize(dat)-1 ?
    zVecClear(v) : zVecCopyNC( *zArrayElem(&dat->va,i), v );
}

/* zIPSecAccLinear
 * - acceleration at section on linear interpolation.
 */
zVec zIPSecAccLinear(zIPData *dat, int i, zVec v)
{
  register int j;

  if( i <= 0 || i >= zIPSize(dat)-1 ) return zVecClear(v);
  for( j=0; j<zVecSizeNC(v); j++ )
    zVecSetElem( v, j,
      zVecArrayElem(&dat->va,i,j) > zVecArrayElem(&dat->va,i-1,j) ? HUGE_VAL :
      ( zVecArrayElem(&dat->va,i,j) < zVecArrayElem(&dat->va,i-1,j) ? -HUGE_VAL : 0 ) );
  return v;
}

/* methods */
static zIPCom _zm_ip_com_linear = {
  zIPVecLinear,
  zIPVelLinear,
  zIPAccLinear,
  zIPSecVelLinear,
  zIPSecAccLinear,
};

/* zIPCreateLinear
 * - create linear interpolator (dummy).
 */
bool zIPCreateLinear(zIP *ip, zSeq *seq)
{
  register int i;

  if( !zIPDataAlloc( &ip->dat, seq ) ) return false;
  for( i=1; i<zIPSize(&ip->dat); i++ ){
    zVecSubNC( zIPSecVec(&ip->dat,i), zIPSecVec(&ip->dat,i-1), *zArrayElem(&ip->dat.va,i-1) );
    zVecDivDRC( *zArrayElem(&ip->dat.va,i-1), zIPDelta(&ip->dat,i) );
  }
  ip->com = &_zm_ip_com_linear;
  return true;
}
