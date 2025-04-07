/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_lagrange - interpolation: Lagrange's interpolation.
 */

#include <zm/zm_ip.h>

/* vector on Lagrange interpolation */
static zVec _zIPVecLagrange(const zIPData *dat, double t, zVec v)
{
  int i, j;
  double p;

  zVecZero( v );
  for( i=0; i<zIPDataSize(dat); i++ ){
    for( p=1, j=0; j<zIPDataSize(dat); j++ )
      if( i != j )
        p *= t - zIPDataTime(dat,j);
    zVecCatDRC( v, p, *zArrayElem(&dat->va,i) );
  }
  return v;
}

/* velocity on Lagrange interpolation */
static zVec _zIPVelLagrange(const zIPData *dat, double t, zVec v)
{
  int i, j, k;
  double p;

  zVecZero( v );
  for( i=0; i<zIPDataSize(dat); i++ )
    for( j=i+1; j<zIPDataSize(dat); j++ ){
      for( p=1, k=0; k<zIPDataSize(dat); k++ )
        if( k != i && k != j )
          p *= t - zIPDataTime(dat,k);
      zVecCatNCDRC( v, p, *zArrayElem(&dat->va,i) );
      zVecCatNCDRC( v, p, *zArrayElem(&dat->va,j) );
    }
  return v;
}

/* acceleration on Lagrange interpolation */
static zVec _zIPAccLagrange(const zIPData *dat, double t, zVec v)
{
  int i, j, k, l;
  double p;

  zVecZero( v );
  for( i=0; i<zIPDataSize(dat); i++ )
    for( j=i+1; j<zIPDataSize(dat); j++ )
      for( k=j+1; k<zIPDataSize(dat); k++ ){
        for( p=2, l=0; l<zIPDataSize(dat); l++ )
          if( l != i && l != j && l != k )
            p *= t - zIPDataTime(dat,l);
        zVecCatNCDRC( v, p, *zArrayElem(&dat->va,i) );
        zVecCatNCDRC( v, p, *zArrayElem(&dat->va,j) );
        zVecCatNCDRC( v, p, *zArrayElem(&dat->va,k) );
      }
  return v;
}

/* velocity at section on Lagrange interpolation */
static zVec _zIPSecVelLagrange(const zIPData *dat, int i, zVec v)
{
  return _zIPVelLagrange( dat, zIPDataTime(dat,i), v );
}

/* acceleration at section on Lagrange interpolation */
static zVec _zIPSecAccLagrange(const zIPData *dat, int i, zVec v)
{
  return _zIPAccLagrange( dat, zIPDataTime(dat,i), v );
}

/* methods */
static zIPCom _zm_ip_com_lagrange = {
  _zIPVecLagrange,
  _zIPVelLagrange,
  _zIPAccLagrange,
  _zIPSecVelLagrange,
  _zIPSecAccLagrange,
};

/* create a Lagrange interpolator */
bool zIPCreateLagrange(zIP *ip, const zSeq *seq)
{
  int i, j;
  double a;

  if( !zIPDataAlloc( &ip->dat, seq ) ) return false;
  for( i=0; i<zIPSize(ip); i++ ){
    for( a=1, j=0; j<zIPSize(ip); j++ )
      if( i != j )
        a *= zIPTime(ip,i) - zIPTime(ip,j);
     zVecDivNC( zIPSecVec(ip,i), a, *zArrayElem(&ip->dat.va,i) );
  }
  ip->com = &_zm_ip_com_lagrange;
  return true;
}
