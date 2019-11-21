/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_lagrange - interpolation: Lagrange's interpolation.
 */

#include <zm/zm_ip.h>

/* vector on Lagrange interpolation. */
static zVec _zIPVecLagrange(zIPData *dat, double t, zVec v)
{
  register int i, j;
  double p;

  zVecZero( v );
  for( i=0; i<zIPSize(dat); i++ ){
    for( p=1, j=0; j<zIPSize(dat); j++ )
      if( i != j )
        p *= t - zIPTime(dat,j);
    zVecCatDRC( v, p, *zArrayElem(&dat->va,i) );
  }
  return v;
}

/* velocity on Lagrange interpolation. */
static zVec _zIPVelLagrange(zIPData *dat, double t, zVec v)
{
  register int i, j, k;
  double p;

  zVecZero( v );
  for( i=0; i<zIPSize(dat); i++ )
    for( j=i+1; j<zIPSize(dat); j++ ){
      for( p=1, k=0; k<zIPSize(dat); k++ )
        if( k != i && k != j )
          p *= t - zIPTime(dat,k);
      zVecCatNCDRC( v, p, *zArrayElem(&dat->va,i) );
      zVecCatNCDRC( v, p, *zArrayElem(&dat->va,j) );
    }
  return v;
}

/* acceleration on Lagrange interpolation. */
static zVec _zIPAccLagrange(zIPData *dat, double t, zVec v)
{
  register int i, j, k, l;
  double p;

  zVecZero( v );
  for( i=0; i<zIPSize(dat); i++ )
    for( j=i+1; j<zIPSize(dat); j++ )
      for( k=j+1; k<zIPSize(dat); k++ ){
        for( p=2, l=0; l<zIPSize(dat); l++ )
          if( l != i && l != j && l != k )
            p *= t - zIPTime(dat,l);
        zVecCatNCDRC( v, p, *zArrayElem(&dat->va,i) );
        zVecCatNCDRC( v, p, *zArrayElem(&dat->va,j) );
        zVecCatNCDRC( v, p, *zArrayElem(&dat->va,k) );
      }
  return v;
}

/* velocity at section on Lagrange interpolation. */
static zVec _zIPSecVelLagrange(zIPData *dat, int i, zVec v)
{
  return _zIPVelLagrange( dat, zIPTime(dat,i), v );
}

/* acceleration at section on Lagrange interpolation. */
static zVec _zIPSecAccLagrange(zIPData *dat, int i, zVec v)
{
  return _zIPAccLagrange( dat, zIPTime(dat,i), v );
}

/* methods */
static zIPCom _zm_ip_com_lagrange = {
  _zIPVecLagrange,
  _zIPVelLagrange,
  _zIPAccLagrange,
  _zIPSecVelLagrange,
  _zIPSecAccLagrange,
};

/* create a Lagrange interpolator. */
bool zIPCreateLagrange(zIP *ip, zSeq *seq)
{
  register int i, j;
  double a;

  if( !zIPDataAlloc( &ip->dat, seq ) ) return false;
  for( i=0; i<zIPSize(&ip->dat); i++ ){
    for( a=1, j=0; j<zIPSize(&ip->dat); j++ )
      if( i != j )
        a *= zIPTime(&ip->dat,i) - zIPTime(&ip->dat,j);
     zVecDivNC( zIPSecVec(&ip->dat,i), a, *zArrayElem(&ip->dat.va,i) );
  }
  ip->com = &_zm_ip_com_lagrange;
  return true;
}
