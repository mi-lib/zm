/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_lagrange - interpolation: Lagrange's interpolation.
 */

#include <zm/zm_ip.h>

static zVec zIPVecLagrange(zIPData *dat, double t, zVec v);
static zVec zIPVelLagrange(zIPData *dat, double t, zVec v);
static zVec zIPAccLagrange(zIPData *dat, double t, zVec v);
static zVec zIPSecVelLagrange(zIPData *dat, int i, zVec v);
static zVec zIPSecAccLagrange(zIPData *dat, int i, zVec v);

/* zIPVecLagrange
 * - vector on Lagrange interpolation.
 */
zVec zIPVecLagrange(zIPData *dat, double t, zVec v)
{
  register int i, j;
  double p;

  zVecClear( v );
  for( i=0; i<zIPSize(dat); i++ ){
    for( p=1, j=0; j<zIPSize(dat); j++ )
      if( i != j )
        p *= t - zIPTime(dat,j);
    zVecCatDRC( v, p, *zArrayElem(&dat->va,i) );
  }
  return v;
}

/* zIPVelLagrange
 * - velocity on Lagrange interpolation.
 */
zVec zIPVelLagrange(zIPData *dat, double t, zVec v)
{
  register int i, j, k;
  double p;

  zVecClear( v );
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

/* zIPAccLagrange
 * - acceleration on Lagrange interpolation.
 */
zVec zIPAccLagrange(zIPData *dat, double t, zVec v)
{
  register int i, j, k, l;
  double p;

  zVecClear( v );
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

/* zIPSecVelLagrange
 * - velocity at section on Lagrange interpolation.
 */
zVec zIPSecVelLagrange(zIPData *dat, int i, zVec v)
{
  return zIPVelLagrange( dat, zIPTime(dat,i), v );
}

/* zIPSecAccLagrange
 * - acceleration at section on Lagrange interpolation.
 */
zVec zIPSecAccLagrange(zIPData *dat, int i, zVec v)
{
  return zIPAccLagrange( dat, zIPTime(dat,i), v );
}

/* methods */
static zIPCom _zm_ip_com_lagrange = {
  zIPVecLagrange,
  zIPVelLagrange,
  zIPAccLagrange,
  zIPSecVelLagrange,
  zIPSecAccLagrange,
};

/* zIPCreateLagrange
 * - create Lagrange interpolator.
 */
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
