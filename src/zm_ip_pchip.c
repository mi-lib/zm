/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_pchip - interpolation: Piecewise Cubic Hermite Interporating Polynomial.
 */

#include <zm/zm_ip.h>

/* cubic Hermite base functions and derivatives */

static double _zPCHIPBase1(double t){ return ( 3 - 2*t )*t*t; }
static double _zPCHIPBase2(double t){ return ( t - 1 )*t*t; }
static double _zPCHIPBaseVel1(double t){ return 6*t*( 1 - t ); }
static double _zPCHIPBaseVel2(double t){ return ( 3*t - 2 )*t; }
static double _zPCHIPBaseAcc1(double t){ return 6*( 1 - 2*t ); }
static double _zPCHIPBaseAcc2(double t){ return 2*( 3*t - 1 ); }

/* value on Ferguson curve. */
double zFergusonVal(double t, double term, double x0, double v0, double x1, double v1)
{
  double dt1, dt2;

  dt1 = 1 - ( dt2 = t / term );
  return _zPCHIPBase1(dt1) * x0 + _zPCHIPBase1(dt2) * x1
         + term * ( -_zPCHIPBase2(dt1) * v0 + _zPCHIPBase2(dt2) * v1 );
}

/* vector on PCHIP interpolation */
static zVec _zIPVecPCHIP(const zIPData *dat, double t, zVec v)
{
  int i;
  double dt, dt1, dt2;

  zVecZero( v );
  i = zIPSeg( dat, t );
  dt = zIPDelta(dat,i+1);
  dt1 = ( zIPTime(dat,i+1) - t ) / dt;
  dt2 = ( t - zIPTime(dat,i) ) / dt;
  zVecCatDRC( v, _zPCHIPBase1(dt1), zIPSecVec(dat,i  ) );
  zVecCatDRC( v, _zPCHIPBase1(dt2), zIPSecVec(dat,i+1) );
  zVecCatDRC( v,-dt*_zPCHIPBase2(dt1), *zArrayElem(&dat->va,i  ) );
  zVecCatDRC( v, dt*_zPCHIPBase2(dt2), *zArrayElem(&dat->va,i+1) );
  return v;
}

/* velocity on PCHIP interpolation */
static zVec _zIPVelPCHIP(const zIPData *dat, double t, zVec v)
{
  int i;
  double dt, dt1, dt2;

  zVecZero( v );
  i = zIPSeg( dat, t );
  dt = zIPDelta(dat,i+1);
  dt1 = ( zIPTime(dat,i+1) - t ) / dt;
  dt2 = ( t - zIPTime(dat,i) ) / dt;
  zVecCatDRC( v, _zPCHIPBaseVel1(dt1), zIPSecVec(dat,i  ) );
  zVecCatDRC( v, _zPCHIPBaseVel1(dt2), zIPSecVec(dat,i+1) );
  zVecCatDRC( v,-dt*_zPCHIPBaseVel2(dt1), *zArrayElem(&dat->va,i  ) );
  zVecCatDRC( v, dt*_zPCHIPBaseVel2(dt2), *zArrayElem(&dat->va,i+1) );
  return v;
}

/* acceleration on PCHIP interpolation */
static zVec _zIPAccPCHIP(const zIPData *dat, double t, zVec v)
{
  int i;
  double dt, dt1, dt2;

  zVecZero( v );
  i = zIPSeg( dat, t );
  dt = zIPDelta(dat,i+1);
  dt1 = ( zIPTime(dat,i+1) - t ) / dt;
  dt2 = ( t - zIPTime(dat,i) ) / dt;
  zVecCatDRC( v, _zPCHIPBaseAcc1(dt1), zIPSecVec(dat,i  ) );
  zVecCatDRC( v, _zPCHIPBaseAcc1(dt2), zIPSecVec(dat,i+1) );
  zVecCatDRC( v,-dt*_zPCHIPBaseAcc2(dt1), *zArrayElem(&dat->va,i  ) );
  zVecCatDRC( v, dt*_zPCHIPBaseAcc2(dt2), *zArrayElem(&dat->va,i+1) );
  return v;
}

/* velocity at section on PCHIP interpolation */
static zVec _zIPSecVelPCHIP(const zIPData *dat, int i, zVec v)
{
  return zVecCopy( *zArrayElem(&dat->va,i), v );
}

/* acceleration at section on PCHIP interpolation */
static zVec _zIPSecAccPCHIP(const zIPData *dat, int i, zVec v)
{
  zVecMul( zIPSecVec(dat,i), -6, v );
  return zVecCatDRC( v, -zIPDelta(dat,i+1)*4, *zArrayElem(&dat->va,i) );
}

/* methods */
static zIPCom _zm_ip_com_pchip = {
  _zIPVecPCHIP,
  _zIPVelPCHIP,
  _zIPAccPCHIP,
  _zIPSecVelPCHIP,
  _zIPSecAccPCHIP,
};

/* initialize gradient vectors by three-point cubic interpolation */
static void _zIPInitPCHIP(zIP *ip)
{
  int i, n;
  double dt1, dt2, dt3;

  dt1 = zIPDelta(&ip->dat,1);
  dt2 = zIPDelta(&ip->dat,2);
  dt3 = dt1 + dt2;
  zVecZero( *zArrayElem(&ip->dat.va,0) );
  zVecCatDRC( *zArrayElem(&ip->dat.va,0),-(dt1+dt3)/(dt1*dt3), zIPSecVec(&ip->dat,0) );
  zVecCatDRC( *zArrayElem(&ip->dat.va,0),      dt3 /(dt1*dt2), zIPSecVec(&ip->dat,1) );
  zVecCatDRC( *zArrayElem(&ip->dat.va,0),     -dt1 /(dt2*dt3), zIPSecVec(&ip->dat,2) );
  n = zIPSize(&ip->dat) - 1;
  for( i=1; i<n; i++ ){
    dt1 = zIPDelta(&ip->dat,i  );
    dt2 = zIPDelta(&ip->dat,i+1);
    dt3 = dt1 + dt2;
    zVecZero( *zArrayElem(&ip->dat.va,i) );
    zVecCatDRC( *zArrayElem(&ip->dat.va,i),     -dt2 /(dt1*dt3), zIPSecVec(&ip->dat,i-1) );
    zVecCatDRC( *zArrayElem(&ip->dat.va,i), (dt2-dt1)/(dt1*dt2), zIPSecVec(&ip->dat,i  ) );
    zVecCatDRC( *zArrayElem(&ip->dat.va,i),      dt1 /(dt2*dt3), zIPSecVec(&ip->dat,i+1) );
  }
  dt1 = zIPDelta(&ip->dat,n-1);
  dt2 = zIPDelta(&ip->dat,n  );
  dt3 = dt1 + dt2;
  zVecZero( *zArrayElem(&ip->dat.va,n) );
  zVecCatDRC( *zArrayElem(&ip->dat.va,n),      dt2 /(dt1*dt3), zIPSecVec(&ip->dat,n-2) );
  zVecCatDRC( *zArrayElem(&ip->dat.va,n),     -dt3 /(dt1*dt2), zIPSecVec(&ip->dat,n-2) );
  zVecCatDRC( *zArrayElem(&ip->dat.va,n), (dt2+dt3)/(dt2*dt3), zIPSecVec(&ip->dat,n  ) );
}

/* modify gradient vectors */
static void _zIPModifyPCHIP(zIP *ip)
{
  int i, j, m;
  double d, a, b, l;

  m = zVecSizeNC( zIPSecVec(&ip->dat,0) );
  for( i=0; i<zIPSize(&ip->dat)-1; i++ ){
    for( j=0; j<m; j++ ){
      d = ( zIPSecVal(&ip->dat,i+1,j) - zIPSecVal(&ip->dat,i,j) ) / zIPDelta(&ip->dat,i+1);
      if( zIsTiny( d ) ){
        zVecSetElemNC( *zArrayElem(&ip->dat.va,i  ), j, 0 );
        zVecSetElemNC( *zArrayElem(&ip->dat.va,i+1), j, 0 );
        continue;
      }
      a = zVecElemNC( *zArrayElem(&ip->dat.va,i  ), j ) / d;
      b = zVecElemNC( *zArrayElem(&ip->dat.va,i+1), j ) / d;
      if( ( l = sqrt( a*a + b*b ) ) >= 3 ){
        zVecSetElemNC( *zArrayElem(&ip->dat.va,i  ), j, a * d * 3 / l );
        zVecSetElemNC( *zArrayElem(&ip->dat.va,i+1), j, b * d * 3 / l );
      }
    }
  }
}

/* create a PCHIP interpolator */
bool zIPCreatePCHIP(zIP *ip, const zSeq *seq)
{
  if( !zIPDataAlloc( &ip->dat, seq ) ) return false;
  _zIPInitPCHIP( ip );
  _zIPModifyPCHIP( ip );
  ip->com = &_zm_ip_com_pchip;
  return true;
}
