/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_spline - interpolation: cubic spline interpolation.
 */

#include <zm/zm_ip.h>
#include <zm/zm_le.h>

/* cubic Hermite base functions and derivatives */

static double _zFergusonBase1(double t){ return ( 3 - 2*t )*t*t; }
static double _zFergusonBase2(double t){ return ( t - 1 )*t*t; }
static double _zFergusonBaseVel1(double t){ return 6*t*( 1 - t ); }
static double _zFergusonBaseVel2(double t){ return ( 3*t - 2 )*t; }
static double _zFergusonBaseAcc1(double t){ return 6*( 1 - 2*t ); }
static double _zFergusonBaseAcc2(double t){ return 2*( 3*t - 1 ); }

/* value on Ferguson curve. */
double zFergusonVal(double t, double term, double x0, double v0, double x1, double v1)
{
  double dt1, dt2;

  dt1 = 1 - ( dt2 = t / term );
  return _zFergusonBase1(dt1) * x0 + _zFergusonBase1(dt2) * x1
         + term * ( -_zFergusonBase2(dt1) * v0 + _zFergusonBase2(dt2) * v1 );
}

/* vector on spline interpolation */
static zVec _zIPVecSpline(const zIPData *dat, double t, zVec v)
{
  int i;
  double dt, dt1, dt2;

  i = zIPDataSeg( dat, t );
  dt = zIPDataDeltaTime(dat,i+1);
  dt1 = ( zIPDataTime(dat,i+1) - t ) / dt;
  dt2 = ( t - zIPDataTime(dat,i) ) / dt;
  zVecZero( v );
  zVecCatDRC( v, _zFergusonBase1(dt1), zIPDataSecVec(dat,i  ) );
  zVecCatDRC( v, _zFergusonBase1(dt2), zIPDataSecVec(dat,i+1) );
  zVecCatDRC( v,-dt*_zFergusonBase2(dt1), *zArrayElem(&dat->va,i  ) );
  zVecCatDRC( v, dt*_zFergusonBase2(dt2), *zArrayElem(&dat->va,i+1) );
  return v;
}

/* velocity on spline interpolation */
static zVec _zIPVelSpline(const zIPData *dat, double t, zVec v)
{
  int i;
  double dt, dt1, dt2;

  i = zIPDataSeg( dat, t );
  dt = zIPDataDeltaTime(dat,i+1);
  dt1 = ( zIPDataTime(dat,i+1) - t ) / dt;
  dt2 = ( t - zIPDataTime(dat,i) ) / dt;
  zVecZero( v );
  zVecCatDRC( v,-_zFergusonBaseVel1(dt1)/dt, zIPDataSecVec(dat,i  ) );
  zVecCatDRC( v, _zFergusonBaseVel1(dt2)/dt, zIPDataSecVec(dat,i+1) );
  zVecCatDRC( v, _zFergusonBaseVel2(dt1), *zArrayElem(&dat->va,i  ) );
  zVecCatDRC( v, _zFergusonBaseVel2(dt2), *zArrayElem(&dat->va,i+1) );
  return v;
}

/* acceleration on spline interpolation */
static zVec _zIPAccSpline(const zIPData *dat, double t, zVec v)
{
  int i;
  double dt, dt1, dt2;

  i = zIPDataSeg( dat, t );
  dt = zIPDataDeltaTime(dat,i+1);
  dt1 = ( zIPDataTime(dat,i+1) - t ) / dt;
  dt2 = ( t - zIPDataTime(dat,i) ) / dt;
  zVecZero( v );
  zVecCatDRC( v, _zFergusonBaseAcc1(dt1)/_zSqr(dt), zIPDataSecVec(dat,i  ) );
  zVecCatDRC( v, _zFergusonBaseAcc1(dt2)/_zSqr(dt), zIPDataSecVec(dat,i+1) );
  zVecCatDRC( v,-_zFergusonBaseAcc2(dt1)/dt, *zArrayElem(&dat->va,i  ) );
  zVecCatDRC( v, _zFergusonBaseAcc2(dt2)/dt, *zArrayElem(&dat->va,i+1) );
  return v;
}

/* velocity at a section on spline interpolation */
static zVec _zIPSecVelSpline(const zIPData *dat, int i, zVec v)
{
  return i >= zIPDataSize(dat) ?
    zVecZero( v ) : zVecCopyNC( *zArrayElem(&dat->va,i), v );
}

/* acceleration at a section on spline interpolation */
static zVec _zIPSecAccSpline(const zIPData *dat, int i, zVec v)
{
  zVecZero( v );
  if( i >= zIPDataSize(dat) ) return v;
  if( i == zIPDataSize(dat)-1 ){
    zVecCatNCDRC( v, 6.0/zIPDataDeltaTime(dat,i), zIPDataSecVec(dat,i-1) );
    zVecCatNCDRC( v,-6.0/zIPDataDeltaTime(dat,i), zIPDataSecVec(dat,i) );
    zVecCatNCDRC( v,-4.0, *zArrayElem(&dat->va,i) );
    zVecCatNCDRC( v,-2.0, *zArrayElem(&dat->va,i-1) );
    return zVecDivDRC( v, zIPDataDeltaTime(dat,i) );
  }
  zVecCatNCDRC( v, 6.0/zIPDataDeltaTime(dat,i+1), zIPDataSecVec(dat,i+1) );
  zVecCatNCDRC( v,-6.0/zIPDataDeltaTime(dat,i+1), zIPDataSecVec(dat,i) );
  zVecCatNCDRC( v,-4.0, *zArrayElem(&dat->va,i) );
  zVecCatNCDRC( v,-2.0, *zArrayElem(&dat->va,i+1) );
  return zVecDivDRC( v, zIPDataDeltaTime(dat,i+1) );
}

/* find coefficients of a cubic polynomial of a segment of a spline interpolation */
static bool _zIPDataSplineCoeff(const zIPData *dat, int i, zVec a, zVec b, zVec c, zVec d)
{
  double dt1, dt2, dt3;

  if( !zVecSizeEqual(*zArrayElem(&dat->va,0),a) ||
      !zVecSizeEqual(*zArrayElem(&dat->va,0),b) ||
      !zVecSizeEqual(*zArrayElem(&dat->va,0),c) ||
      !zVecSizeEqual(*zArrayElem(&dat->va,0),d) ){
    ZRUNERROR( ZM_ERR_IP_SIZEMISMATCH );
    return false;
  }
  /* temporary */
  dt1 = zIPDataDeltaTime(dat,i+1);
  dt2 = _zSqr( dt1 );
  dt3 = dt1 * dt2;
  zVecSub( zIPDataSecVec(dat,i), zIPDataSecVec(dat,i+1), c );
  zVecAdd( *zArrayElem(&dat->va,i), *zArrayElem(&dat->va,i+1), d );
  /* a */
  zVecMulNC( c, 2/dt3, a );
  zVecCatNCDRC( a, 1.0/dt2, d );
  /* b */
  zVecMulNC( c,-3/dt2, b );
  zVecCatNCDRC( b, -1.0/dt1, d );
  zVecCatNCDRC( b, -1.0/dt1, *zArrayElem(&dat->va,i) );
  /* c */
  zVecCopyNC( *zArrayElem(&dat->va,i), c );
  /* d */
  zVecCopy( zIPDataSecVec(dat,i), d );
  return true;
}

/* methods */
static zIPCom _zm_ip_com_spline = {
  _zIPVecSpline,
  _zIPVelSpline,
  _zIPAccSpline,
  _zIPSecVelSpline,
  _zIPSecAccSpline,
};

/* set a fixed edge condition */
static void _zIPFixEdgeSpline(zIP *ip, zVec b, zVec d, int i, int j, const zVec v)
{
  zVecSetElemNC( b, i, 1 );
  zVecSetElemNC( d, i, zVecElemNC(v,j) );
}

/* set a free edge condition */
static void _zIPFreeEdgeSpline(zIP *ip, zVec a, zVec b, zVec d, int i, int iv, int j)
{
  double r;

  r = 1.0 / zIPDeltaTime(ip,iv);
  zVecSetElemNC( a, i, 2 * r );
  zVecSetElemNC( b, i,     r );
  zVecSetElemNC( d, i, 3.0*(zIPSecVal(ip,iv,j)-zIPSecVal(ip,iv-1,j)) / zSqr( zIPDeltaTime(ip,iv) ) );
}

/* create a spline interpolator */
bool zIPCreateSpline(zIP *ip, const zSeq *seq, int etype1, const zVec v1, int etype2, const zVec v2)
{
  int i, j, n;
  zVec a, b, c, d, v;
  double r1, r2;
  bool result = false;

  if( !zIPDataAlloc( &ip->dat, seq ) ) return false;
  n = zIPSize(ip);
  a = zVecAlloc( n );
  b = zVecAlloc( n );
  c = zVecAlloc( n );
  d = zVecAlloc( n );
  v = zVecAlloc( n );
  if( !a || !b || !c || !d || !v ) goto TERMINATE;

  for( i=0; i<zVecSizeNC(zListHead(seq)->data.v); i++ ){
    /* create tridiagonal equation */
    for( j=1; j<n-1; j++ ){
      r1 = 1.0 / zIPDeltaTime(ip,j);
      r2 = 1.0 / zIPDeltaTime(ip,j+1);
      zVecSetElemNC( a, j, r1 );
      zVecSetElemNC( b, j, 2 * ( r1 + r2 ) );
      zVecSetElemNC( c, j, r2 );
      r1 = 3.0 / zSqr( zIPDeltaTime(ip,j) );
      r2 = 3.0 / zSqr( zIPDeltaTime(ip,j+1) );
      zVecSetElemNC( d, j,
       -zIPSecVal(ip,j-1,i)* r1
       +zIPSecVal(ip,j  ,i)*(r1-r2)
       +zIPSecVal(ip,j+1,i)* r2 );
    }
    /* setting of the edge type at the beginning point */
    switch( etype1 ){
    case ZSPLINE_FIX_EDGE:  _zIPFixEdgeSpline( ip, b, d, 0, i, v1 );    break;
    case ZSPLINE_FREE_EDGE: _zIPFreeEdgeSpline( ip, b, c, d, 0, 1, i ); break;
    default:
      ZRUNERROR( ZM_ERR_IP_INVALID_EDGETYPE );
      goto TERMINATE;
    }
    /* setting of the edge type at the termination point */
    switch( etype2 ){
    case ZSPLINE_FIX_EDGE:  _zIPFixEdgeSpline( ip, b, d, n-1, i, v2 );      break;
    case ZSPLINE_FREE_EDGE: _zIPFreeEdgeSpline( ip, a, b, d, n-1, n-1, i ); break;
    default:
      ZRUNERROR( ZM_ERR_IP_INVALID_EDGETYPE );
      goto TERMINATE;
    }
    /* solving tridiagonal equation */
    zLETridiagSolveDST( a, b, c, d, v );
    for( j=0; j<zIPSize(ip); j++ )
      zVecArrayElem(&ip->dat.va,j,i) = zVecElemNC(v,j);
  }
  result = true;
  ip->com = &_zm_ip_com_spline;

 TERMINATE:
  zVecFreeAtOnce( 5, a, b, c, d, v );
  return result;
}

/* find coefficients of a cubic polynomial of a segment of a spline interpolation */
bool zIPSplineCoeff(const zIP *ip, int i, zVec a, zVec b, zVec c, zVec d)
{
  return _zIPDataSplineCoeff( &ip->dat, i, a, b, c, d );
}

/* PCHIP: Piecewise Cubic Hermite Interporating Polynomial interpolator */

/* acceleration at section on PCHIP interpolation */
static zVec _zIPSecAccPCHIP(const zIPData *dat, int i, zVec v)
{
  zVecMul( zIPDataSecVec(dat,i), -6, v );
  return zVecCatDRC( v, -zIPDataDeltaTime(dat,i+1)*4, *zArrayElem(&dat->va,i) );
}

/* methods */
static zIPCom _zm_ip_com_pchip = {
  _zIPVecSpline,
  _zIPVelSpline,
  _zIPAccSpline,
  _zIPSecVelSpline,
  _zIPSecAccPCHIP,
};

/* initialize gradient vectors by three-point cubic interpolation */
static void _zIPInitPCHIP(zIP *ip)
{
  int i, n;
  double dt1, dt2, dt3;

  dt1 = zIPDeltaTime(ip,1);
  dt2 = zIPDeltaTime(ip,2);
  dt3 = dt1 + dt2;
  zVecZero( *zArrayElem(&ip->dat.va,0) );
  zVecCatDRC( *zArrayElem(&ip->dat.va,0),-(dt1+dt3)/(dt1*dt3), zIPSecVec(ip,0) );
  zVecCatDRC( *zArrayElem(&ip->dat.va,0),      dt3 /(dt1*dt2), zIPSecVec(ip,1) );
  zVecCatDRC( *zArrayElem(&ip->dat.va,0),     -dt1 /(dt2*dt3), zIPSecVec(ip,2) );
  n = zIPSize(ip) - 1;
  for( i=1; i<n; i++ ){
    dt1 = zIPDeltaTime(ip,i  );
    dt2 = zIPDeltaTime(ip,i+1);
    dt3 = dt1 + dt2;
    zVecZero( *zArrayElem(&ip->dat.va,i) );
    zVecCatDRC( *zArrayElem(&ip->dat.va,i),     -dt2 /(dt1*dt3), zIPSecVec(ip,i-1) );
    zVecCatDRC( *zArrayElem(&ip->dat.va,i), (dt2-dt1)/(dt1*dt2), zIPSecVec(ip,i  ) );
    zVecCatDRC( *zArrayElem(&ip->dat.va,i),      dt1 /(dt2*dt3), zIPSecVec(ip,i+1) );
  }
  dt1 = zIPDeltaTime(ip,n-1);
  dt2 = zIPDeltaTime(ip,n  );
  dt3 = dt1 + dt2;
  zVecZero( *zArrayElem(&ip->dat.va,n) );
  zVecCatDRC( *zArrayElem(&ip->dat.va,n),      dt2 /(dt1*dt3), zIPSecVec(ip,n-2) );
  zVecCatDRC( *zArrayElem(&ip->dat.va,n),     -dt3 /(dt1*dt2), zIPSecVec(ip,n-2) );
  zVecCatDRC( *zArrayElem(&ip->dat.va,n), (dt2+dt3)/(dt2*dt3), zIPSecVec(ip,n  ) );
}

/* modify gradient vectors */
static void _zIPModifyPCHIP(zIP *ip)
{
  int i, j, m;
  double d, a, b, l;

  m = zVecSizeNC( zIPSecVec(ip,0) );
  for( i=0; i<zIPSize(ip)-1; i++ ){
    for( j=0; j<m; j++ ){
      d = ( zIPSecVal(ip,i+1,j) - zIPSecVal(ip,i,j) ) / zIPDeltaTime(ip,i+1);
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
