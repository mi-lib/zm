/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_spline - interpolation: spline interpolation.
 */

#include <zm/zm_ip.h>
#include <zm/zm_le.h>

/* vector on spline interpolation */
static zVec _zIPVecSpline(const zIPData *dat, double t, zVec v)
{
  int i;
  double r1, r2, r12, r22;

  i = zIPDataSeg( dat, t );
  r1 = ( zIPDataTime(dat,i  ) - t ) / zIPDataDeltaTime(dat,i+1);
  r2 = ( zIPDataTime(dat,i+1) - t ) / zIPDataDeltaTime(dat,i+1);
  r12 = _zSqr( r1 );
  r22 = _zSqr( r2 );
  zVecZero( v );
  zVecCatNCDRC( v, r22*(3-2*r2), zIPDataSecVec( dat, i ) );
  zVecCatNCDRC( v, r12*(3+2*r1), zIPDataSecVec( dat, i+1 ) );
  zVecCatNCDRC( v,-zIPDataDeltaTime(dat,i+1)*r22*(r2-1), *zArrayElem(&dat->va,i) );
  zVecCatNCDRC( v,-zIPDataDeltaTime(dat,i+1)*r12*(r1+1), *zArrayElem(&dat->va,i+1) );
  return v;
}

/* velocity on spline interpolation */
static zVec _zIPVelSpline(const zIPData *dat, double t, zVec v)
{
  int i;
  double r1, r2;

  i = zIPDataSeg( dat, t );
  r1 = ( zIPDataTime(dat,i  ) - t ) / zIPDataDeltaTime(dat,i+1);
  r2 = ( zIPDataTime(dat,i+1) - t ) / zIPDataDeltaTime(dat,i+1);
  zVecZero( v );
  zVecCatNCDRC( v,-6*r2*(1-r2)/zIPDataDeltaTime(dat,i+1), zIPDataSecVec( dat, i ) );
  zVecCatNCDRC( v,-6*r1*(1+r1)/zIPDataDeltaTime(dat,i+1), zIPDataSecVec( dat, i+1 ) );
  zVecCatNCDRC( v, r2*(3*r2-2), *zArrayElem(&dat->va,i) );
  zVecCatNCDRC( v, r1*(3*r1+2), *zArrayElem(&dat->va,i+1) );
  return v;
}

/* acceleration on spline interpolation */
static zVec _zIPAccSpline(const zIPData *dat, double t, zVec v)
{
  int i;
  double t1, t2;

  i = zIPDataSeg( dat, t );
  t1 = t - zIPDataTime(dat,i);
  t2 = t - zIPDataTime(dat,i+1);
  zVecZero( v );
  zVecCatNCDRC( v, 6+12*t2/zIPDataDeltaTime(dat,i+1), zIPDataSecVec(dat,i) );
  zVecCatNCDRC( v, 2*zIPDataDeltaTime(dat,i+1)+6*t2, *zArrayElem(&dat->va,i) );
  zVecCatNCDRC( v, 6-12*t1/zIPDataDeltaTime(dat,i+1), zIPDataSecVec(dat,i+1) );
  zVecCatNCDRC( v,-2*zIPDataDeltaTime(dat,i+1)+6*t1, *zArrayElem(&dat->va,i+1) );
  return zVecDivDRC( v, zSqr(zIPDataDeltaTime(dat,i+1)) );
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
