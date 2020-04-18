/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_spline - interpolation: spline interpolation.
 */

#include <zm/zm_ip.h>
#include <zm/zm_le.h>

/* vector on spline interpolation */
static zVec _zIPVecSpline(zIPData *dat, double t, zVec v)
{
  register int i;
  double r1, r2;

  i = zIPSeg( dat, t );
  r1 = ( zIPTime(dat,i  ) - t ) / zIPDelta(dat,i+1);
  r2 = ( zIPTime(dat,i+1) - t ) / zIPDelta(dat,i+1);
  zVecZero( v );
  zVecCatNCDRC( v, zSqr(r2)*(3-2*r2), zIPSecVec( dat, i ) );
  zVecCatNCDRC( v, zSqr(r1)*(3+2*r1), zIPSecVec( dat, i+1 ) );
  zVecCatNCDRC( v,-zIPDelta(dat,i+1)*zSqr(r2)*(r2-1), *zArrayElem(&dat->va,i) );
  zVecCatNCDRC( v,-zIPDelta(dat,i+1)*zSqr(r1)*(r1+1), *zArrayElem(&dat->va,i+1) );
  return v;
}

/* velocity on spline interpolation */
static zVec _zIPVelSpline(zIPData *dat, double t, zVec v)
{
  register int i;
  double r1, r2;

  i = zIPSeg( dat, t );
  r1 = ( zIPTime(dat,i  ) - t ) / zIPDelta(dat,i+1);
  r2 = ( zIPTime(dat,i+1) - t ) / zIPDelta(dat,i+1);
  zVecZero( v );
  zVecCatNCDRC( v,-6*r2*(1-r2)/zIPDelta(dat,i+1), zIPSecVec( dat, i ) );
  zVecCatNCDRC( v,-6*r1*(1+r1)/zIPDelta(dat,i+1), zIPSecVec( dat, i+1 ) );
  zVecCatNCDRC( v, r2*(3*r2-2), *zArrayElem(&dat->va,i) );
  zVecCatNCDRC( v, r1*(3*r1+2), *zArrayElem(&dat->va,i+1) );
  return v;
}

/* acceleration on spline interpolation */
static zVec _zIPAccSpline(zIPData *dat, double t, zVec v)
{
  register int i;
  double t1, t2;

  i = zIPSeg( dat, t );
  t1 = t - zIPTime(dat,i);
  t2 = t - zIPTime(dat,i+1);
  zVecZero( v );
  zVecCatNCDRC( v, 6+12*t2/zIPDelta(dat,i+1), zIPSecVec(dat,i) );
  zVecCatNCDRC( v, 2*zIPDelta(dat,i+1)+6*t2, *zArrayElem(&dat->va,i) );
  zVecCatNCDRC( v, 6-12*t1/zIPDelta(dat,i+1), zIPSecVec(dat,i+1) );
  zVecCatNCDRC( v,-2*zIPDelta(dat,i+1)+6*t1, *zArrayElem(&dat->va,i+1) );
  return zVecDivDRC( v, zSqr(zIPDelta(dat,i+1)) );
}

/* velocity at a section on spline interpolation */
static zVec _zIPSecVelSpline(zIPData *dat, int i, zVec v)
{
  return i >= zIPSize(dat) ?
    zVecZero( v ) : zVecCopyNC( *zArrayElem(&dat->va,i), v );
}

/* acceleration at a section on spline interpolation */
static zVec _zIPSecAccSpline(zIPData *dat, int i, zVec v)
{
  zVecZero( v );
  if( i >= zIPSize(dat) ) return v;
  if( i == zIPSize(dat)-1 ){
    zVecCatNCDRC( v, 6.0/zIPDelta(dat,i), zIPSecVec(dat,i-1) );
    zVecCatNCDRC( v,-6.0/zIPDelta(dat,i), zIPSecVec(dat,i) );
    zVecCatNCDRC( v,-4.0, *zArrayElem(&dat->va,i) );
    zVecCatNCDRC( v,-2.0, *zArrayElem(&dat->va,i-1) );
    return zVecDivDRC( v, zIPDelta(dat,i) );
  }
  zVecCatNCDRC( v, 6.0/zIPDelta(dat,i+1), zIPSecVec(dat,i+1) );
  zVecCatNCDRC( v,-6.0/zIPDelta(dat,i+1), zIPSecVec(dat,i) );
  zVecCatNCDRC( v,-4.0, *zArrayElem(&dat->va,i) );
  zVecCatNCDRC( v,-2.0, *zArrayElem(&dat->va,i+1) );
  return zVecDivDRC( v, zIPDelta(dat,i+1) );
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
static void _zIPFixEdgeSpline(zIP *ip, zVec b, zVec d, int i, int j, zVec v)
{
  zVecSetElemNC( b, i, 1 );
  zVecSetElemNC( d, i, zVecElemNC(v,j) );
}

/* set a free edge condition */
static void _zIPFreeEdgeSpline(zIP *ip, zVec a, zVec b, zVec d, int i, int iv, int j)
{
  double r;

  r = 1.0 / zIPDelta(&ip->dat,iv);
  zVecSetElemNC( a, i, 2 * r );
  zVecSetElemNC( b, i,     r );
  zVecSetElemNC( d, i, 3.0*(zIPSecVal(&ip->dat,iv,j)-zIPSecVal(&ip->dat,iv-1,j)) / zSqr( zIPDelta(&ip->dat,iv) ) );
}

/* create a spline interpolator */
bool zIPCreateSpline(zIP *ip, zSeq *seq, int etype1, zVec v1, int etype2, zVec v2)
{
  register int i, j;
  int n;
  zVec a, b, c, d, v;
  double r1, r2;
  bool result = false;

  if( !zIPDataAlloc( &ip->dat, seq ) ) return false;
  n = zIPSize(&ip->dat);
  a = zVecAlloc( n );
  b = zVecAlloc( n );
  c = zVecAlloc( n );
  d = zVecAlloc( n );
  v = zVecAlloc( n );
  if( !a || !b || !c || !d || !v ) goto TERMINATE;

  for( i=0; i<zVecSizeNC(zListHead(seq)->data.v); i++ ){
    /* create tridiagonal equation */
    for( j=1; j<n-1; j++ ){
      r1 = 1.0 / zIPDelta(&ip->dat,j);
      r2 = 1.0 / zIPDelta(&ip->dat,j+1);
      zVecSetElemNC( a, j, r1 );
      zVecSetElemNC( b, j, 2 * ( r1 + r2 ) );
      zVecSetElemNC( c, j, r2 );
      r1 = 3.0 / zSqr( zIPDelta(&ip->dat,j) );
      r2 = 3.0 / zSqr( zIPDelta(&ip->dat,j+1) );
      zVecSetElemNC( d, j,
       -zIPSecVal(&ip->dat,j-1,i)* r1
       +zIPSecVal(&ip->dat,j  ,i)*(r1-r2)
       +zIPSecVal(&ip->dat,j+1,i)* r2 );
    }
    /* setting of the edge type at the beginning point */
    switch( etype1 ){
    case ZSPLINE_FIX_EDGE:  _zIPFixEdgeSpline( ip, b, d, 0, i, v1 );    break;
    case ZSPLINE_FREE_EDGE: _zIPFreeEdgeSpline( ip, b, c, d, 0, 1, i ); break;
    default:
      ZRUNERROR( ZM_ERR_IP_INVTYPE );
      goto TERMINATE;
    }
    /* setting of the edge type at the termination point */
    switch( etype2 ){
    case ZSPLINE_FIX_EDGE:  _zIPFixEdgeSpline( ip, b, d, n-1, i, v2 );      break;
    case ZSPLINE_FREE_EDGE: _zIPFreeEdgeSpline( ip, a, b, d, n-1, n-1, i ); break;
    default:
      ZRUNERROR( ZM_ERR_IP_INVTYPE );
      goto TERMINATE;
    }
    /* solving tridiagonal equation */
    zTridiagSolveDST( a, b, c, d, v );
    for( j=0; j<zIPSize(&ip->dat); j++ )
      zVecArrayElem(&ip->dat.va,j,i) = zVecElemNC(v,j);
  }
  result = true;
  ip->com = &_zm_ip_com_spline;

 TERMINATE:
  zVecFreeAO( 5, a, b, c, d, v );
  return result;
}
