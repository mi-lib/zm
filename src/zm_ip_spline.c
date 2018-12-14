/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_spline - interpolation: spline interpolation.
 */

#include <zm/zm_ip.h>
#include <zm/zm_le.h>

static zVec zIPVecSpline(zIPData *dat, double t, zVec v);
static zVec zIPVelSpline(zIPData *dat, double t, zVec v);
static zVec zIPAccSpline(zIPData *dat, double t, zVec v);
static zVec zIPSecVelSpline(zIPData *dat, int i, zVec v);
static zVec zIPSecAccSpline(zIPData *dat, int i, zVec v);

/* zIPVecSpline
 * - vector on spline interpolation.
 */
zVec zIPVecSpline(zIPData *dat, double t, zVec v)
{
  register int i;
  double r1, r2;

  i = zIPSeg( dat, t );
  r1 = ( zIPTime(dat,i  ) - t ) / zIPDelta(dat,i+1);
  r2 = ( zIPTime(dat,i+1) - t ) / zIPDelta(dat,i+1);
  zVecClear( v );
  zVecCatNCDRC( v, zSqr(r2)*(3-2*r2), zIPSecVec( dat, i ) );
  zVecCatNCDRC( v, zSqr(r1)*(3+2*r1), zIPSecVec( dat, i+1 ) );
  zVecCatNCDRC( v,-zIPDelta(dat,i+1)*zSqr(r2)*(r2-1), *zArrayElem(&dat->va,i) );
  zVecCatNCDRC( v,-zIPDelta(dat,i+1)*zSqr(r1)*(r1+1), *zArrayElem(&dat->va,i+1) );
  return v;
}

/* zIPVelSpline
 * - velocity on spline interpolation.
 */
zVec zIPVelSpline(zIPData *dat, double t, zVec v)
{
  register int i;
  double r1, r2;

  i = zIPSeg( dat, t );
  r1 = ( zIPTime(dat,i  ) - t ) / zIPDelta(dat,i+1);
  r2 = ( zIPTime(dat,i+1) - t ) / zIPDelta(dat,i+1);
  zVecClear( v );
  zVecCatNCDRC( v,-6*r2*(1-r2)/zIPDelta(dat,i+1), zIPSecVec( dat, i ) );
  zVecCatNCDRC( v,-6*r1*(1+r1)/zIPDelta(dat,i+1), zIPSecVec( dat, i+1 ) );
  zVecCatNCDRC( v, r2*(3*r2-2), *zArrayElem(&dat->va,i) );
  zVecCatNCDRC( v, r1*(3*r1+2), *zArrayElem(&dat->va,i+1) );
  return v;
}

/* zIPAccSpline
 * - acceleration on spline interpolation.
 */
zVec zIPAccSpline(zIPData *dat, double t, zVec v)
{
  register int i;
  double t1, t2;

  i = zIPSeg( dat, t );
  t1 = t - zIPTime(dat,i);
  t2 = t - zIPTime(dat,i+1);
  zVecClear( v );
  zVecCatNCDRC( v, 6+12*t2/zIPDelta(dat,i+1), zIPSecVec(dat,i) );
  zVecCatNCDRC( v, 2*zIPDelta(dat,i+1)+6*t2, *zArrayElem(&dat->va,i) );
  zVecCatNCDRC( v, 6-12*t1/zIPDelta(dat,i+1), zIPSecVec(dat,i+1) );
  zVecCatNCDRC( v,-2*zIPDelta(dat,i+1)+6*t1, *zArrayElem(&dat->va,i+1) );
  return zVecDivDRC( v, zSqr(zIPDelta(dat,i+1)) );
}

/* zIPSecVelSpline
 * - velocity at section on spline interpolation.
 */
zVec zIPSecVelSpline(zIPData *dat, int i, zVec v)
{
  return i >= zIPSize(dat) ?
    zVecClear( v ) : zVecCopyNC( *zArrayElem(&dat->va,i), v );
}

/* zIPSecAccSpline
 * - acceleration at section on spline interpolation.
 */
zVec zIPSecAccSpline(zIPData *dat, int i, zVec v)
{
  zVecClear( v );
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
  zIPVecSpline,
  zIPVelSpline,
  zIPAccSpline,
  zIPSecVelSpline,
  zIPSecAccSpline,
};

/* zIPCreateSpline
 * - create spline interpolator.
 */
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
      zVecSetElem( a, j, r1 );
      zVecSetElem( b, j, 2 * ( r1 + r2 ) );
      zVecSetElem( c, j, r2 );
      r1 = 3.0 / zSqr( zIPDelta(&ip->dat,j) );
      r2 = 3.0 / zSqr( zIPDelta(&ip->dat,j+1) );
      zVecSetElem( d, j,
       -zIPSecVal(&ip->dat,j-1,i)* r1
       +zIPSecVal(&ip->dat,j  ,i)*(r1-r2)
       +zIPSecVal(&ip->dat,j+1,i)* r2 );
    }
    /* setting of the edge type at the beginning point */
    switch( etype1 ){
    case ZSPLINE_FIX_EDGE:
      zVecSetElem( b, 0, 1 );
      zVecSetElem( d, 0, zVecElem(v1,i) );
      break;
    case ZSPLINE_FREE_EDGE:
      r1 = 1.0 / zIPDelta(&ip->dat,1);
      zVecSetElem( b, 0, 2 * r1 );
      zVecSetElem( c, 0,     r1 );
      zVecSetElem( d, 0, 3.0*(zIPSecVal(&ip->dat,1,i)-zIPSecVal(&ip->dat,0,i))/zSqr( zIPDelta(&ip->dat,1) ) );
      break;
    default:
      ZRUNERROR( ZM_ERR_IP_INVTYPE );
      goto TERMINATE;
    }
    /* setting of the edge type at the termination point */
    switch( etype2 ){
    case ZSPLINE_FIX_EDGE:
      zVecSetElem( b, n-1, 1 );
      zVecSetElem( d, n-1, zVecElem(v2,i) );
      break;
    case ZSPLINE_FREE_EDGE:
      r2 = 1.0 / zIPDelta(&ip->dat,n-1);
      zVecSetElem( a, n-1, 2 * r2 );
      zVecSetElem( b, n-1,     r2 );
      zVecSetElem( d, n-1, 3.0*(zIPSecVal(&ip->dat,n-1,i)-zIPSecVal(&ip->dat,n-2,i))/zSqr(zIPDelta(&ip->dat,n-1)) );
      break;
    default:
      ZRUNERROR( ZM_ERR_IP_INVTYPE );
      goto TERMINATE;
    }
    /* solving tridiagonal equation */
    zTridiagSolveDST( a, b, c, d, v );
    for( j=0; j<zIPSize(&ip->dat); j++ )
      zVecArrayElem(&ip->dat.va,j,i) = zVecElem(v,j);
  }
  result = true;
  ip->com = &_zm_ip_com_spline;

 TERMINATE:
  zVecFreeAO( 5, a, b, c, d, v );
  return result;
}
