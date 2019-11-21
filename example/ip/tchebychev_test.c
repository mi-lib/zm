#include <zm/zm_ip.h>

static double _zIPScaleChebyshev(double t, double tmin, double tmax)
{
  return ( 2 * t - ( tmax + tmin ) ) / ( tmax - tmin );
}

/* vector on Chebyshev interpolation. */
static zVec _zIPVecChebyshev(zIPData *dat, double t, zVec v)
{
  register int i;
  double x, f1, f2, f;

  zVecZero( v );
  f1 = 1;
  f2 = x = _zIPScaleChebyshev( t, zIPTime(dat,0), zIPTime(dat,zIPSize(dat)-1) );
  zVecCopyNC( *zArrayElem(&dat->va,0), v );
  zVecCatNCDRC( v, x, *zArrayElem(&dat->va,1) );
  for( i=2; i<zIPSize(dat); i++ ){
    f = 2 * x * f2 - f1;
    zVecCatDRC( v, f, *zArrayElem(&dat->va,i) );
    f1 = f2;
    f2 = f;
  }
  return v;
}

/* velocity on Chebyshev interpolation. */
static zVec _zIPVelChebyshev(zIPData *dat, double t, zVec v)
{
  register int i;
  double x, f1, f2, f;
  double dx, df1, df2, df;

  zVecZero( v );
  f1 = 1;
  f2 = x = _zIPScaleChebyshev( t, zIPTime(dat,0), zIPTime(dat,zIPSize(dat)-1) );
  df1 = 0;
  df2 = dx = 2 / ( zIPTime(dat,zIPSize(dat)-1) - zIPTime(dat,0) );
  zVecMul( *zArrayElem(&dat->va,1), df2, v );
  for( i=2; i<zIPSize(dat); i++ ){
    f = 2 * x * f2 - f1;
    df = 2 * ( dx * f2 + x * df2 ) - df1;
    zVecCatDRC( v, df, *zArrayElem(&dat->va,i) );
    f1 = f2;
    f2 = f;
    df1 = df2;
    df2 = df;
  }
  return v;
}

/* acceleration on Chebyshev interpolation. */
static zVec _zIPAccChebyshev(zIPData *dat, double t, zVec v)
{
  register int i;
  double x, f1, f2, f;
  double dx, df1, df2, df;
  double ddf1, ddf2, ddf;

  zVecZero( v );
  f1 = 1;
  f2 = x = _zIPScaleChebyshev( t, zIPTime(dat,0), zIPTime(dat,zIPSize(dat)-1) );
  df1 = 0;
  df2 = dx = 2 / ( zIPTime(dat,zIPSize(dat)-1) - zIPTime(dat,0) );
  ddf1 = ddf2 = 0;
  for( i=2; i<zIPSize(dat); i++ ){
    f = 2 * x * f2 - f1;
    df = 2 * ( dx * f2 + x * df2 ) - df1;
    ddf = 2 * ( 2 * dx * df2 + x * ddf2 ) - ddf1;
    zVecCatDRC( v, ddf, *zArrayElem(&dat->va,i) );
    f1 = f2;
    f2 = f;
    df1 = df2;
    df2 = df;
    ddf1 = ddf2;
    ddf2 = ddf;
  }
  return v;
}

/* velocity at section on Chebyshev interpolation. */
static zVec _zIPSecVelChebyshev(zIPData *dat, int i, zVec v)
{
  return _zIPVelChebyshev( dat, zIPTime(dat,i), v );
}

/* acceleration at section on Chebyshev interpolation. */
static zVec _zIPSecAccChebyshev(zIPData *dat, int i, zVec v)
{
  return _zIPAccChebyshev( dat, zIPTime(dat,i), v );
}

/* methods */
static zIPCom _zm_ip_com_tchebychev = {
  _zIPVecChebyshev,
  _zIPVelChebyshev,
  _zIPAccChebyshev,
  _zIPSecVelChebyshev,
  _zIPSecAccChebyshev,
};

/* create Chebyshev interpolator. */
bool zIPCreateChebyshev(zIP *ip, zSeq *seq)
{
  register int i, j;
  zSeqListCell *cp;
  zMat r, p, v;
  double x, f1, f2, f;
  bool ret = true;

  if( !zIPDataAlloc( &ip->dat, seq ) ) return false;
  cp = zListHead(seq);
  r = zMatAllocSqr( zIPSize(&ip->dat) );
  p = zMatAlloc( zIPSize(&ip->dat), zVecSizeNC(cp->data.v) );
  v = zMatAlloc( zIPSize(&ip->dat), zVecSizeNC(cp->data.v) );
  if( !r || !p || !v ){
    ZALLOCERROR();
    ret = false;
    goto TERMINATE;
  }
  for( i=0; i<zIPSize(&ip->dat); i++, cp=zListCellPrev(cp) ){
    x = _zIPScaleChebyshev( zIPTime(&ip->dat,i), zIPTime(&ip->dat,0), zIPTime(&ip->dat,zIPSize(&ip->dat)-1) );
    zMatElemNC(r,i,0) = f1 = 1;
    zMatElemNC(r,i,1) = f2 = x;
    for( j=2; j<zIPSize(&ip->dat); j++ ){
      zMatElemNC(r,i,j) = f = 2 * x * f2 - f1;
      f1 = f2;
      f2 = f;
    }
    zMatPutRowNC( p, i, cp->data.v );
  }
  zMulInvMatMat( r, p, v );
  for( i=0; i<zIPSize(&ip->dat); i++ )
    zMatGetRowNC( v, i, *zArrayElem(&ip->dat.va,i) );

  ip->com = &_zm_ip_com_tchebychev;
 TERMINATE:
  zMatFreeAO( 3, r, p, v );
  return ret;
}



int main(int argc, char *argv[])
{
  zSeq seq;
  zIP ip;
  zVec v;
  int point_num;
  double t, tmax;
  register int i;
  /* example data array */
  double tp[] = { 0, 2, 3, 5, 6 };
  double vp1[] = { 0, 3, 1,-2, 4 };
  double vp2[] = { 1, 2,-3, 4, 2 };

  /* creation of x-values and y-values vector */
  point_num = sizeof(tp) / sizeof(double);
  zSeqInit( &seq );
  for( t=0, i=0; i<point_num; i++ ){
    v = zVecCreateList( 2, vp1[i], vp2[i] );
    zSeqEnqueue( &seq, v, tp[i]-t );
    t = tp[i];
  }
  /* creation of Chebyshev interpolator */
  if( !zIPCreateChebyshev( &ip, &seq ) ) return 1;

  /* value, velocity, acceleration */
  tmax = tp[point_num-1];
  v = zVecAlloc( 2 );
  for( i=0; ; i++ ){
    t = tp[0] + 0.1 * i;
    zIPVec( &ip, t, v );
    printf( "%f %f %f ", t, zVecElem(v,0), zVecElem(v,1) );
    zIPVel( &ip, t, v );
    printf( "%f %f ", zVecElem(v,0), zVecElem(v,1) );
    zIPAcc( &ip, t, v );
    printf( "%f %f\n", zVecElem(v,0), zVecElem(v,1) );
    if( t > tmax ) break;
  }

  /* destruction of instances */
  zIPDestroy( &ip );
  zVecFree( v );
  return 0;
}
