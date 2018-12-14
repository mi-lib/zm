#include <zm/zm_ip.h>

static zVec zIPVecTchebychev(zIPData *dat, double t, zVec v);
static zVec zIPVelTchebychev(zIPData *dat, double t, zVec v);
static zVec zIPAccTchebychev(zIPData *dat, double t, zVec v);
static zVec zIPSecVelTchebychev(zIPData *dat, int i, zVec v);
static zVec zIPSecAccTchebychev(zIPData *dat, int i, zVec v);

double _zIPScaleTchebychev(double t, double tmin, double tmax)
{
  return ( 2 * t - ( tmax + tmin ) ) / ( tmax - tmin );
}

/* zIPVecTchebychev
 * - vector on Tchebychev interpolation.
 */
zVec zIPVecTchebychev(zIPData *dat, double t, zVec v)
{
  register int i;
  double x, f1, f2, f;

  zVecClear( v );
  x = _zIPScaleTchebychev( t, zIPTime(dat,0), zIPTime(dat,zIPSize(dat)-1) );
  f1 = 1;
  f2 = x;
  zVecCopy( *zArrayElem(&dat->va,0), v );
  zVecCatNCDRC( v, x, *zArrayElem(&dat->va,1) );
  for( i=2; i<zIPSize(dat); i++ ){
    f = 2 * x * f2 - f1;
    zVecCatDRC( v, f, *zArrayElem(&dat->va,i) );
    f1 = f2;
    f2 = f;
  }
  return v;
}

/* zIPVelTchebychev
 * - velocity on Tchebychev interpolation.
 */
zVec zIPVelTchebychev(zIPData *dat, double t, zVec v)
{
  return NULL;
}

/* zIPAccTchebychev
 * - acceleration on Tchebychev interpolation.
 */
zVec zIPAccTchebychev(zIPData *dat, double t, zVec v)
{
  return NULL;
}

/* zIPSecVelTchebychev
 * - velocity at section on Tchebychev interpolation.
 */
zVec zIPSecVelTchebychev(zIPData *dat, int i, zVec v)
{
  return NULL;
}

/* zIPSecAccTchebychev
 * - acceleration at section on Tchebychev interpolation.
 */
zVec zIPSecAccTchebychev(zIPData *dat, int i, zVec v)
{
  return NULL;
}

/* methods */
static zIPCom _zm_ip_com_tchebychev = {
  zIPVecTchebychev,
  zIPVelTchebychev,
  zIPAccTchebychev,
  zIPSecVelTchebychev,
  zIPSecAccTchebychev,
};

/* zIPCreateTchebychev
 * - create Tchebychev interpolator.
 */
bool zIPCreateTchebychev(zIP *ip, zSeq *seq)
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
    x = _zIPScaleTchebychev( zIPTime(&ip->dat,i), zIPTime(&ip->dat,0), zIPTime(&ip->dat,zIPSize(&ip->dat)-1) );
    zMatElem(r,i,0) = f1 = 1;
    zMatElem(r,i,1) = f2 = x;
    for( j=2; j<zIPSize(&ip->dat); j++ ){
      zMatElem(r,i,j) = f = 2 * x * f2 - f1;
      f1 = f2;
      f2 = f;
    }
    zMatSetRowNC( p, i, cp->data.v );
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
  /* creation of spline interpolator */
  zIPCreateTchebychev( &ip, &seq );

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
