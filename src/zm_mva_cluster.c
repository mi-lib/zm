/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_mva_cluster - multivariate analysis analysis: clustering.
 */

#include <zm/zm_mva.h>

/* ********************************************************** */
/* vector cluster class.
 * ********************************************************** */

/* create and initialize a vector cluster. */
zVecCluster *zVecClusterCreate(zVecCluster *c, int core_size)
{
  zListInit( zVecClusterSampleList(c) );
  if( !( c->core = zVecAlloc( core_size ) ) ){
    zVecClusterDestroy( c );
    return NULL;
  }
  c->var = HUGE_VAL;
  c->_sil = NULL;
  return c;
}

/* destroy a vector cluster. */
void zVecClusterDestroy(zVecCluster *c)
{
  zVecAddrListDestroy( zVecClusterSampleList(c) );
  zVecFree( c->core );
  zFree( c->_sil );
}

/* the maximum silhouette in a cluster. */
double zVecClusterMaxSilhouette(const zVecCluster *c)
{
  return c && c->_sil ? c->_sil[zListSize(zVecClusterSampleList(c))-1] : -HUGE_VAL;
}

/* print a vector cluster to a file. */
void zVecClusterFPrint(FILE *fp, const zVecCluster *c)
{
  int i = 0;
  zVecListCell *vc;

  fprintf( fp, "%d members\n", zListSize(zVecClusterSampleList(c)) );
  zListForEach( zVecClusterSampleList(c), vc ){
    fprintf( fp, " %d: ", i++ );
    zVecFPrint( fp, vc->data );
  }
  fprintf( fp, " core: " );
  zVecFPrint( fp, c->core );
}

/* print data of vectors of a vector cluster to a file. */
void zVecClusterValueFPrint(FILE *fp, const zVecCluster *c)
{
  zVecListFPrint( fp, zVecClusterSampleList(c) );
}

/* ********************************************************** */
/* methods for core and error computation
 * ********************************************************** */

/* a default function to compute error of a sample for clustering. */
static zVec _zVecClusterErrorDefault(const zVecClusterMethod *method, const zVec p, const zVec core, void *dummy, zVec err)
{
  return zVecSub( p, core, err );
}

/* a default function to compute the distance between two samples for clustering. */
static double _zVecClusterDistDefault(const zVecClusterMethod *method, const zVec p1, const zVec p2, void *dummy)
{
  return zVecDist( p1, p2 );
}

/* a default function to find core (typically, mean) for clustering. */
static zVec _zVecClusterCoreDefault(const zVecClusterMethod *method, const zVecAddrList *pl, void *dummy, zVec core)
{
  return zVecListMean( pl, core );
}

/* a function to find medoid as the core for clustering. */
static double _zVecClusterDistSum(const zVecClusterMethod *method, const zVecAddrList *pl, zVec v)
{
  zVecAddrListCell *pp;
  double sum = 0;

  zListForEach( pl, pp )
    sum += zVecClusterMethodDistF( method, pp->data, v );
  return sum;
}
static double _zVecClusterDistAve(const zVecClusterMethod *method, const zVecAddrList *pl, zVec v)
{
  return _zVecClusterDistSum( method, pl, v ) / zListSize(pl);
}
static zVec _zVecClusterCoreMedoid(const zVecClusterMethod *method, const zVecAddrList *pl, void *dummy, zVec core)
{
  zVecAddrListCell *cp, *mp;
  double d, dmin = HUGE_VAL;

  mp = zListTail( pl );
  zListForEach( pl, cp ){
    if( ( d = _zVecClusterDistSum( method, pl, cp->data ) ) < dmin ){
      dmin = d;
      mp = cp;
    }
  }
  return zVecCopy( mp->data, core );
}

/* a default loaded mean computation function for clustering. */
static zVec _zVecClusterLoadedMeanDefault(const zVecClusterMethod *method, const zVecAddrList *pl, const double load[], double n, void *dummy, zVec mean)
{
  zVecListCell *vc;
  int i = 0;

  zVecZero( mean );
  zListForEach( pl, vc ){
    zVecCatDRC( mean, load[i], vc->data );
    i++;
  }
  return zVecDivDRC( mean, n );
}

/* a function to compute linear-sum error of a sample for clustering. */
zVec zVecClusterErrorLS(const zVecClusterMethod *method, const zVec p, const zVec core, void *dummy, zVec err)
{
  double xm, e;

  xm = zVecElem(p,zVecSizeNC(p)-1);
  zVecSetElem( p, zVecSizeNC(p)-1, 1 );
  e = zVecInnerProd(core,p) - xm;
  zVecSetElem( p, zVecSizeNC(p)-1, xm );
  zVecSetElem( err, 0, e );
  return err;
}

/* a function to find core of linear-sum for clustering. */
zVec zVecClusterCoreLS(const zVecClusterMethod *method, const zVecAddrList *pl, void *dummy, zVec core)
{
  zVecAddrListCell *vc;
  zMat c;
  zVec b;
  double xm;

  c = zMatAllocSqr( zVecSizeNC(zListTail(pl)->data) );
  b = zVecAlloc( zVecSizeNC(zListTail(pl)->data) );
  if( !c || !b ){
    ZALLOCERROR();
    core = NULL;
    goto TERMINATE;
  }
  zListForEach( pl, vc ){
    xm = zVecElem(vc->data,zVecSizeNC(vc->data)-1);
    zVecSetElem( vc->data, zVecSizeNC(vc->data)-1, 1 );
    zMatAddDyadNC( c, vc->data, vc->data );
    zVecCatNCDRC( b, xm, vc->data );
    zVecSetElem( vc->data, zVecSizeNC(vc->data)-1, xm );
  }
  zLESolveGauss( c, b, core );
 TERMINATE:
  zMatFree( c );
  zVecFree( b );
  return core;
}

/* a loaded-mean function of linear-sum for clustering. */
zVec zVecClusterLoadedMeanLS(const zVecClusterMethod *method, const zVecAddrList *pl, const double load[], double n, void *dummy, zVec mean)
{
  zVecAddrListCell *vc;
  zMat c;
  zVec b, g;
  double xm;
  int i = 0;

  c = zMatAllocSqr( zVecSizeNC(zListTail(pl)->data) );
  b = zVecAlloc( zVecSizeNC(zListTail(pl)->data) );
  g = zVecAlloc( zVecSizeNC(zListTail(pl)->data) );
  if( !c || !b || !g ){
    ZALLOCERROR();
    mean = NULL;
    goto TERMINATE;
  }
  zListForEach( pl, vc ){
    xm = zVecElem(vc->data,zVecSizeNC(vc->data)-1);
    zVecSetElem( vc->data, zVecSizeNC(vc->data)-1, 1 );
    zVecMulNC( vc->data, load[i], g );
    zMatAddDyadNC( c, g, vc->data );
    zVecCatNCDRC( b, load[i]*xm, vc->data );
    zVecSetElem( vc->data, zVecSizeNC(vc->data)-1, xm );
    i++;
  }
  zLESolveGauss( c, b, mean );
 TERMINATE:
  zMatFree( c );
  zVecFree( b );
  zVecFree( g );
  return mean;
}

/* initialize methods for clustering. */
void zVecClusterMethodInit(zVecClusterMethod *method)
{
  /* core of a cluster */
  method->core_fp = NULL;
  method->core_util = NULL;
  method->core_size = 0;
  /* error between samples */
  method->error_fp = NULL;
  method->error = NULL;
  method->error_util = NULL;
  /* distance between samples */
  method->dist_fp = NULL;
  method->dist_util = NULL;
  /* loaded mean of a cluster */
  method->lm_fp = NULL;
  method->lm_util = NULL;
  method->_lerr = NULL;
}

/* set a function to find core of a cluster. */
zVecClusterMethod *zVecClusterMethodSetCoreFunc(zVecClusterMethod *method, int core_size, zVec (* fp)(const zVecClusterMethod*,const zVecAddrList*,void*,zVec), void *util)
{
  if( core_size <= 0 ){
    ZRUNERROR( ZM_ERR_MVA_INVALID_SIZE, core_size );
    return NULL;
  }
  method->core_size = core_size;
  method->core_fp = fp ? fp : _zVecClusterCoreDefault;
  method->core_util = util;
  return method;
}

/* set a function to compute error of a sample in a cluster. */
zVecClusterMethod *zVecClusterMethodSetErrorFunc(zVecClusterMethod *method, int error_size, zVec (* fp)(const zVecClusterMethod*,const zVec,const zVec,void*,zVec), void *util)
{
  if( method->error ) zVecFree( method->error );
  if( error_size <= 0 ){
    ZRUNERROR( ZM_ERR_MVA_INVALID_SIZE, error_size );
    return NULL;
  }
  method->error_fp = fp ? fp : _zVecClusterErrorDefault;
  method->error_util = util;
  return ( method->error = zVecAlloc( error_size ) ) ? method : NULL;
}

/* set a function to compute the distance between two samples in a cluster. */
zVecClusterMethod *zVecClusterMethodSetDistFunc(zVecClusterMethod *method, double (* fp)(const zVecClusterMethod*,const zVec,const zVec,void*), void *util)
{
  method->dist_fp = fp ? fp : _zVecClusterDistDefault;
  method->dist_util = util;
  return method;
}

/* set a function to compute the loaded mean of samples (for GMM). */
zVecClusterMethod *zVecClusterMethodSetLoadedMeanFunc(zVecClusterMethod *method, zVec (* fp)(const zVecClusterMethod*,const zVecAddrList*,const double[],double,void*,zVec), void *util)
{
  if( !method->core_fp ){
    ZRUNERROR( ZM_ERR_MVA_NOCOREFUNC );
    return NULL;
  }
  if( !method->error_fp ){
    ZRUNERROR( ZM_ERR_MVA_NOERRORFUNC );
    return NULL;
  }
  if( !method->error_fp ){
    ZRUNERROR( ZM_ERR_MVA_NODISTFUNC );
    return NULL;
  }
  if( method->_lerr ) zVecFree( method->_lerr );
  method->lm_fp = fp ? fp : _zVecClusterLoadedMeanDefault;
  method->lm_util = util;
  return ( method->_lerr = zVecAlloc( zVecSize(method->error) ) ) ? method : NULL;
}

/* set default methods for clustering. */
zVecClusterMethod *zVecClusterMethodSetDefault(zVecClusterMethod *method, int size)
{
  zVecClusterMethodInit( method );
  return zVecClusterMethodSetCoreFunc( method, size, NULL, NULL ) &&
         zVecClusterMethodSetErrorFunc( method, size, NULL, NULL ) &&
         zVecClusterMethodSetDistFunc( method, NULL, NULL ) &&
         zVecClusterMethodSetLoadedMeanFunc( method, NULL, NULL ) ? method : NULL;
}

/* set methods of linear-sum for clustering. */
zVecClusterMethod *zVecClusterMethodSetLS(zVecClusterMethod *method, void *util)
{
  return zVecClusterMethodSetCoreFunc( method, method->core_size, zVecClusterCoreLS, util ) &&
         zVecClusterMethodSetErrorFunc( method, 1, zVecClusterErrorLS, util ) &&
         zVecClusterMethodSetLoadedMeanFunc( method, zVecClusterLoadedMeanLS, util ) ? method : NULL;
}

/* copy methods for clustering. */
zVecClusterMethod *zVecClusterMethodCopy(const zVecClusterMethod *src, zVecClusterMethod *dest)
{
  return zVecClusterMethodSetCoreFunc( dest, src->core_size, src->core_fp, src->core_util ) &&
         zVecClusterMethodSetErrorFunc( dest, zVecSize(src->error), src->error_fp, src->error_util ) &&
         zVecClusterMethodSetDistFunc( dest, src->dist_fp, src->dist_util ) &&
         zVecClusterMethodSetLoadedMeanFunc( dest, src->lm_fp, src->lm_util ) ? dest : NULL;
}

/* destroy methods for clustering. */
void zVecClusterMethodDestroy(zVecClusterMethod *method)
{
  zVecFree( method->error );
  zVecFree( method->_lerr );
  zVecClusterMethodInit( method );
}

/* ********************************************************** */
/* multiple vecter clusters class.
 * ********************************************************** */

/* initialize multiple vector clusters. */
zVecMCluster *zVecMClusterInit(zVecMCluster *mc, int size)
{
  zListInit( zVecMClusterClusterList(mc) );
  return zVecClusterMethodSetDefault( &mc->method, size ) ? mc : NULL;
}

/* allocate multiple vector clusters. */
zVecMCluster *zVecMClusterAlloc(zVecMCluster *mc, int n)
{
  zVecClusterListCell *cc;
  int i;

  zListInit( zVecMClusterClusterList(mc) );
  if( !mc->method.core_fp ){
    ZRUNERROR( ZM_ERR_MVA_NOCOREFUNC );
    return NULL;
  }
  for( i=0; i<n; i++ ){
    if( !( cc = zAlloc( zVecClusterListCell, 1 ) ) ){
      ZALLOCERROR();
      goto ERR;
    }
    if( !zVecClusterCreate( &cc->data, mc->method.core_size ) ){
      ZALLOCERROR();
      free( cc );
      goto ERR;
    }
    zListInsertHead( zVecMClusterClusterList(mc), cc );
  }
  return mc;

 ERR:
  zVecMClusterDestroy( mc );
  return NULL;
}

/* move a multiple cluster to another. */
bool zVecMClusterMove(zVecMCluster *src, zVecMCluster *dest)
{
  if( !zVecClusterMethodCopy( &src->method, &dest->method ) ) return false;
  zListMove( zVecMClusterClusterList(src), zVecMClusterClusterList(dest) );
  zVecClusterMethodDestroy( &src->method );
  return true;
}

/* destroy multiple vector clusters. */
void zVecMClusterDestroy(zVecMCluster *mc)
{
  zVecClusterListCell *cc;

  while( !zListIsEmpty(zVecMClusterClusterList(mc)) ){
    zListDeleteHead( zVecMClusterClusterList(mc), &cc );
    zVecClusterDestroy( &cc->data );
    free( cc );
  }
  zVecClusterMethodDestroy( &mc->method );
}

/* evenness of clusters. */
double zVecMClusterEvenness(const zVecMCluster *mc)
{
  zVecClusterListCell *cp;
  int size_min = INT_MAX, size_max = 0;

  zListForEach( zVecMClusterClusterList(mc), cp ){
    if( zListSize(zVecClusterSampleList(&cp->data)) < size_min )
      size_min = zListSize(zVecClusterSampleList(&cp->data));
    if( zListSize(zVecClusterSampleList(&cp->data)) > size_max )
      size_max = zListSize(zVecClusterSampleList(&cp->data));
  }
  return (double)size_max / size_min;
}

/* print multiple vector clusters to a file. */
void zVecMClusterFPrint(FILE *fp, const zVecMCluster *mc)
{
  int i = 0;
  zVecClusterListCell *cc;

  zListForEach( zVecMClusterClusterList(mc), cc ){
    fprintf( fp, "#cluster[%d] : ", i++ );
    zVecClusterFPrint( fp, &cc->data );
  }
}

/* print data of multiple vector clusters to files. */
void zVecMClusterValueFPrint(FILE *fp[], const zVecMCluster *mc)
{
  zVecClusterListCell *cc;
  int i = 0;

  zListForEach( zVecMClusterClusterList(mc), cc )
    zVecClusterValueFPrint( fp[i++], &cc->data );
}

/* print cores of each cluster of multiple vector clusters to files. */
void zVecMClusterCoreFPrint(FILE *fp[], const zVecMCluster *mc)
{
  zVecClusterListCell *cc;
  int i = 0;

  zListForEach( zVecMClusterClusterList(mc), cc )
    zVecValueFPrint( fp[i++], cc->data.core );
}

/* print vectors in a set of clusters to files with a common basename. */
bool zVecMClusterValuePrintFile(const zVecMCluster *mc, const char *basename)
{
  zVecClusterListCell *cp;
  char filename[BUFSIZ];
  FILE *fp;
  int i = 0;

  zListForEach( zVecMClusterClusterList(mc), cp ){
    sprintf( filename, "%s%d", basename, i++ );
    if( !( fp = fopen( filename, "w" ) ) ){
      ZOPENERROR( filename );
      return false;
    }
    zVecClusterValueFPrint( fp, &cp->data );
    fclose( fp );
  }
  return true;
}

/* ********************************************************** */
/* clustering based on K-means family
 * ********************************************************** */

#if DEBUG
/* for debug */
static void _zVecMClusterKMeansValuePrintFile(const zVecMCluster *mc, int step)
{
  FILE *fp;
  char filename[BUFSIZ];
  zVecClusterListCell *cc;
  int i;

  i = 0;
  zListForEach( zVecMClusterClusterList(mc), cc ){
    sprintf( filename, "%d_%d", step, i );
    fp = fopen( filename, "w" );
    zVecClusterValueFPrint( fp, &cc->data );
    fclose( fp );
    sprintf( filename, "%d_%dm", step, i );
    fp = fopen( filename, "w" );
    zVecValueFPrint( fp, cc->data.core );
    fclose( fp );
  }
}
#endif /* DEBUG */

/* compute variance of a vector cluster. */
static double _zVecMClusterVar(const zVecMCluster *mc, zVecCluster *c)
{
  zVecListCell *pc;

  c->var = 0;
  zListForEach( zVecClusterSampleList(c), pc ){
    zVecMClusterErrorF( mc, pc->data, c->core );
    c->var += zVecSqrNorm( mc->method.error );
  }
  return ( c->var /= zListSize(zVecClusterSampleList(c)) );
}

/* assign all points to core points of clusters. */
static bool _zVecMClusterKMeansInitCluster(zVecMCluster *mc, const zVecAddrList *points)
{
  zVecListCell *pp;
  double d, dmin;
  zVecClusterListCell *cp, *ccp;

  zListForEach( points, pp ){
    dmin = HUGE_VAL;
    ccp = NULL;
    zListForEach( zVecMClusterClusterList(mc), cp ){
      if( pp->data == zListTail(zVecClusterSampleList(&cp->data))->data ) goto CONTINUE;
      if( ( d = zVecMClusterDistF( mc, pp->data, zListTail(zVecClusterSampleList(&cp->data))->data ) ) < dmin ){
        dmin = d;
        ccp = cp;
      }
    }
    if( !zVecAddrListInsertHead( zVecClusterSampleList(&ccp->data), pp->data ) ) return false;
    CONTINUE: ;
  }
  return true;
}

/* arrange initial candidates of clusters for K-means based on KKZ method. */
static bool _zVecMClusterKMeansInitKKZ(zVecMCluster *mc, const zVecAddrList *points)
{
  zVecListCell *pp, *pp_max = NULL;
  zVecClusterListCell *cp, *ccp;
  zVec mean;
  double d, dmax, dmin;

  if( zListIsEmpty( points ) ){
    ZRUNERROR( ZM_ERR_MVA_EMPTYSET );
    return false;
  }
  mean = zVecAlloc( zVecSizeNC( zListHead(points)->data ) );
  zVecListMean( points, mean );
  /* first core point */
  dmax = 0;
  zListForEach( points, pp ){
    if( ( d = zVecMClusterDistF( mc, pp->data, mean ) ) > dmax ){
      dmax = d;
      pp_max = pp;
    }
  }
  cp = zListTail( zVecMClusterClusterList(mc) );
  if( !zVecAddrListInsertTail( zVecClusterSampleList(&cp->data), pp_max->data ) ) return false;
  zVecFree( mean );
  /* second - k-th mean point */
  for( cp=zListCellNext(cp); cp!=zListRoot(zVecMClusterClusterList(mc)); cp=zListCellNext(cp) ){
    dmax = 0;
    pp_max = NULL;
    zListForEach( points, pp ){
      dmin = HUGE_VAL;
      for( ccp=zListTail(zVecMClusterClusterList(mc)); ccp!=cp; ccp=zListCellNext(ccp) )
        if( ( d = zVecMClusterDistF( mc, pp->data, zListTail(zVecClusterSampleList(&ccp->data))->data ) ) < dmin )
          dmin = d;
      if( dmin > dmax ){
        dmax = dmin;
        pp_max = pp;
      }
    }
    if( !zVecAddrListInsertTail( zVecClusterSampleList(&cp->data), pp_max->data ) ) return false;
  }
  return _zVecMClusterKMeansInitCluster( mc, points );
}

/* arrange initial candidates of clusters for K-means based on K-means++ method. */
static bool _zVecMClusterKMeansInitPP(zVecMCluster *mc, const zVecAddrList *points)
{
  zVecListCell *pp;
  zVecClusterListCell *cp, *ccp;
  zVec score;
  double d, dmin, p, rate;
  int i;
  bool ret = false;

  if( zListIsEmpty( points ) ){
    ZRUNERROR( ZM_ERR_MVA_EMPTYSET );
    return false;
  }
  if( !( score = zVecAlloc( zListSize(points)) ) ) return false;
  /* first core point */
  i = zRandI( 0, zListSize(points)-1 );
  zListItem( points, i, &pp );
  cp = zListTail( zVecMClusterClusterList(mc) );
  if( !zVecAddrListInsertTail( zVecClusterSampleList(&cp->data), pp->data ) ) goto TERMINATE;
  /* second - k-th core point */
  for( cp=zListCellNext(cp); cp!=zListRoot(zVecMClusterClusterList(mc)); cp=zListCellNext(cp) ){
    i = 0;
    zListForEach( points, pp ){
      for( dmin=HUGE_VAL, ccp=zListTail(zVecMClusterClusterList(mc)); ccp!=cp; ccp=zListCellNext(ccp) ){
        if( ( d = zVecMClusterDistF( mc, pp->data, zListTail(zVecClusterSampleList(&ccp->data))->data ) ) < dmin )
          dmin = d;
      }
      zVecSetElemNC( score, i++, zSqr(dmin) );
    }
    zVecDivDRC( score, zVecElemSum( score ));
    /* roulette */
    p = zRandF( 0, 1 );
    rate = 0;
    i = 0;
    zListForEach( points, pp ){
      if( ( rate += zVecElemNC(score,i) ) >= p ) break;
      i++;
    }
    if( !zVecAddrListInsertTail( zVecClusterSampleList(&cp->data), pp->data ) ) goto TERMINATE;
  }
  ret = _zVecMClusterKMeansInitCluster( mc, points );
 TERMINATE:
  zVecFree( score );
  return ret;
}

/* recluster tentative clusters for K-means. */
static int _zVecMClusterKMeansRecluster(zVecMCluster *mc)
{
  zVecClusterListCell *cc1, *cc2, *cc;
  zVecListCell *pc, *pc_prev;
  bool ismoved;
  double d, dmin;
  int i, iter = 0;

  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    /* compute core of clusters */
    zListForEach( zVecMClusterClusterList(mc), cc1 ){
      zVecMClusterCoreF( mc, zVecClusterSampleList(&cc1->data), cc1->data.core );
      _zVecMClusterVar( mc, &cc1->data );
    }
    ismoved = false;
    zListForEach( zVecMClusterClusterList(mc), cc1 ){
      /* for each tentative cluster */
      zListForEach( zVecClusterSampleList(&cc1->data), pc ){ /* for each point */
        dmin = HUGE_VAL;
        cc = cc1;
        zListForEach( zVecMClusterClusterList(mc), cc2 ){ /* find the nearest core point */
          zVecMClusterErrorF( mc, pc->data, cc2->data.core );
          if( ( d = zVecSqrNorm( mc->method.error ) ) < dmin ){
            dmin = d;
            cc = cc2;
          }
        }
        if( cc != cc1 ){ /* move the point into nearer cluster */
          ismoved = true;
          pc_prev = zListCellPrev(pc);
          zListPurge( zVecClusterSampleList(&cc1->data), pc );
          zListInsertHead( zVecClusterSampleList(&cc->data), pc );
          pc = pc_prev;
        }
      }
    }
    if( ismoved == false ) return i; /* clusters settled */
  }
  ZITERWARN( iter );
  return iter;
}

/* clustering of vectors by K-means with KKZ initialization. */
int zVecMClusterKMeansKKZ(zVecMCluster *mc, const zVecAddrList *points, int k)
{
  if( !zVecMClusterAlloc( mc, k ) ) return -1;
  if( !_zVecMClusterKMeansInitKKZ( mc, points ) ) return -1;
  return _zVecMClusterKMeansRecluster( mc );
}

/* clustering of vectors by K-means++. */
int zVecMClusterKMeans(zVecMCluster *mc, const zVecAddrList *points, int k)
{
  if( !zVecMClusterAlloc( mc, k ) ) return -1;
  if( !_zVecMClusterKMeansInitPP( mc, points ) ) return -1;
  return _zVecMClusterKMeansRecluster( mc );
}

/* clustering of vectors by K-medoids. */
int zVecMClusterKMedoids(zVecMCluster *mc, const zVecAddrList *points, int k)
{
  zVecMClusterSetCoreFunc( mc, zVecSizeNC(zListTail(points)->data), _zVecClusterCoreMedoid, NULL );
  return zVecMClusterKMeans( mc, points, k );
}

/* ********************************************************** */
/* silhouette analysis
 * ********************************************************** */

/* compute the mean silhouette of a set of clusters. */
double zVecMClusterMeanSilhouette(zVecMCluster *mc)
{
  zVecClusterListCell *cp, *ccp;
  zVecAddrListCell *pp;
  double a, b, b_tmp, score = 0;
  int i, n = 0;

  zListForEach( zVecMClusterClusterList(mc), cp ){
    zFree( cp->data._sil );
    if( !( cp->data._sil = zAlloc( double, zListSize(zVecClusterSampleList(&cp->data)) ) ) ){
      ZALLOCERROR();
      return NAN;
    }
    i = 0;
    zListForEach( zVecClusterSampleList(&cp->data), pp ){
      a = _zVecClusterDistAve( &mc->method, zVecClusterSampleList(&cp->data), pp->data ); /* intra-cluster distance */
      b = HUGE_VAL;
      zListForEach( zVecMClusterClusterList(mc), ccp ){
        if( ccp == cp ) continue;
        if( ( b_tmp = _zVecClusterDistAve( &mc->method, zVecClusterSampleList(&ccp->data), pp->data ) ) < b )
          b = b_tmp; /* inter-cluster distance */
      }
      score += cp->data._sil[i++] = ( b - a ) / zMax( a, b ); /* silhouette of the sample */
    }
    zDataSort( cp->data._sil, zListSize(zVecClusterSampleList(&cp->data)) );
    n += zListSize(zVecClusterSampleList(&cp->data));
  }
  return score / n;
}

/* print silhouettes of a vector cluster to a file. */
static bool _zVecClusterSilhouetteFPrint(FILE *fp, const zVecCluster *c, int offset)
{
  int i;

  if( !c->_sil ){
    ZRUNWARN( ZM_WARN_MVA_NOSILHOUETTE );
    return false;
  }
  for( i=0; i<zListSize(zVecClusterSampleList(c)); i++ ){
    fprintf( fp, "0, %d\n", offset + i );
    fprintf( fp, "%f, %d\n\n", c->_sil[i], offset + i );
  }
  return true;
}

/* print silhouettes of a set of vector clusters to files. */
bool zVecMClusterSilhouettePrintFile(const zVecMCluster *mc, const char *basename)
{
  zVecClusterListCell *cp;
  int i = 0, offset = 0;
  char filename[BUFSIZ];
  FILE *fp;

  zListForEach( zVecMClusterClusterList(mc), cp ){
    sprintf( filename, "%s%d", basename, i++ );
    if( !( fp = fopen( filename, "w" ) ) ){
      ZOPENERROR( filename );
      return false;
    }
    if( !_zVecClusterSilhouetteFPrint( fp, &cp->data, offset ) ) return false;
    offset += zListSize(zVecClusterSampleList(&cp->data));
  }
  return true;
}

/* ********************************************************** */
/*! \brief clustering based on X-means
 *//* ******************************************************* */

#define Z_XMEANS_MINSIZE 10

/* merge subclusters with original. */
static void _zVecMClusterMerge(zVecMCluster *mc, zVecClusterListCell *vc, zVecMCluster *submc)
{
  zVecClusterListCell *prev, *next;

  prev = zListCellPrev(vc);
  next = zListCellNext(vc);
  zVecClusterDestroy( &vc->data );
  free( vc );
  zListCellBind( prev, zListTail(zVecMClusterClusterList(submc)) );
  zListCellBind( zListHead(zVecMClusterClusterList(submc)), next );
  zListSize(zVecMClusterClusterList(mc)) += zListSize(zVecMClusterClusterList(submc)) - 1;
  zListInit( zVecMClusterClusterList(submc) );
  zVecMClusterDestroy( submc );
}

/* initialize clusters for X-means. */
static bool _zVecMClusterXMeansInit(zVecMCluster *mc, const zVecAddrList *points)
{
  zVecClusterListCell *vc;
  zVecListCell *pc;

  vc = zListHead(zVecMClusterClusterList(mc));
  zListForEach( points, pc )
    if( !zVecAddrListInsertHead( zVecClusterSampleList(&vc->data), pc->data ) ) return false;
  zVecMClusterCoreF( mc, zVecClusterSampleList(&vc->data), vc->data.core );
  _zVecMClusterVar( mc, &vc->data );
  return true;
}

/* find maximum error in a cluster from the mean value. */
static double _zVecClusterMaxError(const zVecCluster *c, const zVecClusterMethod *method)
{
  zVecListCell *vc;
  double dmax = 0, d;

  zListForEach( zVecClusterSampleList(c), vc ){
    zVecClusterMethodErrorF( method, vc->data, c->core );
    if( ( d = zVecNorm( method->error ) ) > dmax ) dmax = d;
  }
  return dmax;
}

/* check if the subclusters better classify samples than the original based on hyperdensity of clusters. */
static bool _zVecMClusterXMeansCheckDensity(const zVecClusterMethod *method, const zVecCluster *c, const zVecCluster *c1, const zVecCluster *c2)
{
  double rm, rm1, rm2;

  /* original cluster */
  rm = pow( _zVecClusterMaxError( c, method ), zVecSizeNC(method->error) );
  /* bidivided clusters */
  rm1 = pow( _zVecClusterMaxError( c1, method ), zVecSizeNC(method->error) );
  rm2 = pow( _zVecClusterMaxError( c2, method ), zVecSizeNC(method->error) );
  return rm > rm1 + rm2 ? true : false;
}

/* probability density function based on Gaussian distribution. */
static double _zVecClusterXMeansPDF(const zVecClusterMethod *method, const zVecCluster *c, const zVec x)
{
  zVecClusterMethodErrorF( method, x, c->core );
  return exp( -0.5*zVecSqrNorm(method->error)/c->var ) / sqrt(c->var);
}

/* check if the subclusters better classify samples than the original based on BIC. */
static bool _zVecMClusterXMeansCheckBIC(const zVecClusterMethod *method, const zVecCluster *c, const zVecCluster *c1, const zVecCluster *c2)
{
  zVecListCell *vc;
  int n;
  double ls = 0, bic1, bic2;

  n = zListSize(zVecClusterSampleList(c));
  zListForEach( zVecClusterSampleList(c1), vc ){
    ls += log( _zVecClusterXMeansPDF( method, c1, vc->data )
             + _zVecClusterXMeansPDF( method, c2, vc->data ) );
  }
  zListForEach( zVecClusterSampleList(c2), vc ){
    ls += log( _zVecClusterXMeansPDF( method, c1, vc->data )
             + _zVecClusterXMeansPDF( method, c2, vc->data ) );
  }
  bic1 = n * ( log(c->var) + 1 );
  bic2 = 2*( n*log(2) - ls ) + method->core_size * log(n);
  return bic1 < bic2 ? true : false;
}

/* test if the subclusters better classify samples than the original. */
static bool _zVecMClusterXMeansTest(const zVecMCluster *mc, const zVecCluster *c, const zVecMCluster *submc, bool (* testfunc)(const zVecClusterMethod*,const zVecCluster*,const zVecCluster*,const zVecCluster*))
{
  zVecCluster *c1, *c2;

  c1 = &zListHead(zVecMClusterClusterList(submc))->data;
  c2 = &zListTail(zVecMClusterClusterList(submc))->data;
  if( 10 * zListSize(zVecClusterSampleList(c1)) < zListSize(zVecClusterSampleList(c2)) ||
      10 * zListSize(zVecClusterSampleList(c2)) < zListSize(zVecClusterSampleList(c1)) )
    return false; /* avoid too crisp cluster */
  return testfunc( &mc->method, c, c1, c2 );
}

/* internal function for recursive all of X-means method. */
static int _zVecMClusterXMeans(zVecMCluster *mc, zVecClusterListCell *cc, bool (* testfunc)(const zVecClusterMethod*,const zVecCluster*,const zVecCluster*,const zVecCluster*))
{
  zVecMCluster submc;
  int iter, iter1, iter2;

  /* ignore a too small cluster */
  if( zListSize(zVecClusterSampleList(&cc->data)) < Z_XMEANS_MINSIZE ) return 0;
  if( !zVecMClusterInit( &submc, mc->method.core_size ) ||
      !zVecMClusterMethodCopy( mc, &submc ) ) return 0;

  iter = zVecMClusterKMeans( &submc, zVecClusterSampleList(&cc->data), 2 );
  if( _zVecMClusterXMeansTest( mc, &cc->data, &submc, testfunc ) ){
    zVecMClusterDestroy( &submc );
    return 0;
  }
  iter1 = _zVecMClusterXMeans( &submc, zListHead(zVecMClusterClusterList(&submc)), testfunc );
  iter2 = _zVecMClusterXMeans( &submc, zListTail(zVecMClusterClusterList(&submc)), testfunc );
  _zVecMClusterMerge( mc, cc, &submc );
  return iter + iter1 + iter2;
}

/* cluster vectors by X-means method. */
int zVecMClusterXMeans(zVecMCluster *mc, const zVecAddrList *points, bool (* testfunc)(const zVecClusterMethod*,const zVecCluster*,const zVecCluster*,const zVecCluster*))
{
  if( !zVecMClusterAlloc( mc, 1 ) ) return -1;
  if( !_zVecMClusterXMeansInit( mc, points ) ) return -1;
  return _zVecMClusterXMeans( mc, zListHead(zVecMClusterClusterList(mc)), testfunc );
}

/* cluster vectors by X-means method based on hyperdensity. */
int zVecMClusterXMeansDensity(zVecMCluster *mc, const zVecAddrList *points)
{
  return zVecMClusterXMeans( mc, points, _zVecMClusterXMeansCheckDensity );
}

/* cluster vectors by X-means method based on BIC. */
int zVecMClusterXMeansBIC(zVecMCluster *mc, const zVecAddrList *points)
{
  return zVecMClusterXMeans( mc, points, _zVecMClusterXMeansCheckBIC );
}
