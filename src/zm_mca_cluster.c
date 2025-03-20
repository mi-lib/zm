/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_mca_cluster - multiple classification analysis: clustering.
 */

#include <zm/zm_mca.h>

/* ********************************************************** */
/* vector cluster class.
 * ********************************************************** */

/* create and initialize a vector cluster. */
zCluster *zClusterCreate(zCluster *c, int coresize)
{
  zListInit( zClusterSampleList(c) );
  if( !( c->core = zVecAlloc( coresize ) ) ){
    zClusterDestroy( c );
    return NULL;
  }
  c->var = HUGE_VAL;
  c->_sil = NULL;
  return c;
}

/* destroy a vector cluster. */
void zClusterDestroy(zCluster *c)
{
  zVecAddrListDestroy( zClusterSampleList(c) );
  zVecFree( c->core );
  zFree( c->_sil );
}

/* the maximum silhouette in a cluster. */
double zClusterMaxSilhouette(zCluster *c)
{
  return c && c->_sil ? c->_sil[zListSize(zClusterSampleList(c))-1] : -HUGE_VAL;
}

/* print a vector cluster to a file. */
void zClusterFPrint(FILE *fp, zCluster *c)
{
  int i = 0;
  zVecListCell *vc;

  fprintf( fp, "%d members\n", zListSize(zClusterSampleList(c)) );
  zListForEach( zClusterSampleList(c), vc ){
    fprintf( fp, " %d: ", i++ );
    zVecFPrint( fp, vc->data );
  }
  fprintf( fp, " core: " );
  zVecFPrint( fp, c->core );
}

/* print data of vectors of a vector cluster to a file. */
void zClusterValueFPrint(FILE *fp, zCluster *c)
{
  zVecListFPrint( fp, zClusterSampleList(c) );
}

/* ********************************************************** */
/* methods for core and error computation
 * ********************************************************** */

/* a default function to compute error of a sample for clustering. */
static zVec _zClusterErrorDefault(zClusterMethod *method, zVec p, zVec core, void *dummy, zVec err)
{
  return zVecSub( p, core, err );
}

/* a default function to compute the distance between two samples for clustering. */
static double _zClusterDistDefault(zClusterMethod *method, zVec p1, zVec p2, void *dummy)
{
  return zVecDist( p1, p2 );
}

/* a default function to find core (typically, mean) for clustering. */
static zVec _zClusterCoreDefault(zClusterMethod *method, zVecAddrList *pl, void *dummy, zVec core)
{
  return zVecListMean( pl, core );
}

/* a function to find medoid as the core for clustering. */
static double _zClusterDistSum(zClusterMethod *cm, zVecAddrList *pl, zVec v)
{
  zVecAddrListCell *pp;
  double sum = 0;

  zListForEach( pl, pp )
    sum += zClusterMethodDistF( cm, pp->data, v );
  return sum;
}
static double _zClusterDistAve(zClusterMethod *cm, zVecAddrList *pl, zVec v)
{
  return _zClusterDistSum( cm, pl, v ) / zListSize(pl);
}
static zVec _zClusterCoreMedoid(zClusterMethod *cm, zVecAddrList *pl, void *dummy, zVec core)
{
  zVecAddrListCell *cp, *mp;
  double d, dmin = HUGE_VAL;

  mp = zListTail( pl );
  zListForEach( pl, cp ){
    if( ( d = _zClusterDistSum( cm, pl, cp->data ) ) < dmin ){
      dmin = d;
      mp = cp;
    }
  }
  return zVecCopy( mp->data, core );
}

/* a default loaded mean computation function for clustering. */
static zVec _zClusterMeanLoadedDefault(zClusterMethod *method, zVecAddrList *pl, double load[], double n, void *dummy, zVec mean)
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

/* initialize methods for clustering. */
void zClusterMethodInit(zClusterMethod *method)
{
  /* core of a cluster */
  method->core_size = 0;
  method->core_fp = NULL;
  method->core_util = NULL;
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
zClusterMethod *zClusterMethodSetCoreFunc(zClusterMethod *method, int size, zVec (* fp)(zClusterMethod*,zVecAddrList*,void*,zVec), void *util)
{
  if( size <= 0 ){
    ZRUNERROR( ZM_ERR_MCA_INVALID_SIZE, size );
    return NULL;
  }
  method->core_size = size;
  method->core_fp = fp ? fp : _zClusterCoreDefault;
  method->core_util = util;
  return method;
}

/* set a function to compute error of a sample in a cluster. */
zClusterMethod *zClusterMethodSetErrorFunc(zClusterMethod *method, int size, zVec (* fp)(zClusterMethod*,zVec,zVec,void*,zVec), void *util)
{
  if( method->error ) zVecFree( method->error );
  if( size <= 0 ){
    ZRUNERROR( ZM_ERR_MCA_INVALID_SIZE, size );
    return NULL;
  }
  method->error_fp = fp ? fp : _zClusterErrorDefault;
  method->error_util = util;
  return ( method->error = zVecAlloc( size ) ) ? method : NULL;
}

/* set a function to compute the distance between two samples in a cluster. */
zClusterMethod *zClusterMethodSetDistFunc(zClusterMethod *method, double (* fp)(zClusterMethod*,zVec,zVec,void*), void *util)
{
  method->dist_fp = fp ? fp : _zClusterDistDefault;
  method->dist_util = util;
  return method;
}

/* set a function to compute the loaded mean of samples (for GMM). */
zClusterMethod *zClusterMethodSetLoadedMeanFunc(zClusterMethod *method, zVec (* fp)(zClusterMethod*,zVecAddrList*,double[],double,void*,zVec), void *util)
{
  if( !method->core_fp ){
    ZRUNERROR( ZM_ERR_MCA_NOCOREFUNC );
    return NULL;
  }
  if( !method->error_fp ){
    ZRUNERROR( ZM_ERR_MCA_NOERRORFUNC );
    return NULL;
  }
  if( !method->error_fp ){
    ZRUNERROR( ZM_ERR_MCA_NODISTFUNC );
    return NULL;
  }
  if( method->_lerr ) zVecFree( method->_lerr );
  method->lm_fp = fp ? fp : _zClusterMeanLoadedDefault;
  method->lm_util = util;
  return ( method->_lerr = zVecAlloc( zVecSize(method->error) ) ) ? method : NULL;
}

/* create methods for clustering. */
zClusterMethod *zClusterMethodCreate(zClusterMethod *method, int size)
{
  zClusterMethodInit( method );
  return zClusterMethodSetCoreFunc( method, size, NULL, NULL ) &&
         zClusterMethodSetErrorFunc( method, size, NULL, NULL ) &&
         zClusterMethodSetDistFunc( method, NULL, NULL ) &&
         zClusterMethodSetLoadedMeanFunc( method, NULL, NULL ) ? method : NULL;
}

/* copy methods for clustering. */
zClusterMethod *zClusterMethodCopy(zClusterMethod *src, zClusterMethod *dest)
{
  return zClusterMethodSetCoreFunc( dest, src->core_size, src->core_fp, src->core_util ) &&
         zClusterMethodSetErrorFunc( dest, zVecSize(src->error), src->error_fp, src->error_util ) &&
         zClusterMethodSetDistFunc( dest, src->dist_fp, src->dist_util ) &&
         zClusterMethodSetLoadedMeanFunc( dest, src->lm_fp, src->lm_util ) ? dest : NULL;
}

/* destroy methods for clustering. */
void zClusterMethodDestroy(zClusterMethod *method)
{
  zVecFree( method->error );
  zVecFree( method->_lerr );
  zClusterMethodInit( method );
}

/* ********************************************************** */
/* multiple vecter clusters class.
 * ********************************************************** */

/* initialize multiple vector clusters. */
zMCluster *zMClusterInit(zMCluster *mc, int size)
{
  zListInit( zMClusterClusterList(mc) );
  return zClusterMethodCreate( &mc->method, size ) ? mc : NULL;
}

/* allocate multiple vector clusters. */
zMCluster *zMClusterAlloc(zMCluster *mc, int n)
{
  zClusterListCell *cc;
  int i;

  zListInit( zMClusterClusterList(mc) );
  if( !mc->method.core_fp ){
    ZRUNERROR( ZM_ERR_MCA_NOCOREFUNC );
    return NULL;
  }
  for( i=0; i<n; i++ ){
    if( !( cc = zAlloc( zClusterListCell, 1 ) ) ){
      ZALLOCERROR();
      goto ERR;
    }
    if( !zClusterCreate( &cc->data, mc->method.core_size ) ){
      ZALLOCERROR();
      free( cc );
      goto ERR;
    }
    zListInsertHead( zMClusterClusterList(mc), cc );
  }
  return mc;

 ERR:
  zMClusterDestroy( mc );
  return NULL;
}

/* move a multiple cluster to another. */
bool zMClusterMove(zMCluster *src, zMCluster *dest)
{
  if( !zClusterMethodCopy( &src->method, &dest->method ) ) return false;
  zListMove( zMClusterClusterList(src), zMClusterClusterList(dest) );
  zClusterMethodDestroy( &src->method );
  return true;
}

/* destroy multiple vector clusters. */
void zMClusterDestroy(zMCluster *mc)
{
  zClusterListCell *cc;

  while( !zListIsEmpty(zMClusterClusterList(mc)) ){
    zListDeleteHead( zMClusterClusterList(mc), &cc );
    zClusterDestroy( &cc->data );
    free( cc );
  }
  zClusterMethodDestroy( &mc->method );
}

/* evenness of clusters. */
double zMClusterEvenness(zMCluster *mc)
{
  zClusterListCell *cp;
  int size_min = INT_MAX, size_max = 0;

  zListForEach( zMClusterClusterList(mc), cp ){
    if( zListSize(zClusterSampleList(&cp->data)) < size_min )
      size_min = zListSize(zClusterSampleList(&cp->data));
    if( zListSize(zClusterSampleList(&cp->data)) > size_max )
      size_max = zListSize(zClusterSampleList(&cp->data));
  }
  return (double)size_max / size_min;
}

/* print multiple vector clusters to a file. */
void zMClusterFPrint(FILE *fp, zMCluster *mc)
{
  int i = 0;
  zClusterListCell *cc;

  zListForEach( zMClusterClusterList(mc), cc ){
    fprintf( fp, "#cluster[%d] : ", i++ );
    zClusterFPrint( fp, &cc->data );
  }
}

/* print data of multiple vector clusters to files. */
void zMClusterValueFPrint(FILE *fp[], zMCluster *mc)
{
  zClusterListCell *cc;
  int i = 0;

  zListForEach( zMClusterClusterList(mc), cc )
    zClusterValueFPrint( fp[i++], &cc->data );
}

/* print cores of each cluster of multiple vector clusters to files. */
void zMClusterCoreFPrint(FILE *fp[], zMCluster *mc)
{
  zClusterListCell *cc;
  int i = 0;

  zListForEach( zMClusterClusterList(mc), cc )
    zVecValueFPrint( fp[i++], cc->data.core );
}

/* print vectors in a set of clusters to files with a common basename. */
bool zMClusterValuePrintFile(zMCluster *mc, const char *basename)
{
  zClusterListCell *cp;
  char filename[BUFSIZ];
  FILE *fp;
  int i = 0;

  zListForEach( zMClusterClusterList(mc), cp ){
    sprintf( filename, "%s%d", basename, i++ );
    if( !( fp = fopen( filename, "w" ) ) ){
      ZOPENERROR( filename );
      return false;
    }
    zClusterValueFPrint( fp, &cp->data );
    fclose( fp );
  }
  return true;
}

/* ********************************************************** */
/* clustering based on K-means family
 * ********************************************************** */

#if DEBUG
/* for debug */
static void _zMClusterKMeansValuePrintFile(zMCluster *mc, int step)
{
  FILE *fp;
  char filename[BUFSIZ];
  zClusterListCell *cc;
  int i;

  i = 0;
  zListForEach( zMClusterClusterList(mc), cc ){
    sprintf( filename, "%d_%d", step, i );
    fp = fopen( filename, "w" );
    zClusterValueFPrint( fp, &cc->data );
    fclose( fp );
    sprintf( filename, "%d_%dm", step, i );
    fp = fopen( filename, "w" );
    zVecValueFPrint( fp, cc->data.core );
    fclose( fp );
  }
}
#endif /* DEBUG */

/* compute variance of a vector cluster. */
static double _zMClusterVar(zMCluster *mc, zCluster *c)
{
  zVecListCell *pc;

  c->var = 0;
  zListForEach( zClusterSampleList(c), pc ){
    zMClusterErrorF( mc, pc->data, c->core );
    c->var += zVecSqrNorm( mc->method.error );
  }
  return ( c->var /= zListSize(zClusterSampleList(c)) );
}

/* assign all points to core points of clusters. */
static bool _zMClusterKMeansInitCluster(zMCluster *mc, zVecAddrList *points)
{
  zVecListCell *pp;
  double d, dmin;
  zClusterListCell *cp, *ccp;

  zListForEach( points, pp ){
    dmin = HUGE_VAL;
    ccp = NULL;
    zListForEach( zMClusterClusterList(mc), cp ){
      if( pp->data == zListTail(zClusterSampleList(&cp->data))->data ) goto CONTINUE;
      if( ( d = zMClusterDistF( mc, pp->data, zListTail(zClusterSampleList(&cp->data))->data ) ) < dmin ){
        dmin = d;
        ccp = cp;
      }
    }
    if( !zVecAddrListInsertHead( zClusterSampleList(&ccp->data), pp->data ) ) return false;
    CONTINUE: ;
  }
  return true;
}

/* arrange initial candidates of clusters for K-means based on KKZ method. */
static bool _zMClusterKMeansInitKKZ(zMCluster *mc, zVecAddrList *points)
{
  zVecListCell *pp, *pp_max = NULL;
  zClusterListCell *cp, *ccp;
  zVec mean;
  double d, dmax, dmin;

  if( zListIsEmpty( points ) ){
    ZRUNERROR( ZM_ERR_MCA_EMPTYSET );
    return false;
  }
  mean = zVecAlloc( zVecSizeNC( zListHead(points)->data ) );
  zVecListMean( points, mean );
  /* first core point */
  dmax = 0;
  zListForEach( points, pp ){
    if( ( d = zMClusterDistF( mc, pp->data, mean ) ) > dmax ){
      dmax = d;
      pp_max = pp;
    }
  }
  cp = zListTail( zMClusterClusterList(mc) );
  if( !zVecAddrListInsertTail( zClusterSampleList(&cp->data), pp_max->data ) ) return false;
  zVecFree( mean );
  /* second - k-th mean point */
  for( cp=zListCellNext(cp); cp!=zListRoot(zMClusterClusterList(mc)); cp=zListCellNext(cp) ){
    dmax = 0;
    pp_max = NULL;
    zListForEach( points, pp ){
      dmin = HUGE_VAL;
      for( ccp=zListTail(zMClusterClusterList(mc)); ccp!=cp; ccp=zListCellNext(ccp) )
        if( ( d = zMClusterDistF( mc, pp->data, zListTail(zClusterSampleList(&ccp->data))->data ) ) < dmin )
          dmin = d;
      if( dmin > dmax ){
        dmax = dmin;
        pp_max = pp;
      }
    }
    if( !zVecAddrListInsertTail( zClusterSampleList(&cp->data), pp_max->data ) ) return false;
  }
  return _zMClusterKMeansInitCluster( mc, points );
}

/* arrange initial candidates of clusters for K-means based on K-means++ method. */
static bool _zMClusterKMeansInitPP(zMCluster *mc, zVecAddrList *points)
{
  zVecListCell *pp;
  zClusterListCell *cp, *ccp;
  zVec score;
  double d, dmin, p, rate;
  int i;
  bool ret = false;

  if( zListIsEmpty( points ) ){
    ZRUNERROR( ZM_ERR_MCA_EMPTYSET );
    return false;
  }
  if( !( score = zVecAlloc( zListSize(points)) ) ) return false;
  /* first core point */
  i = zRandI( 0, zListSize(points)-1 );
  zListItem( points, i, &pp );
  cp = zListTail( zMClusterClusterList(mc) );
  if( !zVecAddrListInsertTail( zClusterSampleList(&cp->data), pp->data ) ) goto TERMINATE;
  /* second - k-th core point */
  for( cp=zListCellNext(cp); cp!=zListRoot(zMClusterClusterList(mc)); cp=zListCellNext(cp) ){
    i = 0;
    zListForEach( points, pp ){
      for( dmin=HUGE_VAL, ccp=zListTail(zMClusterClusterList(mc)); ccp!=cp; ccp=zListCellNext(ccp) ){
        if( ( d = zMClusterDistF( mc, pp->data, zListTail(zClusterSampleList(&ccp->data))->data ) ) < dmin )
          dmin = d;
      }
      zVecSetElemNC( score, i++, zSqr(dmin) );
    }
    zVecDivDRC( score, zVecSumElem( score ));
    /* roulette */
    p = zRandF( 0, 1 );
    rate = 0;
    i = 0;
    zListForEach( points, pp ){
      if( ( rate += zVecElemNC(score,i) ) >= p ) break;
      i++;
    }
    if( !zVecAddrListInsertTail( zClusterSampleList(&cp->data), pp->data ) ) goto TERMINATE;
  }
  ret = _zMClusterKMeansInitCluster( mc, points );
 TERMINATE:
  zVecFree( score );
  return ret;
}

/* recluster tentative clusters for K-means. */
static int _zMClusterKMeansRecluster(zMCluster *mc)
{
  zClusterListCell *cc1, *cc2, *cc;
  zVecListCell *pc, *pc_prev;
  bool ismoved;
  double d, dmin;
  int i, iter = 0;

  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    /* compute core of clusters */
    zListForEach( zMClusterClusterList(mc), cc1 ){
      zMClusterCoreF( mc, zClusterSampleList(&cc1->data), cc1->data.core );
      _zMClusterVar( mc, &cc1->data );
    }
    ismoved = false;
    zListForEach( zMClusterClusterList(mc), cc1 ){
      /* for each tentative cluster */
      zListForEach( zClusterSampleList(&cc1->data), pc ){ /* for each point */
        dmin = HUGE_VAL;
        cc = cc1;
        zListForEach( zMClusterClusterList(mc), cc2 ){ /* find the nearest core point */
          zMClusterErrorF( mc, pc->data, cc2->data.core );
          if( ( d = zVecSqrNorm( mc->method.error ) ) < dmin ){
            dmin = d;
            cc = cc2;
          }
        }
        if( cc != cc1 ){ /* move the point into nearer cluster */
          ismoved = true;
          pc_prev = zListCellPrev(pc);
          zListPurge( zClusterSampleList(&cc1->data), pc );
          zListInsertHead( zClusterSampleList(&cc->data), pc );
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
int zMClusterKMeansKKZ(zMCluster *mc, zVecAddrList *points, int k)
{
  if( !zMClusterAlloc( mc, k ) ) return -1;
  if( !_zMClusterKMeansInitKKZ( mc, points ) ) return -1;
  return _zMClusterKMeansRecluster( mc );
}

/* clustering of vectors by K-means++. */
int zMClusterKMeans(zMCluster *mc, zVecAddrList *points, int k)
{
  if( !zMClusterAlloc( mc, k ) ) return -1;
  if( !_zMClusterKMeansInitPP( mc, points ) ) return -1;
  return _zMClusterKMeansRecluster( mc );
}

/* clustering of vectors by K-medoids. */
int zMClusterKMedoids(zMCluster *mc, zVecAddrList *points, int k)
{
  zMClusterSetCoreFunc( mc, zVecSizeNC(zListTail(points)->data), _zClusterCoreMedoid, NULL );
  return zMClusterKMeans( mc, points, k );
}

/* ********************************************************** */
/* silhouette analysis
 * ********************************************************** */

/* compute the mean silhouette of a set of clusters. */
double zMClusterMeanSilhouette(zMCluster *mc)
{
  zClusterListCell *cp, *ccp;
  zVecAddrListCell *pp;
  double a, b, b_tmp, score = 0;
  int i, n = 0;

  zListForEach( zMClusterClusterList(mc), cp ){
    zFree( cp->data._sil );
    if( !( cp->data._sil = zAlloc( double, zListSize(zClusterSampleList(&cp->data)) ) ) ){
      ZALLOCERROR();
      return NAN;
    }
    i = 0;
    zListForEach( zClusterSampleList(&cp->data), pp ){
      a = _zClusterDistAve( &mc->method, zClusterSampleList(&cp->data), pp->data ); /* intra-cluster distance */
      b = HUGE_VAL;
      zListForEach( zMClusterClusterList(mc), ccp ){
        if( ccp == cp ) continue;
        if( ( b_tmp = _zClusterDistAve( &mc->method, zClusterSampleList(&ccp->data), pp->data ) ) < b )
          b = b_tmp; /* inter-cluster distance */
      }
      score += cp->data._sil[i++] = ( b - a ) / zMax( a, b ); /* silhouette of the sample */
    }
    zDataSort( cp->data._sil, zListSize(zClusterSampleList(&cp->data)) );
    n += zListSize(zClusterSampleList(&cp->data));
  }
  return score / n;
}

/* print silhouettes of a vector cluster to a file. */
static bool _zClusterSilhouetteFPrint(FILE *fp, zCluster *c, int offset)
{
  int i;

  if( !c->_sil ){
    ZRUNWARN( ZM_WARN_MCA_NOSILHOUETTE );
    return false;
  }
  for( i=0; i<zListSize(zClusterSampleList(c)); i++ ){
    fprintf( fp, "0, %d\n", offset + i );
    fprintf( fp, "%f, %d\n\n", c->_sil[i], offset + i );
  }
  return true;
}

/* print silhouettes of a set of vector clusters to files. */
bool zMClusterSilhouettePrintFile(zMCluster *mc, const char *basename)
{
  zClusterListCell *cp;
  int i = 0, offset = 0;
  char filename[BUFSIZ];
  FILE *fp;

  zListForEach( zMClusterClusterList(mc), cp ){
    sprintf( filename, "%s%d", basename, i++ );
    if( !( fp = fopen( filename, "w" ) ) ){
      ZOPENERROR( filename );
      return false;
    }
    if( !_zClusterSilhouetteFPrint( fp, &cp->data, offset ) ) return false;
    offset += zListSize(zClusterSampleList(&cp->data));
  }
  return true;
}

/* ********************************************************** */
/*! \brief clustering based on X-means
 *//* ******************************************************* */

#define Z_XMEANS_MINSIZE 10

/* merge subclusters with original. */
static void _zMClusterMerge(zMCluster *mc, zClusterListCell *vc, zMCluster *submc)
{
  zClusterListCell *prev, *next;

  prev = zListCellPrev(vc);
  next = zListCellNext(vc);
  zClusterDestroy( &vc->data );
  free( vc );
  zListCellBind( prev, zListTail(zMClusterClusterList(submc)) );
  zListCellBind( zListHead(zMClusterClusterList(submc)), next );
  zListSize(zMClusterClusterList(mc)) += zListSize(zMClusterClusterList(submc)) - 1;
  zListInit( zMClusterClusterList(submc) );
  zMClusterDestroy( submc );
}

/* initialize clusters for X-means. */
static bool _zMClusterXMeansInit(zMCluster *mc, zVecAddrList *points)
{
  zClusterListCell *vc;
  zVecListCell *pc;

  vc = zListHead(zMClusterClusterList(mc));
  zListForEach( points, pc )
    if( !zVecAddrListInsertHead( zClusterSampleList(&vc->data), pc->data ) ) return false;
  zMClusterCoreF( mc, zClusterSampleList(&vc->data), vc->data.core );
  _zMClusterVar( mc, &vc->data );
  return true;
}

/* find maximum error in a cluster from the mean value. */
static double _zClusterMaxError(zCluster *c, zClusterMethod *method)
{
  zVecListCell *vc;
  double dmax = 0, d;

  zListForEach( zClusterSampleList(c), vc ){
    zClusterMethodErrorF( method, vc->data, c->core );
    if( ( d = zVecNorm( method->error ) ) > dmax ) dmax = d;
  }
  return dmax;
}

/* check if the subclusters better classify samples than the original based on hyperdensity of clusters. */
static bool _zMClusterXMeansCheckDensity(zClusterMethod *cm, zCluster *c, zCluster *c1, zCluster *c2)
{
  double rm, rm1, rm2;

  /* original cluster */
  rm = pow( _zClusterMaxError( c, cm ), zVecSizeNC(cm->error) );
  /* bidivided clusters */
  rm1 = pow( _zClusterMaxError( c1, cm ), zVecSizeNC(cm->error) );
  rm2 = pow( _zClusterMaxError( c2, cm ), zVecSizeNC(cm->error) );
  return rm > rm1 + rm2 ? true : false;
}

/* probability density function based on Gaussian distribution. */
static double _zClusterXMeansPDF(zClusterMethod *cm, zCluster *c, zVec x)
{
  zClusterMethodErrorF( cm, x, c->core );
  return exp( -0.5*zVecSqrNorm(cm->error)/c->var ) / sqrt(c->var);
}

/* check if the subclusters better classify samples than the original based on BIC. */
static bool _zMClusterXMeansCheckBIC(zClusterMethod *cm, zCluster *c, zCluster *c1, zCluster *c2)
{
  zVecListCell *vc;
  int n;
  double ls = 0, bic1, bic2;

  n = zListSize(zClusterSampleList(c));
  zListForEach( zClusterSampleList(c1), vc ){
    ls += log( _zClusterXMeansPDF( cm, c1, vc->data )
             + _zClusterXMeansPDF( cm, c2, vc->data ) );
  }
  zListForEach( zClusterSampleList(c2), vc ){
    ls += log( _zClusterXMeansPDF( cm, c1, vc->data )
             + _zClusterXMeansPDF( cm, c2, vc->data ) );
  }
  bic1 = n * ( log(c->var) + 1 );
  bic2 = 2*( n*log(2) - ls ) + cm->core_size * log(n);
  return bic1 < bic2 ? true : false;
}

/* test if the subclusters better classify samples than the original. */
static bool _zMClusterXMeansTest(zMCluster *mc, zCluster *c, zMCluster *submc, bool (* testfunc)(zClusterMethod*,zCluster*,zCluster*,zCluster*))
{
  zCluster *c1, *c2;

  c1 = &zListHead(zMClusterClusterList(submc))->data;
  c2 = &zListTail(zMClusterClusterList(submc))->data;
  if( 10 * zListSize(zClusterSampleList(c1)) < zListSize(zClusterSampleList(c2)) ||
      10 * zListSize(zClusterSampleList(c2)) < zListSize(zClusterSampleList(c1)) )
    return false; /* avoid too crisp cluster */
  return testfunc( &mc->method, c, c1, c2 );
}

/* internal function for recursive all of X-means method. */
static int _zMClusterXMeans(zMCluster *mc, zClusterListCell *cc, bool (* testfunc)(zClusterMethod*,zCluster*,zCluster*,zCluster*))
{
  zMCluster submc;
  int iter, iter1, iter2;

  /* ignore a too small cluster */
  if( zListSize(zClusterSampleList(&cc->data)) < Z_XMEANS_MINSIZE ) return 0;
  if( !zMClusterInit( &submc, mc->method.core_size ) ||
      !zMClusterMethodCopy( mc, &submc ) ) return 0;

  iter = zMClusterKMeans( &submc, zClusterSampleList(&cc->data), 2 );
  if( _zMClusterXMeansTest( mc, &cc->data, &submc, testfunc ) ){
    zMClusterDestroy( &submc );
    return 0;
  }
  iter1 = _zMClusterXMeans( &submc, zListHead(zMClusterClusterList(&submc)), testfunc );
  iter2 = _zMClusterXMeans( &submc, zListTail(zMClusterClusterList(&submc)), testfunc );
  _zMClusterMerge( mc, cc, &submc );
  return iter + iter1 + iter2;
}

/* cluster vectors by X-means method. */
int zMClusterXMeans(zMCluster *mc, zVecAddrList *points, bool (* testfunc)(zClusterMethod*,zCluster*,zCluster*,zCluster*))
{
  if( !zMClusterAlloc( mc, 1 ) ) return -1;
  if( !_zMClusterXMeansInit( mc, points ) ) return -1;
  return _zMClusterXMeans( mc, zListHead(zMClusterClusterList(mc)), testfunc );
}

/* cluster vectors by X-means method based on hyperdensity. */
int zMClusterXMeansDensity(zMCluster *mc, zVecAddrList *points)
{
  return zMClusterXMeans( mc, points, _zMClusterXMeansCheckDensity );
}

/* cluster vectors by X-means method based on BIC. */
int zMClusterXMeansBIC(zMCluster *mc, zVecAddrList *points)
{
  return zMClusterXMeans( mc, points, _zMClusterXMeansCheckBIC );
}
