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
zCluster *zClusterCreate(zCluster *c, int meansize)
{
  zListInit( &c->vl );
  if( !( c->mean = zVecAlloc( meansize ) ) ){
    ZALLOCERROR();
    zClusterDestroy( c );
    return NULL;
  }
  return c;
}

/* destroy a vector cluster. */
void zClusterDestroy(zCluster *c)
{
  zVecAddrListDestroy( &c->vl );
  zVecFree( c->mean );
}

/* print a vector cluster to a file. */
void zClusterFPrint(FILE *fp, zCluster *c)
{
  register int i = 0;
  zVecListCell *vc;

  fprintf( fp, "%d members\n", zListSize(&c->vl) );
  zListForEach( &c->vl, vc ){
    fprintf( fp, " %d: ", i++ );
    zVecFPrint( fp, vc->data );
  }
  fprintf( fp, " mean: " );
  zVecFPrint( fp, c->mean );
}

/* print data of vectors of a vector cluster to a file. */
void zClusterDataFPrint(FILE *fp, zCluster *c)
{
  zVecListFPrint( fp, &c->vl );
}

/* ********************************************************** */
/* methods for mean and error computation
 * ********************************************************** */

/* a default error function for clustering. */
static zVec _zClusterErrorDefault(zVec p, zVec mean, void *dummy, zVec err)
{
  return zVecSub( p, mean, err );
}

/* a default mean computation function for clustering. */
static zVec _zClusterMeanDefault(zVecAddrList *pl, void *dummy, zVec mean)
{
  return zVecListMean( pl, mean );
}

/* a default loaded mean computation function for clustering. */
static zVec _zClusterMeanLoadedDefault(zVecAddrList *pl, double load[], double n, void *dummy, zVec mean)
{
  zVecListCell *vc;
  register int i = 0;

  zVecZero( mean );
  zListForEach( pl, vc ){
    zVecCatDRC( mean, load[i], vc->data );
    i++;
  }
  return zVecDivDRC( mean, n );
}

/* initialize methods for clustering. */
void zClusterMethodInit(zClusterMethod *met)
{
  met->_meansize = 0;
  met->_mean_fp = NULL;
  met->_errorsize = 0;
  met->_error_fp = NULL;
  met->_err = NULL;
  met->_mean_l_fp = NULL;
  met->_err_l = NULL;
}

/* create methods for clustering. */
zClusterMethod *zClusterMethodCreate(zClusterMethod *met, int meansize, zVec (* mean_fp)(zVecAddrList*,void*,zVec), int errorsize, zVec (* error_fp)(zVec,zVec,void*,zVec))
{
  zClusterMethodInit( met );
  met->_meansize = meansize;
  met->_mean_fp = mean_fp ? mean_fp : _zClusterMeanDefault;
  met->_errorsize = errorsize;
  met->_error_fp = error_fp ? error_fp : _zClusterErrorDefault;
  if( !( met->_err = zVecAlloc( errorsize ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  return met;
}

/* assign a loaded mean computation function. */
zClusterMethod *zClusterMethodLoadedCreate(zClusterMethod *met, zVec (* mean_l_fp)(zVecAddrList*,double[],double,void*,zVec))
{
  if( !met->_mean_fp ){
    ZRUNERROR( ZM_ERR_MCA_NOMEAN );
    return NULL;
  }
  met->_mean_l_fp = mean_l_fp ? mean_l_fp : _zClusterMeanLoadedDefault;
  if( !( met->_err_l = zVecAlloc( met->_errorsize ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  return met;
}

/* destroy methods for clustering. */
void zClusterMethodDestroy(zClusterMethod *met)
{
  zVecFree( met->_err );
  zVecFree( met->_err_l );
}

/* ********************************************************** */
/* multiple vecter clusters class.
 * ********************************************************** */

/* initialize multiple vector clusters. */
zMCluster *zMClusterInit(zMCluster *mc, int meansize, zVec (* mean_fp)(zVecAddrList*,void*,zVec), int errorsize, zVec (* error_fp)(zVec,zVec,void*,zVec))
{
  zListInit( &mc->cl );
  return zClusterMethodCreate( &mc->met, meansize, mean_fp, errorsize, error_fp ) != NULL ? mc : NULL;
}

/* allocate multiple vector clusters. */
zMCluster *zMClusterAlloc(zMCluster *mc, int n)
{
  zClusterListCell *cc;
  register int i;

  zListInit( &mc->cl );
  for( i=0; i<n; i++ ){
    if( !( cc = zAlloc( zClusterListCell, 1 ) ) ){
      ZALLOCERROR();
      goto ERR;
    }
    if( !zClusterCreate( &cc->data, mc->met._meansize ) ){
      ZALLOCERROR();
      free( cc );
      goto ERR;
    }
    zListInsertHead( &mc->cl, cc );
  }
  return mc;

 ERR:
  zMClusterDestroy( mc );
  return NULL;
}

/* destroy multiple vector clusters. */
void zMClusterDestroy(zMCluster *mc)
{
  zClusterListCell *cc;

  while( !zListIsEmpty(&mc->cl) ){
    zListDeleteHead( &mc->cl, &cc );
    zClusterDestroy( &cc->data );
    free( cc );
  }
  zClusterMethodDestroy( &mc->met );
}

/* print multiple vector clusters to a file. */
void zMClusterFPrint(FILE *fp, zMCluster *mc)
{
  register int i = 0;
  zClusterListCell *cc;

  zListForEach( &mc->cl, cc ){
    fprintf( fp, "#cluster[%d] : ", i++ );
    zClusterFPrint( fp, &cc->data );
  }
}

/* print data of multiple vector clusters to files. */
void zMClusterDataFPrint(FILE *fp[], zMCluster *mc)
{
  zClusterListCell *cc;
  int i = 0;

  zListForEach( &mc->cl, cc )
    zClusterDataFPrint( fp[i++], &cc->data );
}

/* print means of each cluster of multiple vector clusters to files. */
void zMClusterMeanFPrint(FILE *fp[], zMCluster *mc)
{
  zClusterListCell *cc;
  int i = 0;

  zListForEach( &mc->cl, cc )
    zVecDataFPrint( fp[i++], cc->data.mean );
}

/* ********************************************************** */
/* clustering based on K-means family
 * ********************************************************** */

/* ********************************************************** */
/* K-means
 * ********************************************************** */

#if DEBUG
/* for debug */
static void _zMClusterKMeansDataPrintFile(zMCluster *mc, int step)
{
  FILE *fp;
  char filename[BUFSIZ];
  zClusterListCell *cc;
  int i;

  i = 0;
  zListForEach( &mc->cl, cc ){
    sprintf( filename, "%d_%d", step, i );
    fp = fopen( filename, "w" );
    zClusterDataFPrint( fp, &cc->data );
    fclose( fp );
    sprintf( filename, "%d_%dm", step, i );
    fp = fopen( filename, "w" );
    zVecDataFPrint( fp, cc->data.mean );
    fclose( fp );
  }
}
#endif /* DEBUG */

/* compute variance of a vector cluster. */
static double _zMClusterVar(zCluster *c, zMCluster *mc, void *err_util)
{
  zVecListCell *pc;

  c->var = 0;
  zListForEach( &c->vl, pc ){
    mc->met._error_fp( pc->data, c->mean, err_util, mc->met._err );
    c->var += zVecSqrNorm( mc->met._err );
  }
  return ( c->var /= zListSize(&c->vl) );
}

/* arrange initial candidates of clusters. */
static void _zMClusterKMeansInit(zMCluster *mc, zVecAddrList *points)
{
  zClusterListCell *cc, *ccnn;
  zVecListCell *pc, *pcmax = NULL;
  zVec mean;
  double d, dmax, dmin;

  mean = zVecAlloc( mc->met._meansize );
  zVecListMean( points, mean );
  /* first core point */
  dmax = 0;
  zListForEach( points, pc ){
    if( ( d = zVecDist( pc->data, mean ) ) > dmax ){
      dmax = d;
      pcmax = pc;
    }
  }
  cc = zListTail( &mc->cl );
  if( !zVecAddrListInsertTail( &cc->data.vl, pcmax->data ) ) return;
  zVecFree( mean );
  /* second - kth core point */
  for( cc=zListCellNext(cc); cc!=zListRoot(&mc->cl); cc=zListCellNext(cc) ){
    dmax = 0;
    pcmax = NULL;
    zListForEach( points, pc ){
      dmin = HUGE_VAL;
      for( ccnn=zListTail(&mc->cl); ccnn!=cc; ccnn=zListCellNext(ccnn) )
        if( ( d = zVecDist( pc->data, zListTail(&ccnn->data.vl)->data ) ) < dmin )
          dmin = d;
      if( dmin > dmax ){
        dmax = dmin;
        pcmax = pc;
      }
    }
    if( !zVecAddrListInsertTail( &cc->data.vl, pcmax->data ) ) return;
  }
  /* assign all points to k points */
  zListForEach( points, pc ){
    dmin = HUGE_VAL;
    ccnn = NULL;
    zListForEach( &mc->cl, cc ){
      if( pc->data == zListTail(&cc->data.vl)->data ) goto CONTINUE;
      if( ( d = zVecDist( pc->data, zListTail(&cc->data.vl)->data ) ) < dmin ){
        dmin = d;
        ccnn = cc;
      }
    }
    if( !zVecAddrListInsertHead( &ccnn->data.vl, pc->data ) ) return;
    CONTINUE: ;
  }
}

/* recluster tentative clusters. */
static int _zMClusterKMeansRecluster(zMCluster *mc, void *mean_util, void *err_util)
{
  zClusterListCell *cc1, *cc2, *cc;
  zVecListCell *pc, *pc_prev;
  bool ismoved;
  double d, dmin;
  int iter = 0;
  register int i;

  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    /* compute means */
    zListForEach( &mc->cl, cc1 ){
      if( !zListIsEmpty( &cc1->data.vl ) ){
        mc->met._mean_fp( &cc1->data.vl, mean_util, cc1->data.mean );
        _zMClusterVar( &cc1->data, mc, err_util );
      }
    }
    ismoved = false;
    zListForEach( &mc->cl, cc1 ){
      /* for each tentative cluster */
      zListForEach( &cc1->data.vl, pc ){ /* for each point */
        dmin = HUGE_VAL;
        cc = cc1;
        zListForEach( &mc->cl, cc2 ){ /* find the nearest mean point */
          if( zListIsEmpty( &cc2->data.vl ) ) continue;
          mc->met._error_fp( pc->data, cc2->data.mean, err_util, mc->met._err );
          if( ( d = zVecSqrNorm( mc->met._err ) ) < dmin ){
            dmin = d;
            cc = cc2;
          }
        }
        if( cc != cc1 ){ /* move the point into nearer cluster */
          ismoved = true;
          pc_prev = zListCellPrev(pc);
          zListPurge( &cc1->data.vl, pc );
          zListInsertHead( &cc->data.vl, pc );
          pc = pc_prev;
        }
      }
    }
    if( ismoved == false ) return i; /* clusters settled */
  }
  ZITERWARN( iter );
  return iter;
}

/* clustering of vectors by K-means. */
int zMClusterKMeans(zMCluster *mc, zVecAddrList *points, int k, void *mean_util, void *err_util)
{
  if( !zMClusterAlloc( mc, k ) ) return -1;
  _zMClusterKMeansInit( mc, points );
  return _zMClusterKMeansRecluster( mc, mean_util, err_util );
}

/* ********************************************************** */
/* X-means based on filling rate / Baysian Information Criterion
 * ********************************************************** */

#define Z_XMEANS_MINSIZE 10

/* merge subclusters with original. */
static void _zMClusterMerge(zMCluster *mc, zClusterListCell *vc, zMCluster *submc)
{
  zClusterListCell *prev, *next;

  prev = zListCellPrev(vc);
  next = zListCellNext(vc);
  zClusterDestroy( &vc->data );
  free( vc );
  zListCellBind( prev, zListTail(&submc->cl) );
  zListCellBind( zListHead(&submc->cl), next );
  zListSize(&mc->cl) += zListSize(&submc->cl) - 1;
  zListInit( &submc->cl );
  zMClusterDestroy( submc );
}

/* initialize clusters for X-means. */
static void _zMClusterXMeansInit(zMCluster *mc, zVecAddrList *points, void *mean_util, void *err_util)
{
  zClusterListCell *vc;
  zVecListCell *pc;

  vc = zListHead(&mc->cl);
  zListForEach( points, pc )
    if( !zVecAddrListInsertHead( &vc->data.vl, pc->data ) ) break;
  mc->met._mean_fp( &vc->data.vl, mean_util, vc->data.mean );
  _zMClusterVar( &vc->data, mc, err_util );
}

/* find maximum error in a cluster from the mean value. */
static double _zClusterMaxError(zCluster *c, zVec (* error_fp)(zVec,zVec,void*,zVec), void *err_util, zVec err)
{
  zVecListCell *vc;
  double dmax = 0, d;

  zListForEach( &c->vl, vc ){
    error_fp( vc->data, c->mean, err_util, err );
    if( ( d = zVecNorm( err ) ) > dmax ) dmax = d;
  }
  return dmax;
}

/* test if the subclusters better classify the original cluster. */
static bool _zMClusterXMeansTest(zMCluster *mc, zCluster *org, zMCluster *submc, void *err_util)
{
  zCluster *c1, *c2;
  double fr1, fr2;
  double r, r1, r2;

  c1 = &zListHead(&submc->cl)->data;
  c2 = &zListTail(&submc->cl)->data;
  if( zListSize(&c1->vl) < 0.1 * zListSize(&c2->vl) ||
      zListSize(&c2->vl) < 0.1 * zListSize(&c1->vl) )
    return false; /* avoid too crisp cluster */

  /* original cluster */
  r = _zClusterMaxError( org, mc->met._error_fp, err_util, mc->met._err );
  fr1 = log(zListSize(&org->vl)) - mc->met._errorsize*log(r);
  /* bidivided clusters */
  r1 = _zClusterMaxError( c1, mc->met._error_fp, err_util, mc->met._err );
  r2 = _zClusterMaxError( c2, mc->met._error_fp, err_util, mc->met._err );
  fr2 = log( zListSize(&org->vl) )
      - log( pow(r1,mc->met._errorsize) + pow(r2,mc->met._errorsize) );
  return fr1 < fr2 ? true : false;
}

/* cluster vectors by X-means in a recursive way based on filling rate. */
static int _zMClusterXMeans(zMCluster *mc, zClusterListCell *cc, void *mean_util, void *err_util)
{
  zMCluster submc;
  int iter, iter1, iter2;

  /* avoid a too small / sparse cluster */
  if( zListSize(&cc->data.vl) < Z_XMEANS_MINSIZE ) return 0;

  if( !zMClusterInit( &submc, mc->met._meansize, mc->met._mean_fp, mc->met._errorsize, mc->met._error_fp ) )
    return 0;

  iter = zMClusterKMeans( &submc, &cc->data.vl, 2, mean_util, err_util );
  if( !_zMClusterXMeansTest( mc, &cc->data, &submc, err_util ) ){
    zMClusterDestroy( &submc );
    return 0;
  }
  iter1 = _zMClusterXMeans( &submc, zListHead(&submc.cl), mean_util, err_util );
  iter2 = _zMClusterXMeans( &submc, zListTail(&submc.cl), mean_util, err_util );
  _zMClusterMerge( mc, cc, &submc );
  return iter + iter1 + iter2;
}

/* cluster vectors by X-means based on filling rate. */
int zMClusterXMeans(zMCluster *mc, zVecAddrList *points, void *mean_util, void *err_util)
{
  if( !zMClusterAlloc( mc, 1 ) ) return -1;
  _zMClusterXMeansInit( mc, points, mean_util, err_util );
  return _zMClusterXMeans( mc, zListHead(&mc->cl), mean_util, err_util );
}

/* probability density function based on Gaussian distribution. */
static double _zClusterXMeansPDF(zCluster *c, zMCluster *mc, void *err_util, zVec x)
{
  mc->met._error_fp( x, c->mean, err_util, mc->met._err );
  return exp( -0.5*zVecSqrNorm(mc->met._err)/c->var ) / sqrt(c->var);
}

/* test if the subclusters better classify the original based on BIC. */
static bool _zMClusterXMeansBICTest(zMCluster *mc, zCluster *c, zMCluster *submc, void *err_util)
{
  zCluster *c1, *c2;
  zVecListCell *vc;
  int q, n;
  double ls = 0, bic1, bic2;

  q = mc->met._meansize + 1;
  n = zListSize(&c->vl);
  c1 = &zListHead(&submc->cl)->data;
  c2 = &zListTail(&submc->cl)->data;
  zListForEach( &c1->vl, vc ){
    ls += log( _zClusterXMeansPDF( c1, submc, err_util, vc->data )
             + _zClusterXMeansPDF( c2, submc, err_util, vc->data ) );
  }
  zListForEach( &c2->vl, vc ){
    ls += log( _zClusterXMeansPDF( c1, submc, err_util, vc->data )
             + _zClusterXMeansPDF( c2, submc, err_util, vc->data ) );
  }
  bic1 = n * ( log(c->var) + 1 );
  bic2 = 2*( n*log(2) - ls ) + q*log(n);
  return bic1 < bic2 ? true : false;
}

/* cluster vectors by X-means based on BIC. */
static int _zMClusterXMeansBIC(zMCluster *mc, zClusterListCell *cc, void *mean_util, void *err_util)
{
  zMCluster submc;
  int iter, iter1, iter2;

  /* ignore a too small cluster */
  if( zListSize(&cc->data.vl) < Z_XMEANS_MINSIZE ) return 0;
  if( !zMClusterInit( &submc, mc->met._meansize, mc->met._mean_fp, mc->met._errorsize, mc->met._error_fp ) )
    return 0;
  iter = zMClusterKMeans( &submc, &cc->data.vl, 2, mean_util, err_util );
  if( _zMClusterXMeansBICTest( mc, &cc->data, &submc, err_util ) ){
    zMClusterDestroy( &submc );
    return 0;
  }
  iter1 = _zMClusterXMeansBIC( &submc, zListHead(&submc.cl), mean_util, err_util );
  iter2 = _zMClusterXMeansBIC( &submc, zListTail(&submc.cl), mean_util, err_util );
  _zMClusterMerge( mc, cc, &submc );
  return iter + iter1 + iter2;
}

/* cluster vectors based on BIC. */
int zMClusterXMeansBIC(zMCluster *mc, zVecAddrList *points, void *mean_util, void *err_util)
{
  if( !zMClusterAlloc( mc, 1 ) ) return -1;
  _zMClusterXMeansInit( mc, points, mean_util, err_util );
  return _zMClusterXMeansBIC( mc, zListHead(&mc->cl), mean_util, err_util );
}
