/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_mca_cluster - multiple classification analysis:
 * clustering.
 */

#include <zm/zm_mca.h>

/* ********************************************************** */
/* vector cluster class.
 * ********************************************************** */

/* zClusterCreate
 * - create and initialize a vector cluster.
 */
zCluster *zClusterCreate(zCluster *c, int meansize)
{
  zListInit( &c->vl );
  if( ( c->mean = zVecAlloc( meansize ) ) == NULL ){
    ZALLOCERROR();
    zClusterDestroy( c );
    return NULL;
  }
  return c;
}

/* zClusterDestroy
 * - destroy a vector cluster.
 */
void zClusterDestroy(zCluster *c)
{
  zVecListDestroy( &c->vl, false );
  zVecFree( c->mean );
}

/* zClusterFWrite
 * - output a vector cluster to a file.
 */
void zClusterFWrite(FILE *fp, zCluster *c)
{
  register int i = 0;
  zVecListCell *vc;

  fprintf( fp, "%d members\n", zListNum(&c->vl) );
  zListForEach( &c->vl, vc ){
    fprintf( fp, " %d: ", i++ );
    zVecFWrite( fp, vc->data );
  }
  fprintf( fp, " mean: " );
  zVecFWrite( fp, c->mean );
}

/* zClusterFWrite
 * - output data of vectors of a vector cluster to a file.
 */
void zClusterDataFWrite(FILE *fp, zCluster *c)
{
  zVecListFWrite( fp, &c->vl );
}

/* ********************************************************** */
/* methods for mean and error computation
 * ********************************************************** */

static zVec _zClusterErrorDefault(zVec p, zVec mean, void *dummy, zVec err);
static zVec _zClusterMeanDefault(zVecList *pl, void *dummy, zVec mean);
static zVec _zClusterMeanLoadedDefault(zVecList *pl, double load[], double n, void *dummy, zVec mean);

/* (static)
 * _zClusterErrorDefault
 * - a default error function for clustering.
 */
zVec _zClusterErrorDefault(zVec p, zVec mean, void *dummy, zVec err)
{
  return zVecSub( p, mean, err );
}

/* (static)
 * _zClusterMeanDefault
 * - a default mean computation function for clustering.
 */
zVec _zClusterMeanDefault(zVecList *pl, void *dummy, zVec mean)
{
  return zVecListAve( pl, mean );
}

/* (static)
 * _zClusterMeanLoadedDefault
 * - a default loaded mean computation function for clustering.
 */
zVec _zClusterMeanLoadedDefault(zVecList *pl, double load[], double n, void *dummy, zVec mean)
{
  zVecListCell *vc;
  register int i = 0;

  zVecClear( mean );
  zListForEach( pl, vc ){
    zVecCatDRC( mean, load[i], vc->data );
    i++;
  }
  return zVecDivDRC( mean, n );
}

/* zClusterMethodInit
 * - initialize methods for clustering.
 */
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

/* zClusterMethodCreate
 * - create methods for clustering.
 */
zClusterMethod *zClusterMethodCreate(zClusterMethod *met, int meansize, zVec (* mean_fp)(zVecList*,void*,zVec), int errorsize, zVec (* error_fp)(zVec,zVec,void*,zVec))
{
  zClusterMethodInit( met );
  met->_meansize = meansize;
  met->_mean_fp = mean_fp ? mean_fp : _zClusterMeanDefault;
  met->_errorsize = errorsize;
  met->_error_fp = error_fp ? error_fp : _zClusterErrorDefault;
  if( ( met->_err = zVecAlloc( errorsize ) ) == NULL ){
    ZALLOCERROR();
    return NULL;
  }
  return met;
}

/* zClusterMethodCreateMeanLoaded
 * - assign a loaded mean computation function.
 */
zClusterMethod *zClusterMethodLoadedCreate(zClusterMethod *met, zVec (* mean_l_fp)(zVecList*,double[],double,void*,zVec))
{
  if( met->_mean_fp == NULL ){
    ZRUNERROR( ZM_ERR_MCA_NOMEAN );
    return NULL;
  }
  met->_mean_l_fp = mean_l_fp ? mean_l_fp : _zClusterMeanLoadedDefault;
  if( ( met->_err_l = zVecAlloc( met->_errorsize ) ) == NULL ){
    ZALLOCERROR();
    return NULL;
  }
  return met;
}

/* zClusterMethodDestroy
 * - destroy methods for clustering.
 */
void zClusterMethodDestroy(zClusterMethod *met)
{
  zVecFree( met->_err );
  zVecFree( met->_err_l );
}

/* ********************************************************** */
/* multiple vecter clusters class.
 * ********************************************************** */

/* zMClusterInit
 * - initialize multiple vector clusters.
 */
zMCluster *zMClusterInit(zMCluster *mc, int meansize, zVec (* mean_fp)(zVecList*,void*,zVec), int errorsize, zVec (* error_fp)(zVec,zVec,void*,zVec))
{
  zListInit( &mc->cl );
  return zClusterMethodCreate( &mc->met, meansize, mean_fp, errorsize, error_fp ) != NULL ? mc : NULL;
}

/* zMClusterAlloc
 * - allocate multiple vector clusters.
 */
zMCluster *zMClusterAlloc(zMCluster *mc, int n)
{
  zClusterListCell *cc;
  register int i;

  zListInit( &mc->cl );
  for( i=0; i<n; i++ ){
    if( ( cc = zAlloc( zClusterListCell, 1 ) ) == NULL ){
      ZALLOCERROR();
      goto ERR;
    }
    if( zClusterCreate( &cc->data, mc->met._meansize ) == NULL ){
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

/* zMClusterDestroy
 * - destroy multiple vector clusters.
 */
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

/* zMClusterFWrite
 * - output multiple vector clusters to a file.
 */
void zMClusterFWrite(FILE *fp, zMCluster *mc)
{
  register int i = 0;
  zClusterListCell *cc;

  zListForEach( &mc->cl, cc ){
    fprintf( fp, "#cluster[%d] : ", i++ );
    zClusterFWrite( fp, &cc->data );
  }
}

/* zMClusterDataFWrite
 * - output data of multiple vector clusters to files.
 */
void zMClusterDataFWrite(FILE *fp[], zMCluster *mc)
{
  zClusterListCell *cc;
  int i = 0;

  zListForEach( &mc->cl, cc )
    zClusterDataFWrite( fp[i++], &cc->data );
}

/* zMClusterMeanFWrite
 * - output means of each cluster of multiple vector clusters to files.
 */
void zMClusterMeanFWrite(FILE *fp[], zMCluster *mc)
{
  zClusterListCell *cc;
  int i = 0;

  zListForEach( &mc->cl, cc )
    zVecDataFWrite( fp[i++], cc->data.mean );
}

/* ********************************************************** */
/* clustering based on K-means family
 * ********************************************************** */

/* ********************************************************** */
/* K-means
 * ********************************************************** */

#if DEBUG
/* for debug */
static void _zMClusterKMeansDataWriteFile(zMCluster *mc, int step);
void _zMClusterKMeansDataWriteFile(zMCluster *mc, int step)
{
  FILE *fp;
  char filename[BUFSIZ];
  zClusterListCell *cc;
  int i;

  i = 0;
  zListForEach( &mc->cl, cc ){
    sprintf( filename, "%d_%d", step, i );
    fp = fopen( filename, "w" );
    zClusterDataFWrite( fp, &cc->data );
    fclose( fp );
    sprintf( filename, "%d_%dm", step, i );
    fp = fopen( filename, "w" );
    zVecDataFWrite( fp, cc->data.mean );
    fclose( fp );
  }
}
#endif /* DEBUG */

static double _zMClusterVar(zCluster *c, zMCluster *mc, void *err_util);
static void _zMClusterKMeansInit(zMCluster *mc, zVecList *points);
static int _zMClusterKMeansRecluster(zMCluster *mc, void *mean_util, void *err_util);

/* (static)
 * _zMClusterVar
 * - compute variance of a vector cluster.
 */
double _zMClusterVar(zCluster *c, zMCluster *mc, void *err_util)
{
  zVecListCell *pc;

  c->var = 0;
  zListForEach( &c->vl, pc ){
    mc->met._error_fp( pc->data, c->mean, err_util, mc->met._err );
    c->var += zVecSqrNorm( mc->met._err );
  }
  return ( c->var /= zListNum(&c->vl) );
}

/* (static)
 * _zMClusterKMeansInit
 * - arrange initial candidates of clusters.
 */
void _zMClusterKMeansInit(zMCluster *mc, zVecList *points)
{
  zClusterListCell *cc, *ccnn;
  zVecListCell *pc, *pcmax = NULL;
  zVec mean;
  double d, dmax, dmin;

  mean = zVecAlloc( mc->met._meansize );
  zVecListAve( points, mean );
  /* first core point */
  dmax = 0;
  zListForEach( points, pc ){
    if( ( d = zVecDist( pc->data, mean ) ) > dmax ){
      dmax = d;
      pcmax = pc;
    }
  }
  cc = zListTail( &mc->cl );
  if( zVecListInsertTail( &cc->data.vl, pcmax->data, false ) == NULL )
    return;
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
    if( zVecListInsertTail( &cc->data.vl, pcmax->data, false ) == NULL )
      return;
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
    if( zVecListInsertHead( &ccnn->data.vl, pc->data, false ) == NULL )
      return;
    CONTINUE: ;
  }
}

/* (static)
 * _zMClusterKMeansRecluster
 * - recluster tentative clusters.
 */
int _zMClusterKMeansRecluster(zMCluster *mc, void *mean_util, void *err_util)
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

/* zMClusterKMeans
 * - clustering of vectors by K-means.
 */
int zMClusterKMeans(zMCluster *mc, zVecList *points, int k, void *mean_util, void *err_util)
{
  if( zMClusterAlloc( mc, k ) == NULL ) return -1;
  _zMClusterKMeansInit( mc, points );
  return _zMClusterKMeansRecluster( mc, mean_util, err_util );
}

/* ********************************************************** */
/* X-means based on filling rate / Baysian Information Criterion
 * ********************************************************** */

#define Z_XMEANS_MINSIZE 10

static void _zMClusterMerge(zMCluster *mc, zClusterListCell *vc, zMCluster *submc);
static void _zMClusterXMeansInit(zMCluster *mc, zVecList *points, void *mean_util, void *err_util);

static double _zClusterMaxError(zCluster *c, zVec (* error_fp)(zVec,zVec,void*,zVec), void *err_util, zVec err);
static bool _zMClusterXMeansTest(zMCluster *mc, zCluster *org, zMCluster *submc, void *err_util);
static int _zMClusterXMeans(zMCluster *mc, zClusterListCell *cc, void *mean_util, void *err_util);

static double _zClusterXMeansPDF(zCluster *c, zMCluster *mc, void *err_util, zVec x);
static bool _zMClusterXMeansBICTest(zMCluster *mc, zCluster *c, zMCluster *submc, void *err_util);
static int _zMClusterXMeansBIC(zMCluster *mc, zClusterListCell *cc, void *mean_util, void *err_util);

/* (static)
 * _zMClusterMerge
 * - merge subclusters with original.
 */
void _zMClusterMerge(zMCluster *mc, zClusterListCell *vc, zMCluster *submc)
{
  zClusterListCell *prev, *next;

  prev = zListCellPrev(vc);
  next = zListCellNext(vc);
  zClusterDestroy( &vc->data );
  free( vc );
  zListCellBind( prev, zListTail(&submc->cl) );
  zListCellBind( zListHead(&submc->cl), next );
  zListNum(&mc->cl) += zListNum(&submc->cl) - 1;
  zListInit( &submc->cl );
  zMClusterDestroy( submc );
}

/* (static)
 * _zMClusterXMeansInit
 * - initialize clusters for X-means.
 */
void _zMClusterXMeansInit(zMCluster *mc, zVecList *points, void *mean_util, void *err_util)
{
  zClusterListCell *vc;
  zVecListCell *pc;

  vc = zListHead(&mc->cl);
  zListForEach( points, pc ){
    if( zVecListInsertHead( &vc->data.vl, pc->data, false ) == NULL )
      break;
  }
  mc->met._mean_fp( &vc->data.vl, mean_util, vc->data.mean );
  _zMClusterVar( &vc->data, mc, err_util );
}

/* (static)
 * _zClusterMaxError
 * - find maximum error in a cluster from mean.
 */
double _zClusterMaxError(zCluster *c, zVec (* error_fp)(zVec,zVec,void*,zVec), void *err_util, zVec err)
{
  zVecListCell *vc;
  double dmax = 0, d;

  zListForEach( &c->vl, vc ){
    error_fp( vc->data, c->mean, err_util, err );
    if( ( d = zVecNorm( err ) ) > dmax ) dmax = d;
  }
  return dmax;
}

/* (static)
 * _zMClusterXMeansTest
 * - test if the subclusters better classify the original cluster.
 */
bool _zMClusterXMeansTest(zMCluster *mc, zCluster *org, zMCluster *submc, void *err_util)
{
  zCluster *c1, *c2;
  double fr1, fr2;
  double r, r1, r2;

  c1 = &zListHead(&submc->cl)->data;
  c2 = &zListTail(&submc->cl)->data;
  if( zListNum(&c1->vl) < 0.1 * zListNum(&c2->vl) ||
      zListNum(&c2->vl) < 0.1 * zListNum(&c1->vl) )
    return false; /* avoid too crisp cluster */

  /* original cluster */
  r = _zClusterMaxError( org, mc->met._error_fp, err_util, mc->met._err );
  fr1 = log(zListNum(&org->vl)) - mc->met._errorsize*log(r);
  /* bidivided clusters */
  r1 = _zClusterMaxError( c1, mc->met._error_fp, err_util, mc->met._err );
  r2 = _zClusterMaxError( c2, mc->met._error_fp, err_util, mc->met._err );
  fr2 = log( zListNum(&org->vl) )
      - log( pow(r1,mc->met._errorsize) + pow(r2,mc->met._errorsize) );
  return fr1 < fr2 ? true : false;
}

/* (static)
 * _zMClusterXMeans
 * - cluster vectors by X-means in a recursive way based on filling rate.
 */
int _zMClusterXMeans(zMCluster *mc, zClusterListCell *cc, void *mean_util, void *err_util)
{
  zMCluster submc;
  int iter, iter1, iter2;

  /* avoid a too small / sparse cluster */
  if( zListNum(&cc->data.vl) < Z_XMEANS_MINSIZE ) return 0;

  if( zMClusterInit( &submc, mc->met._meansize, mc->met._mean_fp, mc->met._errorsize, mc->met._error_fp ) == NULL )
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

/* zMClusterXMeans
 * - cluster vectors by X-means based on filling rate.
 */
int zMClusterXMeans(zMCluster *mc, zVecList *points, void *mean_util, void *err_util)
{
  if( zMClusterAlloc( mc, 1 ) == NULL ) return -1;
  _zMClusterXMeansInit( mc, points, mean_util, err_util );
  return _zMClusterXMeans( mc, zListHead(&mc->cl), mean_util, err_util );
}

/* (static)
 * _zClusterXMeansPDF
 * - probability density function based on Gaussian distribution.
 */
double _zClusterXMeansPDF(zCluster *c, zMCluster *mc, void *err_util, zVec x)
{
  mc->met._error_fp( x, c->mean, err_util, mc->met._err );
  return exp( -0.5*zVecSqrNorm(mc->met._err)/c->var ) / sqrt(c->var);
}

/* (static)
 * _zMClusterXMeansBICTest
 * - test if the subclusters better classify the original based on BIC.
 */
bool _zMClusterXMeansBICTest(zMCluster *mc, zCluster *c, zMCluster *submc, void *err_util)
{
  zCluster *c1, *c2;
  zVecListCell *vc;
  int q, n;
  double ls = 0, bic1, bic2;

  q = mc->met._meansize + 1;
  n = zListNum(&c->vl);
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

/* (static)
 *_zMClusterXMeansBIC
 * - cluster vectors by X-means based on BIC.
 */
int _zMClusterXMeansBIC(zMCluster *mc, zClusterListCell *cc, void *mean_util, void *err_util)
{
  zMCluster submc;
  int iter, iter1, iter2;

  /* ignore a too small cluster */
  if( zListNum(&cc->data.vl) < Z_XMEANS_MINSIZE ) return 0;
  if( zMClusterInit( &submc, mc->met._meansize, mc->met._mean_fp, mc->met._errorsize, mc->met._error_fp ) == NULL )
    return 0;
  iter = zMClusterKMeans( &submc, &cc->data.vl, 2, mean_util, err_util );
  if( _zMClusterXMeansBICTest( mc, &cc->data, &submc, err_util ) == true ){
    zMClusterDestroy( &submc );
    return 0;
  }
  iter1 = _zMClusterXMeansBIC( &submc, zListHead(&submc.cl), mean_util, err_util );
  iter2 = _zMClusterXMeansBIC( &submc, zListTail(&submc.cl), mean_util, err_util );
  _zMClusterMerge( mc, cc, &submc );
  return iter + iter1 + iter2;
}

/* zMClusterXMeansBIC
 * - cluster vectors based on BIC.
 */
int zMClusterXMeansBIC(zMCluster *mc, zVecList *points, void *mean_util, void *err_util)
{
  if( zMClusterAlloc( mc, 1 ) == NULL ) return -1;
  _zMClusterXMeansInit( mc, points, mean_util, err_util );
  return _zMClusterXMeansBIC( mc, zListHead(&mc->cl), mean_util, err_util );
}
