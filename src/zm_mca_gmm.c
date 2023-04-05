/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_mca_gmm - multiple classification analysis:
 * Gaussian mixture model
 */

#include <zm/zm_mca.h>

/* ********************************************************** */
/* Gaussian model unit class
 * ********************************************************** */

/* initialize a unit Gaussian model. */
zGMMUnit *zGMMUnitInit(zGMMUnit *gu)
{
  gu->core = NULL;
  gu->cov = NULL;
  gu->weight = HUGE_VAL;
  gu->_cov_inv = NULL;
  gu->_ci_e = NULL;
  gu->_cov_det = HUGE_VAL;
  return gu;
}

/* allocate internal vector for core of a Gaussian model. */
static zGMMUnit *_zGMMUnitCoreAlloc(zGMMUnit *gu, int size)
{
  zVecFree( gu->core );
  return ( gu->core = zVecAlloc( size ) ) ? gu : NULL;
}

/* allocate internal vectors and matrices for error of a Gaussian model. */
static zGMMUnit *_zGMMUnitErrorAlloc(zGMMUnit *gu, int size)
{
  zMatFree( gu->cov );
  zMatFree( gu->_cov_inv );
  zVecFree( gu->_ci_e );
  gu->cov = zMatAllocSqr( size );
  gu->_cov_inv = zMatAllocSqr( size );
  gu->_ci_e = zVecAlloc( size );
  return ( gu->cov && gu->_cov_inv && gu->_ci_e ) ? gu : NULL;
}

/* allocate internal vectors and matrices of a unit Gaussian model. */
zGMMUnit *zGMMUnitAlloc(zGMMUnit *gu, int coresize, int errorsize)
{
  if( !_zGMMUnitCoreAlloc( gu, coresize ) ||
      !_zGMMUnitErrorAlloc( gu, errorsize ) ) return NULL;
  gu->weight = 0;
  return gu;
}

/* free internal vectors and matrices of a unit Gaussian model. */
void zGMMUnitFree(zGMMUnit *gu)
{
  zVecFree( gu->core );
  zMatFree( gu->cov );
  zMatFree( gu->_cov_inv );
  zVecFree( gu->_ci_e );
}

/* estimate a unit Gaussian model from a cluster of points. */
static void _zGMMUnitEstim(zGMMUnit *gu, zVecList *points, zClusterMethod *method)
{
  zVecListCell *pc;

  method->core_fp( method, points, method->core_util, gu->core );
  zMatZero( gu->cov );
  zListForEach( points, pc ){
    method->error_fp( method, pc->data, gu->core, method->error_util, method->error );
    zMatAddDyadNC( gu->cov, method->error, method->error );
  }
  zMatDivNCDRC( gu->cov, zListSize(points) );
  if( zIsTiny( gu->_cov_det = zMatDet( gu->cov ) ) )
    zMatZero( gu->_cov_inv );
  else
    zMatInv( gu->cov, gu->_cov_inv );
}

/* Gaussian probability density function. */
static double _zGMMUnitPDF(zGMMUnit *gu, zVec p, zClusterMethod *method)
{
  method->error_fp( method, p, gu->core, method->error_util, method->error );
  zMulMatVec( gu->_cov_inv, method->error, gu->_ci_e );
  return zIsTiny(gu->_cov_det) ? 1.0 :
    exp( -0.5*zVecInnerProd( method->error, gu->_ci_e ) )
      / sqrt( pow(zPIx2,zVecSizeNC(method->error)) * gu->_cov_det );
}

/* a loaded covariance matrix of a Gaussian mixture model. */
static zMat _zGMMUnitLoadedCov(zGMMUnit *gu, zVecList *points, double load[], double nk, zClusterMethod *method)
{
  zVecListCell *pc;
  int i = 0;

  zMatZero( gu->cov );
  zListForEach( points, pc ){
    method->error_fp( method, pc->data, gu->core, method->error_util, method->error );
    zVecMulNC( method->error, load[i], method->_lerr );
    zMatAddDyadNC( gu->cov, method->_lerr, method->error );
    i++;
  }
  zMatDivDRC( gu->cov, nk );
  zMatInv( gu->cov, gu->_cov_inv );
  gu->_cov_det = zMatDet( gu->cov );
  return gu->cov;
}

/* ********************************************************** */
/* Gaussian mixture model class
 * ********************************************************** */

static zGMM *_zGMMCoreAlloc(zGMM *gmm, int size)
{
  zGMMListCell *gc;

  zListForEach( &gmm->glist, gc )
    if( !_zGMMUnitCoreAlloc( &gc->data, size ) ) return NULL;
  return gmm;
}

static zGMM *_zGMMErrorAlloc(zGMM *gmm, int size)
{
  zGMMListCell *gc;

  zListForEach( &gmm->glist, gc )
    if( !_zGMMUnitErrorAlloc( &gc->data, size ) ) return NULL;
  return gmm;
}

/* initialize a Gaussian mixture model. */
zGMM *zGMMInit(zGMM *gmm, int k, int size)
{
  int i;
  zGMMListCell *gc;

  zListInit( &gmm->glist );
  for( i=0; i<k; i++ ){
    if( ( gc = zAlloc( zGMMListCell, 1 ) ) == NULL ){
      ZALLOCERROR();
      goto ERR;
    }
    zGMMUnitInit( &gc->data );
    if( zGMMUnitAlloc( &gc->data, size, size ) == NULL ){
      ZALLOCERROR();
      free( gc );
      goto ERR;
    }
    zListInsertTail( &gmm->glist, gc );
  }
  gmm->log_likelihood = -HUGE_VAL;
  if( zClusterMethodCreate( &gmm->method, size ) == NULL )
    goto ERR;
  return gmm;

 ERR:
  zClusterMethodInit( &gmm->method );
  zGMMDestroy( gmm );
  return NULL;
}

/* set error function for Gaussian mixture model. */
bool zGMMSetErrorFunc(zGMM *gmm, int size, zVec (* error_fp)(zClusterMethod*,zVec,zVec,void*,zVec), void *util)
{
  if( !zClusterMethodSetErrorFunc( &(gmm)->method, size, error_fp, util ) ) return false;
  return _zGMMErrorAlloc( gmm, size ) ? true : false;
}

/* set distance function for Gaussian mixture model. */
bool zGMMSetDistFunc(zGMM *gmm, double (* dist_fp)(zClusterMethod*,zVec,zVec,void*), void *util)
{
  return zClusterMethodSetDistFunc( &(gmm)->method, dist_fp, util ) ? true : false;
}

/* set core function for Gaussian mixture model. */
bool zGMMSetCoreFunc(zGMM *gmm, int size, zVec (* core_fp)(zClusterMethod*,zVecAddrList*,void*,zVec), void *util)
{
  if( !zClusterMethodSetCoreFunc( &(gmm)->method, size, core_fp, util ) ) return false;
  return _zGMMCoreAlloc( gmm, size ) ? true : false;
}

/* set loaded mean function for Gaussian mixture model. */
bool zGMMSetLoadedMeanFunc(zGMM *gmm, zVec (* lm_fp)(zClusterMethod*,zVecAddrList*,double[],double,void*,zVec), void *util)
{
  return zClusterMethodSetLoadedMeanFunc( &(gmm)->method, lm_fp, util ) ? true : false;
}

/* destroy a Gaussian mixture model. */
void zGMMDestroy(zGMM *gmm)
{
  zGMMListCell *gc;

  while( !zListIsEmpty( &gmm->glist ) ){
    zListDeleteTail( &gmm->glist, &gc );
    zGMMUnitFree( &gc->data );
    free( gc );
  }
  zClusterMethodDestroy( &gmm->method );
}

#ifdef DEBUG
/* print a Gaussian mixture model out to a file. */
static void _zGMMFPrint(FILE *fp, zGMM *gmm)
{
  zGMMListCell *gc;

  zListForEach( &gmm->glist, gc ){
    zVecFPrint( fp, gc->data.core );
    zMatFPrint( fp, gc->data.cov );
    printf( "%g\n", gc->data.weight );
  }
}
#endif /* DEBUG */

/* expectation phase of EM algorithm to estimate a Gaussian mixture model. */
static bool _zGMMCreateEMExpect(zGMM *gmm, zVecList *points, zMat pdf, zMat load, zVec load_det, zVec nk)
{
  int i, j;
  zGMMListCell *gc;
  zVecListCell *pc;
  double p, gamma, eps = 0;

  zVecZero( load_det );
  i = 0;
  zListForEach( points, pc ){
    j = 0;
    zListForEach( &gmm->glist, gc ){
      p = _zGMMUnitPDF( &gc->data, pc->data, &gmm->method );
      zMatElemNC(pdf,j,i) = p;
      zVecElemNC(load_det,i) += gc->data.weight * p;
      j++;
    }
    i++;
  }
  zVecZero( nk );
  i = 0;
  zListForEach( points, pc ){
    j = 0;
    zListForEach( &gmm->glist, gc ){
      gamma = gc->data.weight * zMatElemNC(pdf,j,i) / zVecElemNC(load_det,i);
      eps += gamma - zMatElemNC(load,j,i);
      zMatElemNC(load,j,i) = gamma;
      zVecElemNC(nk,j) += gamma;
      j++;
    }
    i++;
  }
  return zIsTiny( eps ) ? true : false;
}

/* log-likelihood of a Gaussian mixture model. */
static bool _zGMMLogLikelihood(zGMM *gmm, zVec load_det)
{
  double l = 0;
  int i;

  for( i=0; i<zVecSizeNC(load_det); i++ )
    l += log( zVecElemNC(load_det,i) );
  if( zIsTiny( gmm->log_likelihood - l ) ) return true;
  gmm->log_likelihood = l;
  return false;
}

/* maximization phase of EM algorithm to estimate a Gaussian mixture model. */
static void _zGMMCreateEMMaximize(zGMM *gmm, zVecList *points, zMat load, zVec nk)
{
  zGMMListCell *gc;
  int j;

  j = 0;
  zListForEach( &gmm->glist, gc ){
    zGMMLoadedMeanF( gmm, points, zMatRowBuf(load,j), zVecElemNC(nk,j), gc->data.core );
    _zGMMUnitLoadedCov( &gc->data, points, zMatRowBuf(load,j), zVecElemNC(nk,j), &gmm->method );
    gc->data.weight = zVecElemNC(nk,j) / zListSize(points);
    j++;
  }
}

/* create a Gaussian mixture model based on EM algorithm. */
zGMM *zGMMCreateEM(zGMM *gmm, zVecList *points)
{
  zMCluster mc;
  zClusterListCell *cc;
  zGMMListCell *gc;
  zMat pdf, load;
  zVec load_det, nk;
  int k, i, iter = 0;

  if( !zMClusterInit( &mc, zVecSize(gmm->method.error) ) ||
      !zClusterMethodCopy( &gmm->method, &mc.method ) )
    return NULL;

  k = zListSize( &gmm->glist );
  load = zMatAlloc( k, zListSize(points) );
  pdf  = zMatAlloc( k, zListSize(points) );
  load_det = zVecAlloc( zListSize(points) );
  nk = zVecAlloc( k );
  if( !load || !pdf || !load_det || !nk ){
    ZALLOCERROR();
    gmm = NULL;
    goto TERMINATE;
  }
  /* initialize Gaussians by K-means++ */
  zMClusterKMeans( &mc, points, k );
  gc = zListTail( &gmm->glist );
  zListForEach( zMClusterClusterList(&mc), cc ){
    _zGMMUnitEstim( &gc->data, zClusterSampleList(&cc->data), &gmm->method );
    gc->data.weight = 1.0 / zListSize(zMClusterClusterList(&mc));
    gc = zListCellNext(gc);
  }
  zMClusterDestroy( &mc );
  /* iteration */
  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    _zGMMCreateEMExpect( gmm, points, pdf, load, load_det, nk );
    if( _zGMMLogLikelihood( gmm, load_det ) )
      goto TERMINATE;
    _zGMMCreateEMMaximize( gmm, points, load, nk );
  }
  ZITERWARN( iter );
 TERMINATE:
  zMatFree( load );
  zMatFree( pdf );
  zVecFree( load_det );
  zVecFree( nk );
  return gmm;
}

/* Akaike's Information Criterion. */
double zGMMAIC(zGMM *gmm)
{
  return 2 * ( gmm->method.core_size * zListSize(&gmm->glist) - gmm->log_likelihood );
}

/* Bayesian Information Criterion. */
double zGMMBIC(zGMM *gmm, zVecAddrList *sample)
{
  return gmm->method.core_size * zListSize(&gmm->glist) * log( zListSize(sample) ) -2 * gmm->log_likelihood;
}
