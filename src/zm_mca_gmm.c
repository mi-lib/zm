/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_mca_gmm - multiple classification analysis:
 * Gaussian mixture model
 */

#include <zm/zm_mca.h>

/* ********************************************************** */
/* unit Gaussian function class
 * ********************************************************** */

/* allocate internal vectors and matrices of a unit Gaussian function. */
zGMMUnit *zGMMUnitAlloc(zGMMUnit *gu, int meansize, int errorsize)
{
  gu->mean = zVecAlloc( meansize );
  gu->cov = zMatAllocSqr( errorsize );
  gu->_cov_inv = zMatAllocSqr( errorsize );
  gu->_ci_e = zVecAlloc( errorsize );
  if( gu->mean == NULL || gu->cov == NULL || gu->_cov_inv == NULL || gu->_ci_e == NULL ){
    zGMMUnitFree( gu );
    return NULL;
  }
  gu->weight = 0;
  gu->active = false;
  return gu;
}

/* free internal vectors and matrices of a unit Gaussian function. */
void zGMMUnitFree(zGMMUnit *gu)
{
  zVecFree( gu->mean );
  zMatFree( gu->cov );
  zMatFree( gu->_cov_inv );
  zVecFree( gu->_ci_e );
}

/* estimate a unit Gaussian function from a cluster of points. */
static void _zGMMUnitEstim(zGMMUnit *gu, zVecList *points, zClusterMethod *method)
{
  zVecListCell *pc;

  method->core_fp( method, points, method->core_util, gu->mean );
  zMatZero( gu->cov );
  zListForEach( points, pc ){
    method->error_fp( method, pc->data, gu->mean, method->error_util, method->error );
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
  method->error_fp( method, p, gu->mean, method->error_util, method->error );
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
    method->error_fp( method, pc->data, gu->mean, method->error_util, method->error );
    zVecMulNC( method->error, load[i], method->lerror );
    zMatAddDyadNC( gu->cov, method->lerror, method->error );
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

/* initialize a Gaussian mixture model. */
zGMM *zGMMInit(zGMM *gmm, int k, int meansize, zVec (* mean_fp)(zClusterMethod*,zVecList*,void*,zVec), void *mean_util, zVec (* l_mean_fp)(zClusterMethod*,zVecList*,double[],double,void*,zVec), void *l_mean_util, int errorsize, zVec (* error_fp)(zClusterMethod*,zVec,zVec,void*,zVec), void *error_util, double (* dist_fp)(zClusterMethod*,zVec,zVec,void*), void *dist_util)
{
  int i;
  zGMMListCell *gc;

  zListInit( &gmm->gl );
  for( i=0; i<k; i++ ){
    if( ( gc = zAlloc( zGMMListCell, 1 ) ) == NULL ){
      ZALLOCERROR();
      goto ERR;
    }
    if( zGMMUnitAlloc( &gc->data, meansize, errorsize ) == NULL ){
      ZALLOCERROR();
      free( gc );
      goto ERR;
    }
    zListInsertTail( &gmm->gl, gc );
  }
  gmm->log_likelihood = -HUGE_VAL;
  if( zClusterMethodCreate( &gmm->method, meansize, errorsize ) == NULL )
    goto ERR;
  if( !zClusterMethodSetCoreFunc( &gmm->method, meansize, mean_fp, mean_util ) ||
      !zClusterMethodSetErrorFunc( &gmm->method, errorsize, error_fp, error_util ) ||
      !zClusterMethodSetDistFunc( &gmm->method, dist_fp, dist_util ) ||
      !zClusterMethodSetLoadedMeanFunc( &gmm->method, l_mean_fp, l_mean_util ) )
    goto ERR;
  return gmm;

 ERR:
  zClusterMethodInit( &gmm->method );
  zGMMDestroy( gmm );
  return NULL;
}

/* destroy a Gaussian mixture model. */
void zGMMDestroy(zGMM *gmm)
{
  zGMMListCell *gc;

  while( !zListIsEmpty( &gmm->gl ) ){
    zListDeleteTail( &gmm->gl, &gc );
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

  zListForEach( &gmm->gl, gc ){
    zVecFPrint( fp, gc->data.mean );
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
    zListForEach( &gmm->gl, gc ){
      if( !gc->data.active ){
        zMatElemNC(pdf,j,i) = 0;
        zVecElemNC(load_det,i) = HUGE_VAL;
      } else{
        p = _zGMMUnitPDF( &gc->data, pc->data, &gmm->method );
        zMatElemNC(pdf,j,i) = p;
        zVecElemNC(load_det,i) += gc->data.weight * p;
      }
      j++;
    }
    i++;
  }
  zVecZero( nk );
  i = 0;
  zListForEach( points, pc ){
    j = 0;
    zListForEach( &gmm->gl, gc ){
      if( gc->data.active ){
        gamma = gc->data.weight * zMatElemNC(pdf,j,i) / zVecElemNC(load_det,i);
        eps += gamma - zMatElemNC(load,j,i);
        zMatElemNC(load,j,i) = gamma;
        zVecElemNC(nk,j) += gamma;
      }
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
  zListForEach( &gmm->gl, gc ){
    if( gc->data.active ){
      gmm->method.lm_fp( &gmm->method, points, zMatRowBuf(load,j), zVecElemNC(nk,j),gmm->method.core_util, gc->data.mean );
      _zGMMUnitLoadedCov( &gc->data, points, zMatRowBuf(load,j), zVecElemNC(nk,j), &gmm->method );
      gc->data.weight = zVecElemNC(nk,j) / zListSize(points);
    }
    j++;
  }
}

/* create a Gaussian mixture model based on EM algorithm. */
zGMM *zGMMCreateEM(zGMM *gmm, zVecList *points, int k)
{
  zMCluster mc;
  zClusterListCell *cc;
  zGMMListCell *gc;
  zMat pdf, load;
  zVec load_det, nk;
  int i, iter = 0;

  if( !zMClusterInit( &mc, gmm->method.core_size, zVecSize(gmm->method.error) ) ||
      !zClusterMethodCopy( &gmm->method, &mc.method ) )
    return NULL;

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
  gc = zListTail( &gmm->gl );
  zListForEach( zMClusterClusterList(&mc), cc ){
    if( zListIsEmpty(zClusterSampleList(&cc->data)) ){
      gc->data.active = false;
    } else{
      gc->data.active = true;
      _zGMMUnitEstim( &gc->data, zClusterSampleList(&cc->data), &gmm->method );
      gc->data.weight = 1.0 / zListSize(zMClusterClusterList(&mc));
    }
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
