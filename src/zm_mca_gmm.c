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

static void _zGMMUnitEstim(zGMMUnit *gu, zVecList *points, zClusterMethod *met, void *mean_util, void *err_util);
static double _zGMMUnitPDF(zGMMUnit *gu, zVec p, zClusterMethod *met, void *err_util);
static zMat _zGMMUnitLoadedCov(zGMMUnit *gu, zVecList *points, double load[], double nk, zClusterMethod *met, void *err_util);

/* zGMMUnitAlloc
 * - allocate internal vectors and matrices of a unit Gaussian function.
 */
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

/* zGMMUnitFree
 * - free internal vectors and matrices of a unit Gaussian function.
 */
void zGMMUnitFree(zGMMUnit *gu)
{
  zVecFree( gu->mean );
  zMatFree( gu->cov );
  zMatFree( gu->_cov_inv );
  zVecFree( gu->_ci_e );
}

/* (static)
 * _zGMMUnitEstim
 * - estimate a unit Gaussian function from a cluster of points.
 */
void _zGMMUnitEstim(zGMMUnit *gu, zVecList *points, zClusterMethod *met, void *mean_util, void *err_util)
{
  zVecListCell *pc;

  met->_mean_fp( points, mean_util, gu->mean );
  zMatClear( gu->cov );
  zListForEach( points, pc ){
    met->_error_fp( pc->data, gu->mean, err_util, met->_err );
    zMatAddDyadNC( gu->cov, met->_err, met->_err );
  }
  zMatDivNCDRC( gu->cov, zListNum(points) );
  if( zIsTiny( gu->_cov_det = zMatDet( gu->cov ) ) )
    zMatClear( gu->_cov_inv );
  else
    zMatInv( gu->cov, gu->_cov_inv );
}

/* (static)
 * _zGMMUnitPDF
 * - Gaussian probability density function.
 */
double _zGMMUnitPDF(zGMMUnit *gu, zVec p, zClusterMethod *met, void *err_util)
{
  met->_error_fp( p, gu->mean, err_util, met->_err );
  zMulMatVec( gu->_cov_inv, met->_err, gu->_ci_e );
  return zIsTiny(gu->_cov_det) ? 1.0 :
    exp( -0.5*zVecInnerProd( met->_err, gu->_ci_e ) )
      / sqrt( pow(zPIx2,met->_errorsize) * gu->_cov_det );
}

/* (static)
 * _zGMMUnitLoadedCov
 * - a loaded covariance matrix of a Gaussian mixture model
 */
zMat _zGMMUnitLoadedCov(zGMMUnit *gu, zVecList *points, double load[], double nk, zClusterMethod *met, void *err_util)
{
  zVecListCell *pc;
  register int i = 0;

  zMatClear( gu->cov );
  zListForEach( points, pc ){
    met->_error_fp( pc->data, gu->mean, err_util, met->_err );
    zVecMulNC( met->_err, load[i], met->_err_l );
    zMatAddDyadNC( gu->cov, met->_err_l, met->_err );
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

static bool _zGMMCreateEMExpect(zGMM *gmm, zVecList *points, zMat pdf, zMat load, zVec load_det, zVec nk, void *err_util);
static bool _zGMMLogLikelihood(zGMM *gmm, zVec load_det);
static void _zGMMCreateEMMaximize(zGMM *gmm, zVecList *points, zMat load, zVec nk, void *mean_util, void *err_util);

/* zGMMInit
 * - initialize a Gaussian mixture model.
 */
zGMM *zGMMInit(zGMM *gmm, int k, int meansize, zVec (* mean_fp)(zVecList*,void*,zVec), zVec (* mean_l_fp)(zVecList*,double[],double,void*,zVec), int errorsize, zVec (* error_fp)(zVec,zVec,void*,zVec))
{
  register int i;
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
  if( zClusterMethodCreate( &gmm->met, meansize, mean_fp, errorsize, error_fp ) == NULL )
    goto ERR;
  if( zClusterMethodLoadedCreate( &gmm->met, mean_l_fp ) == NULL )
    goto ERR;
  return gmm;

 ERR:
  zClusterMethodInit( &gmm->met );
  zGMMDestroy( gmm );
  return NULL;
}

/* zGMMDestroy
 * - destroy a Gaussian mixture model.
 */
void zGMMDestroy(zGMM *gmm)
{
  zGMMListCell *gc;

  while( !zListIsEmpty( &gmm->gl ) ){
    zListDeleteTail( &gmm->gl, &gc );
    zGMMUnitFree( &gc->data );
    free( gc );
  }
  zClusterMethodDestroy( &gmm->met );
}

#ifdef DEBUG
/* (static)
 * _zGMMFWrite
 * - output a Gaussian mixture model to a file.
 */
static void _zGMMFWrite(FILE *fp, zGMM *gmm);
void _zGMMFWrite(FILE *fp, zGMM *gmm)
{
  zGMMListCell *gc;

  zListForEach( &gmm->gl, gc ){
    zVecFWrite( fp, gc->data.mean );
    zMatFWrite( fp, gc->data.cov );
    printf( "%g\n", gc->data.weight );
  }
}
#endif /* DEBUG */

/* (static)
 * _zGMMCreateEMExpect
 * - expectation phase of EM algorithm to estimate a Gaussian mixture model.
 */
bool _zGMMCreateEMExpect(zGMM *gmm, zVecList *points, zMat pdf, zMat load, zVec load_det, zVec nk, void *err_util)
{
  register int i, j;
  zGMMListCell *gc;
  zVecListCell *pc;
  double p, gamma, eps = 0;

  zVecClear( load_det );
  i = 0;
  zListForEach( points, pc ){
    j = 0;
    zListForEach( &gmm->gl, gc ){
      if( !gc->data.active ){
        zMatElem(pdf,j,i) = 0;
        zVecElem(load_det,i) = HUGE_VAL;
      } else{
        p = _zGMMUnitPDF( &gc->data, pc->data, &gmm->met, err_util );
        zMatElem(pdf,j,i) = p;
        zVecElem(load_det,i) += gc->data.weight * p;
      }
      j++;
    }
    i++;
  }
  zVecClear( nk );
  i = 0;
  zListForEach( points, pc ){
    j = 0;
    zListForEach( &gmm->gl, gc ){
      if( gc->data.active ){
        gamma = gc->data.weight * zMatElem(pdf,j,i) / zVecElem(load_det,i);
        eps += gamma - zMatElem(load,j,i);
        zMatElem(load,j,i) = gamma;
        zVecElem(nk,j) += gamma;
      }
      j++;
    }
    i++;
  }
  return zIsTiny( eps ) ? true : false;
}

/* (static)
 * _zGMMLogLikelihood
 * - log-likelihood of a Gaussian mixture model.
 */
bool _zGMMLogLikelihood(zGMM *gmm, zVec load_det)
{
  double l = 0;
  register int i;

  for( i=0; i<zVecSizeNC(load_det); i++ )
    l += log( zVecElem(load_det,i) );
  if( zIsTiny( gmm->log_likelihood - l ) ) return true;
  gmm->log_likelihood = l;
  return false;
}

/* (static)
 * _zGMMCreateEMMaximize
 * - maximization phase of EM algorithm to estimate a Gaussian mixture model.
 */
void _zGMMCreateEMMaximize(zGMM *gmm, zVecList *points, zMat load, zVec nk, void *mean_util, void *err_util)
{
  zGMMListCell *gc;
  register int j;

  j = 0;
  zListForEach( &gmm->gl, gc ){
    if( gc->data.active ){
      gmm->met._mean_l_fp( points, zMatRowBuf(load,j), zVecElem(nk,j), mean_util, gc->data.mean );
      _zGMMUnitLoadedCov( &gc->data, points, zMatRowBuf(load,j), zVecElem(nk,j), &gmm->met, err_util );
      gc->data.weight = zVecElem(nk,j) / zListNum(points);
    }
    j++;
  }
}

/* zGMMCreateEM
 * - create a Gaussian mixture model based on EM algorithm.
 */
zGMM *zGMMCreateEM(zGMM *gmm, zVecList *points, int k, void *mean_util, void *err_util)
{
  zMCluster mc;
  zClusterListCell *cc;
  zGMMListCell *gc;
  zMat pdf, load;
  zVec load_det, nk;
  int iter = 0;
  register int i;

  if( zMClusterInit( &mc, gmm->met._meansize, gmm->met._mean_fp, gmm->met._errorsize, gmm->met._error_fp ) == NULL )
    return NULL;
  load = zMatAlloc( k, zListNum(points) );
  pdf  = zMatAlloc( k, zListNum(points) );
  load_det = zVecAlloc( zListNum(points) );
  nk = zVecAlloc( k );
  if( !load || !pdf || !load_det || !nk ){
    ZALLOCERROR();
    gmm = NULL;
    goto TERMINATE;
  }
  /* initialize Gaussians by K-means */
  zMClusterKMeans( &mc, points, k, mean_util, err_util );
  gc = zListTail( &gmm->gl );
  zListForEach( &mc.cl, cc ){
    if( zListIsEmpty(&cc->data.vl) ){
      gc->data.active = false;
    } else{
      gc->data.active = true;
      _zGMMUnitEstim( &gc->data, &cc->data.vl, &gmm->met, mean_util, err_util );
      gc->data.weight = 1.0 / zListNum(&mc.cl);
    }
    gc = zListCellNext(gc);
  }
  zMClusterDestroy( &mc );
  /* iteration */
  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    _zGMMCreateEMExpect( gmm, points, pdf, load, load_det, nk, err_util );
    if( _zGMMLogLikelihood( gmm, load_det ) )
      goto TERMINATE;
    _zGMMCreateEMMaximize( gmm, points, load, nk, mean_util, err_util );
  }
  ZITERWARN( iter );
 TERMINATE:
  zMatFree( load );
  zMatFree( pdf );
  zVecFree( load_det );
  zVecFree( nk );
  return gmm;
}
