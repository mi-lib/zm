/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_mva_ransac - multivariate analysis analysis: random sample consensus.
 */

#include <zm/zm_mva.h>

/* randomly select test data from original data set. */
static bool _zRANSACSelectRandom(zVecList *sample, zVecList *selected, int n)
{
  zVecListCell *sp;
  int i;

  if( n > zListSize(sample) ){
    ZRUNERROR( ZM_ERR_INVALID_NUMSAMP, n, zListSize(sample) );
    return false;
  }
  zListInit( selected );
  while( zListSize(selected) < n ){
    i = zRandI( 0, zListSize(sample)-1 );
    zListItem( sample, i, &sp );
    zListPurge( sample, sp );
    zListInsertHead( selected, sp );
  }
  return true;
}

/* count number of inliers with respect to a guess. */
static int _zRANSACCountInlier(zVec q, const zVecList *sample, double (* error_fp)(const zVec,const zVec,void*), void *util, double th)
{
  zVecListCell *sp;
  int count = 0;

  zListForEach( sample, sp ){
    if( fabs( error_fp( q, sp->data, util ) ) < th ) count++;
  }
  return count;
}

/* resample test data for model refinement from original data set. */
static int _zRANSACSelectInlier(zVecList *inlier, zVec q, zVecList *sample, double (* error_fp)(const zVec,const zVec,void*), void *util, double th)
{
  zVecListCell *sp, *sp_prev;

  zListInit( inlier );
  zListForEach( sample, sp ){
    if( fabs( error_fp( q, sp->data, util ) ) < th ){
      sp_prev = zListCellPrev( sp );
      zListPurge( sample, sp );
      zListInsertHead( inlier, sp );
      sp = sp_prev;
    }
  }
  return zListSize(inlier);
}

/* apply RANSAC and save inliers out of the original samples. */
zVec zRANSACSaveInlier(zVec q, zVecList *sample, zVec (* fit_fp)(zVec,const zVecList*,void*), double (* error_fp)(const zVec,const zVec,void*), void *util, int ns, int nt, double th, zVecList *inlier)
{
  zVec qt;
  int count, count_prev = 0;

  if( !( qt = zVecAlloc( zVecSizeNC(q) ) ) ) return NULL;
  if( ns < zVecSizeNC(q) )
    ZRUNWARN( ZM_WARN_INSUFFICIENT_SAMPLES, ns, zVecSizeNC(q) );
  /* find a candidate of the most likely model */
  while( --nt >= 0 ){
    if( !_zRANSACSelectRandom( sample, inlier, ns ) ) break;
    fit_fp( qt, inlier, util );
    zListAppend( sample, inlier );
    if( ( count = _zRANSACCountInlier( qt, sample, error_fp, util, th ) ) > count_prev ){
      count_prev = count;
      zVecCopyNC( qt, q );
    }
  }
  zVecFree( qt );
  /* refine the model with the doubled threshold */
  _zRANSACSelectInlier( inlier, q, sample, error_fp, util, th*2 );
  fit_fp( q, inlier, util );
  return q;
}

/* apply RANSAC. */
zVec zRANSAC(zVec q, zVecList *sample, zVec (* fit_fp)(zVec,const zVecList*,void*), double (* error_fp)(const zVec,const zVec,void*), void *util, int ns, int nt, double th)
{
  zVecList inlier;

  q = zRANSACSaveInlier( q, sample, fit_fp, error_fp, util, ns, nt, th, &inlier );
  zListAppend( sample, &inlier );
  return q;
}

#define Z_RANSAC_SAMPLE_MIN_RATE 5
#define Z_RANSAC_SAMPLE_MAX_RATE 0.5
#define Z_RANSAC_TRIAL_RATE      3
#define Z_RANSAC_THRESHOLD_RATE  0.1

/* compute selection number from the size of parameters and the number of samples. */
static int _zRANSACAutoSelectionNumber(const zVec q, const zVecList *sample, double rate)
{
  return zMin( zVecSizeNC(q)*Z_RANSAC_SAMPLE_MIN_RATE, zListSize(sample)*(1-rate)*Z_RANSAC_SAMPLE_MAX_RATE );
}

/* RANSAC with an automatic setting of parameters */
zVec zRANSACAuto(zVec q, zVecList *sample, zVec (* fit_fp)(zVec,const zVecList*,void*), double (* error_fp)(const zVec,const zVec,void*), void *util, double rate, double nl)
{
  int ns;

  ns = _zRANSACAutoSelectionNumber( q, sample, rate );
  return zRANSAC( q, sample, fit_fp, error_fp, util,
    ns, ns*Z_RANSAC_TRIAL_RATE, Z_RANSAC_THRESHOLD_RATE*nl );
}

/* RANSAC with an automatic setting of parameters */
zVec zRANSACSaveInlierAuto(zVec q, zVecList *sample, zVec (* fit_fp)(zVec,const zVecList*,void*), double (* error_fp)(const zVec,const zVec,void*), void *util, double rate, double nl, zVecList *inlier)
{
  int ns;

  ns = _zRANSACAutoSelectionNumber( q, sample, rate );
  return zRANSACSaveInlier( q, sample, fit_fp, error_fp, util,
    ns, ns*Z_RANSAC_TRIAL_RATE, Z_RANSAC_THRESHOLD_RATE*nl, inlier );
}
