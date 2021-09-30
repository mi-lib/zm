/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_data_ransac - data analysis: random sample consensus.
 */

#include <zm/zm_data.h>

/* randomly select test data from original data set. */
static bool _zRANSACSelectRandom(zVecList *sample, zVecList *va, int n)
{
  zVecListCell *sp;
  int i;

  if( n > zListSize(sample) ){
    ZRUNERROR( ZM_ERR_INVALID_NUMSAMP, n, zListSize(sample) );
    return false;
  }
  zListInit( va );
  while( zListSize(va) < n ){
    i = zRandI( 0, zListSize(sample)-1 );
    zListItem( sample, i, &sp );
    zListPurge( sample, sp );
    zListInsertHead( va, sp );
  }
  return true;
}

/* count number of inliers with respect to a guess. */
static int _zRANSACCountInlier(zVec q, zVecList *sample, double (* error_fp)(zVec,zVec,void*), void *util, double th)
{
  zVecListCell *sp;
  int count = 0;

  zListForEach( sample, sp ){
    if( fabs( error_fp( q, sp->data, util ) ) < th ) count++;
  }
  return count;
}

/* resample test data for model refinement from original data set. */
static int _zRANSACSelectInlier(zVecList *va, zVec q, zVecList *sample, double (* error_fp)(zVec,zVec,void*), void *util, double th)
{
  zVecListCell *sp, *sp_prev;

  zListInit( va );
  zListForEach( sample, sp ){
    if( fabs( error_fp( q, sp->data, util ) ) < th ){
      sp_prev = zListCellPrev( sp );
      zListPurge( sample, sp );
      zListInsertHead( va, sp );
      sp = sp_prev;
    }
  }
  return zListSize(va);
}

/* RANSAC: random sample consensus. */
zVec zRANSAC(zVec q, zVecList *sample, zVec (* fit_fp)(zVec,zVecList*,void*), double (* error_fp)(zVec,zVec,void*), void *util, int ns, int nt, double th)
{
  zVec qt;
  zVecList va;
  int count, count_prev = 0;

  if( !( qt = zVecAlloc( zVecSizeNC(q) ) ) ) return NULL;
  if( ns < zVecSizeNC(q) )
    ZRUNWARN( ZM_WARN_INSUFFICIENT_SAMPLES, ns, zVecSizeNC(q) );
  /* find a candidate of the most likely model */
  while( --nt >= 0 ){
    if( !_zRANSACSelectRandom( sample, &va, ns ) ) break;
    fit_fp( qt, &va, util );
    zListAppend( sample, &va );
    if( ( count = _zRANSACCountInlier( qt, sample, error_fp, util, th ) ) > count_prev ){
      count_prev = count;
      zVecCopyNC( qt, q );
    }
  }
  zVecFree( qt );
  /* refine the model with the doubled threshold */
  _zRANSACSelectInlier( &va, q, sample, error_fp, util, th*2 );
  fit_fp( q, &va, util );
  zListAppend( sample, &va );
  return q;
}

#define Z_RANSAC_SAMPLE_MIN_RATE 5
#define Z_RANSAC_SAMPLE_MAX_RATE 0.5
#define Z_RANSAC_TRIAL_RATE      3
#define Z_RANSAC_THRESHOLD_RATE  0.1

/* RANSAC with an automatic setting of parameters */
zVec zRANSACAuto(zVec q, zVecList *sample, zVec (* fit_fp)(zVec,zVecList*,void*), double (* error_fp)(zVec,zVec,void*), void *util, double r, double nl)
{
  int ns;

  ns = zMin( zVecSizeNC(q)*Z_RANSAC_SAMPLE_MIN_RATE,
             zListSize(sample)*(1-r)*Z_RANSAC_SAMPLE_MAX_RATE );
  return zRANSAC( q, sample, fit_fp, error_fp, util,
    ns, ns*Z_RANSAC_TRIAL_RATE, Z_RANSAC_THRESHOLD_RATE*nl );
}
