/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_data_ransac - data analysis: random sampling consensus.
 */

#include <zm/zm_data.h>

/* randomly resample test data from original data set. */
static bool _zRANSACResampleRandom(zVecList *sample, zVecList *va, int n)
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
static int _zRANSACResampleRefine(zVecList *va, zVec q, zVecList *sample, double (* error_fp)(zVec,zVec,void*), void *util, double th)
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

/* RANSAC: random sampling consensus. */
zVec zRANSAC(zVec q, zVecList *sample, zVec (* fit_fp)(zVec,zVecList*,void*), double (* error_fp)(zVec,zVec,void*), void *util, int ns, int nt, double th)
{
  zVec qt;
  zVecList va;
  int count, count_prev = 0;

  if( !( qt = zVecAlloc( zVecSizeNC(q) ) ) ) return NULL;
  if( ns < zVecSizeNC(q) )
    ZRUNWARN( ZM_WARN_INSUFFICIENT_SAMPLES, ns, zVecSizeNC(q) );
  /* find a candidate of the most likely model. */
  while( --nt >= 0 ){
    if( !_zRANSACResampleRandom( sample, &va, ns ) ) break;
    fit_fp( qt, &va, util );
    zListAppend( sample, &va );
    if( ( count = _zRANSACCountInlier( qt, sample, error_fp, util, th ) ) > count_prev ){
      count_prev = count;
      zVecCopyNC( qt, q );
    }
  }
  zVecFree( qt );
  /* refine the model. */
  _zRANSACResampleRefine( &va, q, sample, error_fp, util, th );
  fit_fp( q, &va, util );
  zListAppend( sample, &va );
  return q;
}
