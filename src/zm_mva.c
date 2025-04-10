/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_mva - multivariate analysis analysis.
 */

#include <zm/zm_mva.h>
#include <zm/zm_rand.h>

/* minimum and maximum of all vectors in a list. */
bool zVecListMinMax(const zVecList *list, zVec min, zVec max)
{
  zVecListCell *vc;
  int i;

  if( !zVecSizeEqual( zListHead(list)->data, min ) ||
      !zVecSizeEqual( zListHead(list)->data, max ) ){
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );
    return false;
  }
  zVecSetAll( min, HUGE_VAL );
  zVecSetAll( max,-HUGE_VAL );
  zListForEach( list, vc ){
    for( i=0; i<zVecSizeNC(min); i++ ){
      if( zVecElemNC(min,i) > zVecElemNC(vc->data,i) )
        zVecSetElemNC( min, i, zVecElemNC(vc->data,i) );
      if( zVecElemNC(max,i) < zVecElemNC(vc->data,i) )
        zVecSetElemNC( max, i, zVecElemNC(vc->data,i) );
    }
  }
  return true;
}

/* sum up all vectors in a list. */
zVec zVecListSum(const zVecList *list, zVec sum)
{
  zVecListCell *vc;

  zVecZero( sum );
  zListForEach( list, vc )
    zVecAddDRC( sum, vc->data );
  return sum;
}

/* mean of all vectors in a list. */
zVec zVecListMean(const zVecList *list, zVec mean)
{
  zVecListSum( list, mean );
  return zVecDivDRC( mean, zListSize(list) );
}

/* variance of all vectors in a list. */
double zVecListVar(const zVecList *list, zVec mean)
{
  zVecListCell *vc;
  double var = 0;

  if( !zVecSizeEqual( mean, zListHead(list)->data ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return HUGE_VAL;
  }
  zListForEach( list, vc )
    var += zVecSqrDist( vc->data, mean );
  return var / zListSize(list);
}

/* mean and variance of all vectors in a list. */
double zVecListMeanVar(const zVecList *list, zVec mean)
{
  zVecListMean( list, mean );
  return zVecListVar( list, mean );
}

/* variance-covariance matrix of all vectors in a list. */
zMat zVecListCov(const zVecList *list, zVec mean, zMat cov)
{
  zVec v;
  zVecListCell *vc;
  int s;

  s = zVecSizeNC( zListHead(list)->data );
  if( zVecSizeNC(mean) != s ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return NULL;
  }
  if( !zMatIsSqr(cov) || zMatRowSizeNC(cov) != s ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  if( ( v = zVecAlloc( s ) ) == NULL ){
    ZALLOCERROR();
    return NULL;
  }
  zMatZero( cov );
  zListForEach( list, vc ){
    zVecSub( vc->data, mean, v );
    zMatAddDyadNC( cov, v, v );
  }
  zVecFree( v );
  return zMatDivDRC( cov, zListSize(list) );
}

/* mean and variance-covariance matrix of all vectors in a list. */
zMat zVecListMeanCov(const zVecList *list, zVec mean, zMat cov)
{
  zVecListMean( list, mean );
  return zVecListCov( list, mean, cov );
}

/* principal component analysis. */
int zVecListPCA(const zVecList *points, double cr, zVec mean, zVec score, zMat loading)
{
  zMat cov;
  int s, n;
  double score_th, score_sum;

  s = zVecSize( zListHead(points)->data );
  if( !( cov = zMatAllocSqr(s) ) ){
    ZALLOCERROR();
    n = -1;
    goto TERMINATE;
  }
  if( !zVecListMeanCov( points, mean, cov ) ||
      !zMatSymEigJacobi( cov, score, loading ) ){
    n = -1;
    goto TERMINATE;
  }
  /* contribution ratio */
  score_th = cr * zVecElemSum( score );
  for( score_sum=0, n=0; n<s; n++ ){
    score_sum += zVecElemNC(score,n);
    if( score_sum > score_th ) break;
  }

 TERMINATE:
  zMatFree( cov );
  return n;
}

/* generate vectors generate vectors from normal distribution. */
int zVecListGenRandND(zVecList *vl, int n, zVec mean, zMat cov)
{
  int i = 0, j;
  zVec v, vo;
  zMat b = NULL;
  zIndex idx = NULL;

  if( !zMatIsSqr( cov ) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return 0;
  }
  if( !zMatColVecSizeEqual( cov, mean ) ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH_VEC );
    return 0;
  }
  vo = zVecAlloc( zVecSizeNC(mean) );
  v = zVecAlloc( zVecSizeNC(mean) );
  if( !vo || !v || !zMatDecompCholeskyAlloc( cov, &b, &idx ) ){
    ZALLOCERROR();
    goto TERMINATE;
  }
  zListInit( vl );
  for( i=0; i<n; i++ ){
    for( j=0; j<zVecSizeNC(vo); j++ )
      zVecSetElemNC( vo, j, zRandND0( NULL ) );
    zMulMatVecNC( b, vo, v );
    zVecAddNCDRC( v, mean );
    if( !zVecListInsertHead( vl, v ) ) break;
  }
 TERMINATE:
  zVecFree( vo );
  zVecFree( v );
  zMatFree( b );
  zIndexFree( idx );
  return i;
}
