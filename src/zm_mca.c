/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_mca - multiple classification analysis.
 */

#include <zm/zm_mca.h>

/* sum up all vectors in a list. */
zVec zVecListSum(zVecList *list, zVec sum)
{
  zVecListCell *vc;

  zVecZero( sum );
  zListForEach( list, vc )
    zVecAddDRC( sum, vc->data );
  return sum;
}

/* average of all vectors in a list. */
zVec zVecListAve(zVecList *list, zVec ave)
{
  zVecListSum( list, ave );
  return zVecDivDRC( ave, zListNum(list) );
}

/* variance of all vectors in a list. */
double zVecListVar(zVecList *list, zVec ave)
{
  zVecListCell *vc;
  double var = 0;

  if( !zVecSizeIsEqual( ave, zListHead(list)->data ) ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return HUGE_VAL;
  }
  zListForEach( list, vc )
    var += zVecSqrDist( vc->data, ave );
  return var / zListNum(list);
}

/* average and variance of all vectors in a list. */
double zVecListAveVar(zVecList *list, zVec ave)
{
  zVecListAve( list, ave );
  return zVecListVar( list, ave );
}

/* variance-covariance matrix of all vectors in a list. */
zMat zVecListCov(zVecList *list, zVec ave, zMat cov)
{
  zVec v;
  zVecListCell *vc;
  int s;

  s = zVecSize( zListHead(list)->data );
  if( zVecSizeNC(ave) != s ){
    ZRUNERROR( ZM_ERR_SIZMIS_MATVEC );
    return NULL;
  }
  if( !zMatIsSqr(cov) || zMatRowSizeNC(cov) != s ){
    ZRUNERROR( ZM_ERR_SIZMIS_MAT );
    return NULL;
  }
  if( ( v = zVecAlloc( s ) ) == NULL ){
    ZALLOCERROR();
    return NULL;
  }
  zMatZero( cov );
  zListForEach( list, vc ){
    zVecSub( vc->data, ave, v );
    zMatAddDyadNC( cov, v, v );
  }
  zVecFree( v );
  return zMatDivDRC( cov, zListNum(list) );
}

/* average and variance-covariance matrix of all vectors in a list. */
zMat zVecListAveCov(zVecList *list, zVec ave, zMat cov)
{
  zVecListAve( list, ave );
  return zVecListCov( list, ave, cov );
}

/* principal component analysis. */
int zPCA(zVecList *points, double cr, zVec ave, zVec score, zMat loading)
{
  zMat cov;
  int s, n;
  double score_th, score_sum;

  s = zVecSize( zListHead(points)->data );
  if( ( cov = zMatAllocSqr(s) ) == NULL ){
    ZALLOCERROR();
    n = -1;
    goto TERMINATE;
  }
  if( zVecListAveCov( points, ave, cov ) == NULL ||
      zEigSymJacobi( cov, score, loading ) == NULL ){
    n = -1;
    goto TERMINATE;
  }
  /* contribution ratio */
  score_th = cr * zVecSum( score );
  for( score_sum=0, n=0; n<s; n++ ){
    score_sum += zVecElem(score,n);
    if( score_sum > score_th ) break;
  }

 TERMINATE:
  zMatFree( cov );
  return n;
}
