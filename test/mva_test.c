#include <zm/zm.h>
#include <float.h>

#define TOL 10.0 /* tolerance in the unit of n eps */

/* check if the scores are in descending order, the loading is orthonormal, and cov loading = loading diag(score). */
bool check_pca(zMat cov, zVec score, zMat loading)
{
  int i, n;
  zMat ltl;
  zVec v, e;
  double norm;
  bool result = true;

  n = zVecSizeNC(score);
  for( i=1; i<n; i++ )
    if( zVecElemNC(score,i-1) < zVecElemNC(score,i) ){
      eprintf( "scores not in descending order: %g < %g\n", zVecElemNC(score,i-1), zVecElemNC(score,i) );
      result = false;
    }
  ltl = zMatAllocSqr( n );
  v = zVecAlloc( n );
  e = zVecAlloc( n );
  norm = zMatNorm( cov );
  for( i=0; i<n; i++ ){
    zMatGetCol( loading, i, v );
    zMulMatVec( cov, v, e );
    zVecCatDRC( e, -zVecElemNC(score,i), v );
    if( !( zVecNorm( e ) <= TOL * n * DBL_EPSILON * norm ) ){
      eprintf( "principal component #%d: relative residual = %g\n", i, zVecNorm(e) / norm );
      result = false;
    }
  }
  zMulMatTMat( loading, loading, ltl );
  for( i=0; i<n; i++ ) zMatElemNC(ltl,i,i) -= 1;
  if( !( zMatNorm( ltl ) <= TOL * n * DBL_EPSILON ) ){
    eprintf( "loading not orthonormal: |L^T L - I| = %g\n", zMatNorm(ltl) );
    result = false;
  }
  zMatFree( ltl );
  zVecFree( v );
  zVecFree( e );
  return result;
}

/* the smallest number of components, the cumulative contribution ratio of which is not less than cr. */
int pca_num(zVec score, double cr)
{
  int n;
  double sum = 0, total;

  total = zVecElemSum( score );
  for( n=0; n<zVecSizeNC(score); n++ ){
    if( sum >= cr * total ) break;
    sum += zVecElemNC(score,n);
  }
  return n;
}

/* points on the coordinate axes, the variances along which are 1/3, 4/3 and 3 in ascending order. */
void assert_pca_axis(void)
{
  zVecList points;
  zVec mean, score, p;
  zMat loading, cov;
  double expected[] = { 3.0, 4.0/3, 1.0/3 };
  double cr[] =       { 0.0, 0.5, 0.9, 0.95, 1.0 };
  int num[] =         { 0,   1,   2,   3,    3 };
  int i, k, n;
  bool result = true;

  zListInit( &points );
  p = zVecAlloc( 3 );
  for( k=0; k<3; k++ )
    for( i=-1; i<=1; i+=2 ){
      zVecZero( p );
      zVecSetElemNC( p, k, i*(k+1) );
      zVecListInsertHead( &points, p ); /* a clone of p is inserted */
    }
  zVecFree( p );
  mean = zVecAlloc( 3 );
  score = zVecAlloc( 3 );
  loading = zMatAllocSqr( 3 );
  cov = zMatAllocSqr( 3 );
  zVecListMeanCov( &points, mean, cov );
  for( i=0; i<5; i++ ){
    if( ( n = zVecListPCA( &points, cr[i], mean, score, loading ) ) != num[i] ){
      eprintf( "cr = %g: number of principal components = %d (expected %d)\n", cr[i], n, num[i] );
      result = false;
    }
    if( !check_pca( cov, score, loading ) ) result = false;
    for( k=0; k<3; k++ ){
      if( !( fabs( zVecElemNC(score,k) - expected[k] ) <= TOL * DBL_EPSILON ) ){
        eprintf( "score #%d = %.17g (expected %.17g)\n", k, zVecElemNC(score,k), expected[k] );
        result = false;
      }
      /* the k-th principal axis is the (2-k)-th coordinate axis */
      if( !( fabs( fabs( zMatElemNC(loading,2-k,k) ) - 1 ) <= TOL * DBL_EPSILON ) ){
        eprintf( "principal axis #%d is not the coordinate axis #%d\n", k, 2-k );
        result = false;
      }
    }
  }
  zVecListDestroy( &points );
  zVecFree( mean );
  zVecFree( score );
  zMatFree( loading );
  zMatFree( cov );
  zAssert( zVecListPCA (points on axes), result );
}

#define N_TEST 100
#define N_POINT 200
void assert_pca_random(void)
{
  zVecList points;
  zVec mean, score, p;
  zMat loading, cov;
  int i, j, k, dim, n;
  double cr;
  bool result = true;

  for( k=0; k<N_TEST; k++ ){
    dim = zRandI( 1, 10 );
    mean = zVecAlloc( dim );
    score = zVecAlloc( dim );
    loading = zMatAllocSqr( dim );
    cov = zMatAllocSqr( dim );
    zListInit( &points );
    p = zVecAlloc( dim );
    for( i=0; i<N_POINT; i++ ){ /* components with various scales */
      for( j=0; j<dim; j++ ) zVecSetElemNC( p, j, zRandF(-1,1) * pow( 10, (j*7)%5 - 2 ) );
      zVecListInsertHead( &points, p ); /* a clone of p is inserted */
    }
    zVecFree( p );
    cr = zRandF( 0, 1 );
    n = zVecListPCA( &points, cr, mean, score, loading );
    zVecListMeanCov( &points, mean, cov );
    if( !check_pca( cov, score, loading ) ) result = false;
    if( n != pca_num( score, cr ) ){
      eprintf( "cr = %g: number of principal components = %d (expected %d)\n", cr, n, pca_num( score, cr ) );
      result = false;
    }
    zVecListDestroy( &points );
    zVecFree( mean );
    zVecFree( score );
    zMatFree( loading );
    zMatFree( cov );
  }
  zAssert( zVecListPCA (random), result );
}

void assert_pca_invalid(void)
{
  zVecList points;
  zVec mean, score, p;
  zMat loading;
  bool result = true;

  zListInit( &points );
  mean = zVecAlloc( 2 );
  score = zVecAlloc( 2 );
  loading = zMatAllocSqr( 2 );
  if( zVecListPCA( &points, 0.9, mean, score, loading ) != -1 ) result = false; /* empty set */
  p = zVecCreateList( 3, 1.0, 2.0, 3.0 );
  zVecListInsertHead( &points, p );
  zVecFree( p );
  if( zVecListPCA( &points, 0.9, mean, score, loading ) != -1 ) result = false; /* size mismatch */
  zVecListDestroy( &points );
  zVecFree( mean );
  zVecFree( score );
  zMatFree( loading );
  zAssert( zVecListPCA (invalid input), result );
}

int main(int argc, char *argv[])
{
  zRandInit();
  assert_pca_axis();
  assert_pca_random();
  assert_pca_invalid();
  return EXIT_SUCCESS;
}
