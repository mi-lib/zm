#include <zm/zm_data.h>
#include <zm/zm_le.h>

/* TEST 1: plane, 2: parabola */
#define TEST 2

#if TEST == 1
/* test: plane */

int q_size = 2;

double y_test(double x)
{
  /* plane: y = 5 ( x - 2 ) */
  return 5 * ( x - 2 );
}

double error_case(zVec q, zVec sample, void *util)
{
  return zVecInnerProd( q, sample ) - 1;
}

zVec fit_case(zVec q, zVecList *list, void *util)
{
  zMat m;
  zVec p;
  zVecListCell *sp;

  m = zMatAllocSqr( zVecSizeNC(q) );
  p = zVecAlloc( zVecSizeNC(q) );
  if( !m || !p ){
    q = NULL;
    goto TERMINATE;
  }
  zListForEach( list, sp ){
    zMatAddDyadNC( m, sp->data, sp->data );
    zVecAddDRC( p, sp->data );
  }
  zLESolveGauss( m, p, q );

 TERMINATE:
  zMatFree( m );
  zVecFree( p );
  return q;
}

void print_case(zVec q)
{
  printf( "plot [-10:10] 's' u 1:2 w p lt 2, (%.10g)*x+(%.10g) w l lt 7\n", -zVecElemNC(q,0)/zVecElemNC(q,1), 1.0/zVecElemNC(q,1) );
}

#else

/* test: parabola */

int q_size = 3;

double y_test(double x)
{
  /* parabola: y = ( x - 0.5 )^2 - 1 = x^2 - x - 0.75 */
  return zSqr( x - 0.5 ) - 1;
}

double error_case(zVec q, zVec sample, void *util)
{
  return zVecElemNC(q,0)*zSqr(zVecElemNC(sample,0)) +
         zVecElemNC(q,1)*zVecElemNC(sample,0) +
         zVecElemNC(q,2) - zVecElemNC(sample,1);
}

zVec fit_case(zVec q, zVecList *list, void *util)
{
  zMat m;
  zVec p, pi;
  zVecListCell *sp;

  m = zMatAllocSqr( zVecSizeNC(q) );
  p = zVecAlloc( zVecSizeNC(q) );
  pi = zVecAlloc( zVecSizeNC(q) );
  if( !m || !p || !pi ){
    q = NULL;
    goto TERMINATE;
  }
  zListForEach( list, sp ){
    zVecSetElemList( pi, zSqr(zVecElemNC(sp->data,0)), zVecElemNC(sp->data,0), 1.0 );
    zMatAddDyadNC( m, pi, pi );
    zVecCatDRC( p, zVecElemNC(sp->data,1), pi );
  }
  zLESolveGauss( m, p, q );

 TERMINATE:
  zMatFree( m );
  zVecFree( p );
  zVecFree( pi );
  return q;
}

void print_case(zVec q)
{
  printf( "plot [-10:10] 's' u 1:2 w p lt 2, (%.10g)*x*x+(%.10g)*x+(%.10g) w l lt 7\n", zVecElemNC(q,0), zVecElemNC(q,1), zVecElemNC(q,2) );
}

#endif

void sample_list(zVecList *sample, int n, double r, double nl)
{
  zVec v;
  double x, y;
  int i;
  FILE *fp;

  zListInit( sample );
  fp = fopen( "s", "w" );
  for( i=0; i<n; i++ ){
    v = zVecAlloc( 2 );
    y = y_test( ( x = zRandF( -10, 10 ) ) ) + zRandF(-nl,nl);
    /* an outlier model with a probability of r */
    if( zRandF(0,1) < r ) y += zRandF(-50,50);
    zVecSetElemList( v, x, y );
    zVecValueFPrint( fp, v );
    zVecListInsertHead( sample, v );
  }
  fclose( fp );
}

#define N  100
#define R    0.5
#define NS  20
#define NT  50
#define TH   0.5
#define NL   5.0

int main(int argc, char *argv[])
{
  zVecList sample;
  zVec q;

  zRandInit();
  sample_list( &sample, N, R, NL );
  q = zVecAlloc( q_size );
#if 0
  zRANSAC( q, &sample, fit_case, error_case, NULL, NS, NT, TH );
#else
  zRANSACAuto( q, &sample, fit_case, error_case, NULL, R, NL );
#endif
  print_case( q );
  zVecFree( q );
  zVecListDestroy( &sample );
  return 0;
}
