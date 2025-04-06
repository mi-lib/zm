#include <zm/zm_mva.h>

/* plane */

int q_size_plane = 2;

double y_test_plane(double x)
{
  /* plane: y = 5 ( x - 2 ) */
  return 5 * ( x - 2 );
}

double error_case_plane(const zVec q, const zVec sample, void *util)
{
  return zVecInnerProd( q, sample ) - 1;
}

zVec fit_case_plane(zVec q, const zVecList *list, void *util)
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

/* parabola */

int q_size_parabola = 3;

double y_test_parabola(double x)
{
  /* parabola: y = ( x - 0.5 )^2 - 1 = x^2 - x - 0.75 */
  return zSqr( x - 0.5 ) - 1;
}

double error_case_parabola(const zVec q, const zVec sample, void *util)
{
  return zVecElemNC(q,0)*zSqr(zVecElemNC(sample,0)) +
         zVecElemNC(q,1)*zVecElemNC(sample,0) +
         zVecElemNC(q,2) - zVecElemNC(sample,1);
}

zVec fit_case_parabola(zVec q, const zVecList *list, void *util)
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

void print_case(const zVec q_plane, const zVec q_parabola)
{
  printf( "plot [-10:10] 's' u 1:2 w p lt 2, (%.10g)*x+(%.10g) w l lt 7, (%.10g)*x*x+(%.10g)*x+(%.10g) w l lt 7\n", -zVecElemNC(q_plane,0)/zVecElemNC(q_plane,1), 1.0/zVecElemNC(q_plane,1), zVecElemNC(q_parabola,0), zVecElemNC(q_parabola,1), zVecElemNC(q_parabola,2) );
}

void sample_list_case(FILE *fp, zVecList *sample, double (* y_test)(double), int n, double r, double nl)
{
  zVec v;
  double x, y;
  int i;

  for( i=0; i<n; i++ ){
    v = zVecAlloc( 2 );
    y = y_test( ( x = zRandF( -10, 10 ) ) ) + zRandF(-nl,nl);
    /* an outlier model with a probability of r */
    if( zRandF(0,1) < r ) y += zRandF(-50,50);
    zVecSetElemList( v, x, y );
    zVecValueFPrint( fp, v );
    zVecListInsertHead( sample, v );
  }
}

void sample_list(zVecList *sample, int n, double r, double nl)
{
  FILE *fp;

  zListInit( sample );
  fp = fopen( "s", "w" );
  sample_list_case( fp, sample, y_test_plane, n, r, nl );
  sample_list_case( fp, sample, y_test_parabola, n, r, nl );
  fclose( fp );
}

#define N  300
#define R    0.2
#define NL   5.0

int main(int argc, char *argv[])
{
  zVecList sample;
  zVecList inlier_plane, inlier_parabola;
  zVec q_plane, q_parabola;

  zRandInit();
  sample_list( &sample, N, R, NL );
  q_plane = zVecAlloc( q_size_plane );
  q_parabola = zVecAlloc( q_size_parabola );
  zRANSACSaveInlierAuto( q_plane, &sample, fit_case_plane, error_case_plane, NULL, R, NL, &inlier_plane );
  zRANSACSaveInlierAuto( q_parabola, &sample, fit_case_parabola, error_case_parabola, NULL, R, NL, &inlier_parabola );
  zListAppend( &sample, &inlier_plane );
  zRANSACSaveInlierAuto( q_plane, &sample, fit_case_plane, error_case_plane, NULL, R, NL, &inlier_plane );
  zListAppend( &sample, &inlier_parabola );
  zRANSACSaveInlierAuto( q_parabola, &sample, fit_case_parabola, error_case_parabola, NULL, R, NL, &inlier_parabola );
  print_case( q_plane, q_parabola );

  zVecFree( q_plane );
  zVecFree( q_parabola );
  zListAppend( &sample, &inlier_plane );
  zListAppend( &sample, &inlier_parabola );
  zVecListDestroy( &sample );
  return 0;
}
