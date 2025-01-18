#include <zm/zm_ip.h>

bool alloc_queue(zSeq *seq)
{
  zVec v;
  int i;

  zSeqInit( seq );
  for( i=0; i<3; i++ ){
    if( !( v = zVecAlloc( 1 ) ) ) return false;
    zSeqEnqueue( seq, v, 0 /* dummy */ );
  }
  return true;
}

void approx_by_lag(zIP *ip, zSeq *seq, double (*f)(double), double t, double eps)
{
  zSeqCell *c;

  /* sample 1 */
  c = zListHead( seq );
  c->data.dt = 0;
  zVecSetElem( c->data.v, 0, f( t-eps ) );
  /* sample 2 */
  c = zListCellPrev( c );
  c->data.dt = eps;
  zVecSetElem( c->data.v, 0, f( t ) );
  /* sample 3 */
  c = zListCellPrev( c );
  c->data.dt = eps;
  zVecSetElem( c->data.v, 0, f( t+eps ) );
  /* approximation by Lagrange interpolation */
  zIPCreateLagrange( ip, seq );
}



#define T    1.0
#define DT   0.1
#define EPST 0.01

double test_function(double t)
{
  return 2*t*(t-0.5)*(t-1) - 2 * (t-0.5) + sin(2*zPIx2*t) + exp(2-2*t);
}

double test_function_dif(double t)
{
  return 2*((t-0.5)*(t-1)+t*(t-1)+t*(t-0.5)) - 2 + 2*zPIx2*cos(2*zPIx2*t) - 2* exp(2-2*t);
}

double test_function_dif2(double t)
{
  return 6*(2*t-1) - zSqr(2*zPIx2)*sin(2*zPIx2*t) + 4* exp(2-2*t);
}

int main(int argc, char *argv[])
{
  zSeq seq;
  zIP ip;
  zVec v;
  int i, step;
  double t;

  alloc_queue( &seq );
  v = zVecAlloc( 1 );
  step = T / DT;
  for( i=0; i<=step; i++ ){
    t = i * DT;
    printf( "%g %g %g %g", t, test_function( t ), test_function_dif( t ), test_function_dif2( t ) );

    approx_by_lag( &ip, &seq, test_function, t, EPST );
    printf( " %g",   zVecElemNC( zIPVec( &ip, EPST, v ), 0 ) );
    printf( " %g",   zVecElemNC( zIPVel( &ip, EPST, v ), 0 ) );
    printf( " %g\n", zVecElemNC( zIPAcc( &ip, EPST, v ), 0 ) );
    zIPDestroy( &ip );
  }
  zVecFree( v );
  zSeqFree( &seq );
  return 0;  
}
