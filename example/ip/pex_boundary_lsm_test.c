/* least square method */
#include <zm/zm_ip.h>

bool zPexIPCreateBounderyLSM(zPexIP *pc, double term, double x1, double v1, double a1, double x2, double v2, double a2, int dim, zVec t, zVec x)
{
  register uint i, j, k, n, m;
  zMat a;
  zVec b, v;
  bool result = true;

  if( !zVecSizeIsEqual( t, x ) ){
    ZRUNERROR( "size mismatch of sample point vector" );
    return false;
  }

  n = zVecSizeNC( t );
  if( !zPexIPCreate( pc, term, dim ) ) return false;
  m = dim + 1;
  a = zMatCreateSqr( m );
  b = zVecCreate( m );
  v = zVecCreate( 2 * m );
  if( !a || !b || !v ){
    result = false;
    goto TERMINATE;
  }

  for( i=0; i<n; i++ ){
    for( j=0; j<zVecSizeNC(v); j++ )
      zVecSetElem( v, j, pow( zVecElem(t,i)/term, j ) );
    for( j=3; j<m-3; j++ ){
      for( k=0; k<m; k++ )
        zMatElem(a,j,k) += zVecElem(v,j+k);
      zVecElem(b,j) += zVecElem(x,i) * zVecElem(v,j);
    }
  }
  zMatSetElem( a, 0, 0, 1 );
  zMatSetElem( a, 1, 1, 1/term );
  zMatSetElem( a, 2, 2, 2/zSqr(term) );
  for( i=0; i<m; i++ )
    zMatSetElem( a, m-3, i, 1 );
  for( i=1; i<m; i++ )
    zMatSetElem( a, m-2, i, i );
  for( i=2; i<m; i++ )
    zMatSetElem( a, m-1, i, i*(i-1) );
  zVecSetElem( b, 0, x1 );
  zVecSetElem( b, 1, v1 );
  zVecSetElem( b, 2, a1 );
  zVecSetElem( b, m-3, x2 );
  zVecSetElem( b, m-2, v2 );
  zVecSetElem( b, m-1, a2 );

  zLESolveGauss( a, b, pc->c );

 TERMINATE:
  zMatDestroy( a );
  zVecDestroy( b );
  zVecDestroy( v );
  return result;
}

int main(int argc, char *argv[])
{
  zPexIP pc;
  zVec tvec, xvec;
  double t;

  tvec = zVecCreateList( 3, 1.0, 2.0, 4.0 );
  xvec = zVecCreateList( 3, 4.0, 0.0, 2.0 );
  zPexIPCreateBounderyLSM( &pc, 5.0, 0, 0, 0, 1, 0, 0, 8, tvec, xvec );

  for( t=0; t<=zPexIPTerm(&pc); t+=0.01 )
    printf( "%f\n", zPexIPVal( &pc, t ) );

  zPexIPDestroy( &pc );
  return 0;
}
