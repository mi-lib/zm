#include <zm/zm_ip.h>

void assert_spline_cubicequation(void)
{
  zSeq seq;
  zIP ip;
  zVec v, a, b, c, d, v_cubic;
  const int test_num = 100;
  int point_num, i, k;
  double tp[] = { 0, 2, 3, 5, 6 };
  double vp[] = { 0, 3, 1,-2, 4 };
  double t;
  bool result = true;

  point_num = sizeof(tp) / sizeof(double);
  zSeqInit( &seq );
  for( t=0, i=0; i<point_num; i++ ){
    v = zVecCreateList( 1, vp[i] );
    zSeqEnqueue( &seq, v, tp[i]-t );
    t = tp[i];
  }
  /* cubic spline interpolator */
  v = zVecAlloc( 1 );
  zVecZero( v );
  zIPCreateSpline( &ip, &seq, ZSPLINE_FIX_EDGE, v, ZSPLINE_FIX_EDGE, v );
  a = zVecAlloc( 1 );
  b = zVecAlloc( 1 );
  c = zVecAlloc( 1 );
  d = zVecAlloc( 1 );
  v_cubic = zVecAlloc( 1 );
  for( k=0; k<test_num; k++ ){
    i = zRandI( 0, point_num-2 );
    zIPSplineCoeff( &ip, i, a, b, c, d );
    t = zRandF( tp[i], tp[i+1] );
    zIPVec( &ip, t, v );
    t -= tp[i];
    zVecZero( v_cubic );
    zVecCatNCDRC( v_cubic, t*t*t, a );
    zVecCatNCDRC( v_cubic, t*t,   b );
    zVecCatNCDRC( v_cubic, t,     c );
    zVecAddNCDRC( v_cubic,        d );
    if( !zVecEqual( v, v_cubic, zTOL ) ){
      eprintf( "k=%d, i=%d\n", k, i );
      zVecFPrint( stderr, v );
      zVecFPrint( stderr, v_cubic );
      result = false;
    }
  }
  zVecFreeAtOnce( 6, v, a, b, c, d, v_cubic );
  zIPDestroy( &ip );
  zSeqFree( &seq );
  zAssert( zIPSplineCoeff, result );
}

/* polynomial function fitting a sequence of points sampled from a known polynomial function. */
#define TOL_PEX_IP_LSM 1.0e-10
void assert_pex_ip_lsm(void)
{
  const int dim = 4, n = 30, test_num = 100;
  zPexIP pc;
  zVec c, t, x;
  double term;
  int i, k;
  bool result = true, result_fewpoints;

  c = zVecAlloc( dim + 1 );
  t = zVecAlloc( n );
  x = zVecAlloc( n );
  for( k=0; k<test_num; k++ ){
    term = zRandF( 0.1, 10 );
    zVecRandUniform( c, -1, 1 );
    for( i=0; i<n; i++ ){
      zVecSetElemNC( t, i, zRandF( 0, term ) );
      zVecSetElemNC( x, i, zPexVal( c, zVecElemNC(t,i)/term ) );
    }
    if( !zPexIPCreateLSM( &pc, term, dim, t, x ) ){
      result = false;
      continue;
    }
    if( zVecDist( pc.c, c ) > TOL_PEX_IP_LSM * zVecNorm( c ) ) result = false;
    zPexIPFree( &pc );
  }
  /* the number of points is smaller than the number of coefficients */
  zVecSetSize( t, dim );
  zVecSetSize( x, dim );
  zEchoOff();
  result_fewpoints = !zPexIPCreateLSM( &pc, 1.0, dim, t, x );
  zEchoOn();
  zVecFreeAtOnce( 3, c, t, x );
  zAssert( zPexIPCreateLSM, result );
  zAssert( zPexIPCreateLSM (too few points), result_fewpoints );
}

/* the solution of the equality-constrained least square problem by the Lagrange multipliers, namely,
 * [ V^T V C^T; C O ] [ c; -lambda ] = [ V^T x; d ], where V is the Vandermonde matrix of the points
 * and C c = d are the boundary conditions. */
void pex_ip_boundary_lsm_kkt(double term, double *bc, const zVec t, const zVec x, zVec c)
{
  int i, j, l, m;
  zMat k;
  zVec rhs, sol, v;

  m = zVecSizeNC(c);
  k = zMatAllocSqr( m + 6 );
  rhs = zVecAlloc( m + 6 );
  sol = zVecAlloc( m + 6 );
  v = zVecAlloc( m );
  for( i=0; i<zVecSizeNC(t); i++ ){
    for( j=0; j<m; j++ ) zVecSetElemNC( v, j, pow( zVecElemNC(t,i)/term, j ) );
    for( j=0; j<m; j++ ){
      for( l=0; l<m; l++ ) zMatElemNC(k,j,l) += zVecElemNC(v,j) * zVecElemNC(v,l);
      zVecElemNC(rhs,j) += zVecElemNC(v,j) * zVecElemNC(x,i);
    }
  }
  for( j=0; j<m; j++ ){
    double cc[6];
    cc[0] = j == 0 ? 1 : 0;
    cc[1] = j == 1 ? 1/term : 0;
    cc[2] = j == 2 ? 2/zSqr(term) : 0;
    cc[3] = 1;
    cc[4] = j / term;
    cc[5] = j * ( j - 1 ) / zSqr(term);
    for( l=0; l<6; l++ ){
      zMatSetElemNC( k, m+l, j, cc[l] );
      zMatSetElemNC( k, j, m+l, cc[l] );
    }
  }
  for( l=0; l<6; l++ ) zVecSetElemNC( rhs, m+l, bc[l] );
  zLESolveGauss( k, rhs, sol );
  for( j=0; j<m; j++ ) zVecSetElemNC( c, j, zVecElemNC(sol,j) );
  zMatFree( k );
  zVecFreeAtOnce( 3, rhs, sol, v );
}

/* polynomial function that satisfies boundary conditions and fits a sequence of points. */
#define TOL_PEX_IP_BOUNDARY 1.0e-10
#define TOL_PEX_IP_BOUNDARY_LSM 1.0e-8 /* relative to the squared residual by the Lagrange multipliers */
void assert_pex_ip_boundary_lsm(void)
{
  const int dim = 8, n = 30, test_num = 100;
  zPexIP pc, pc5;
  zVec t, x, c;
  double term, bc[6], err[6], r, r_ref;
  int i, k;
  bool result_bc = true, result_lsm = true, result_dim5, result_smalldim;

  t = zVecAlloc( n );
  x = zVecAlloc( n );
  c = zVecAlloc( dim + 1 );
  for( k=0; k<test_num; k++ ){
    term = zRandF( 0.1, 10 );
    for( i=0; i<6; i++ ) bc[i] = zRandF( -1, 1 );
    for( i=0; i<n; i++ ){
      zVecSetElemNC( t, i, zRandF( 0, term ) );
      zVecSetElemNC( x, i, sin( zPIx2 * zVecElemNC(t,i) / term ) + zRandF( -0.1, 0.1 ) );
    }
    if( !zPexIPCreateBoundaryLSM( &pc, term, bc[0], bc[1], bc[2], bc[3], bc[4], bc[5], dim, t, x ) ){
      result_bc = result_lsm = false;
      continue;
    }
    /* boundary conditions, where the derivatives are evaluated with respect to t/term */
    err[0] = zPexIPVal(&pc,0) - bc[0];
    err[1] = ( zPexIPVel(&pc,0) - bc[1] ) * term;
    err[2] = ( zPexIPAcc(&pc,0) - bc[2] ) * zSqr(term);
    err[3] = zPexIPVal(&pc,term) - bc[3];
    err[4] = ( zPexIPVel(&pc,term) - bc[4] ) * term;
    err[5] = ( zPexIPAcc(&pc,term) - bc[5] ) * zSqr(term);
    if( zRawVecNorm( err, 6 ) > TOL_PEX_IP_BOUNDARY * zVecNorm( pc.c ) ) result_bc = false;
    /* comparison of the squared residual with that of the solution by the Lagrange multipliers.
     * the coefficients are not compared since the solution by the normal equation is less accurate
     * than the residual, which is stationary at the optimum. */
    pex_ip_boundary_lsm_kkt( term, bc, t, x, c );
    for( r=r_ref=0, i=0; i<n; i++ ){
      r += zSqr( zPexIPVal( &pc, zVecElemNC(t,i) ) - zVecElemNC(x,i) );
      r_ref += zSqr( zPexVal( c, zVecElemNC(t,i)/term ) - zVecElemNC(x,i) );
    }
    if( r - r_ref > TOL_PEX_IP_BOUNDARY_LSM * r_ref ) result_lsm = false;
    zPexIPFree( &pc );
  }
  /* no free coefficient */
  result_dim5 = zPexIPCreateBoundaryLSM( &pc, 2.0, 0, 1, 2, 3, 4, 5, 5, t, x ) &&
                zPexIPCreateBoundary( &pc5, 2.0, 0, 1, 2, 3, 4, 5, NULL ) && zVecEqual( pc.c, pc5.c, 0 );
  zPexIPFree( &pc );
  zPexIPFree( &pc5 );
  zEchoOff();
  result_smalldim = !zPexIPCreateBoundaryLSM( &pc, 2.0, 0, 1, 2, 3, 4, 5, 4, t, x );
  zEchoOn();
  zVecFreeAtOnce( 3, t, x, c );
  zAssert( zPexIPCreateBoundaryLSM (boundary conditions), result_bc );
  zAssert( zPexIPCreateBoundaryLSM (least square solution), result_lsm );
  zAssert( zPexIPCreateBoundaryLSM (no free coefficient), result_dim5 );
  zAssert( zPexIPCreateBoundaryLSM (too small dimension), result_smalldim );
}

int main(void)
{
  zRandInit();
  assert_spline_cubicequation();
  assert_pex_ip_lsm();
  assert_pex_ip_boundary_lsm();
  return 0;
}
