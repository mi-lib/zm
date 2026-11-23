#include <zm/zm.h>
#include <float.h>

/* check if an expansion has nonoverlapping nonzero components sorted in increasing order of magnitude,
 * or a single zero component. */
bool check_expansion(int n, const double e[])
{
  double m;
  int i, ex;

  if( n < 1 ) return false;
  if( n == 1 ) return true;
  for( i=0; i<n; i++ )
    if( e[i] == 0 ) return false;
  for( i=1; i<n; i++ ){
    /* the lowest nonzero bit of e[i] is 2^ex */
    m = ldexp( frexp( e[i], &ex ), DBL_MANT_DIG );
    for( ex-=DBL_MANT_DIG; fmod( m, 2 ) == 0; m/=2, ex++ );
    if( fabs( e[i-1] ) >= ldexp( 1.0, ex ) ) return false;
  }
  return true;
}

/* check if an expansion is zero. */
bool is_zero_expansion(int n, const double e[])
{
  return n == 1 && e[0] == 0;
}

int sgn_long(long x){ return x > 0 ? 1 : ( x < 0 ? -1 : 0 ); }

/* random integer points: random, coplanar, perturbed from coplanar ones and collinear ones. */
void create_points_int(long p[4][3], int mode)
{
  long s1, s2;
  int i, j;

  for( j=0; j<3; j++ )
    for( i=0; i<3; i++ ) p[j][i] = zRandI( -10000, 10000 );
  if( mode == 0 ){
    for( i=0; i<3; i++ ) p[3][i] = zRandI( -10000, 10000 );
    return;
  }
  s1 = zRandI( -3, 3 );
  s2 = zRandI( -3, 3 );
  for( i=0; i<3; i++ ) p[3][i] = p[0][i] + s1*(p[1][i]-p[0][i]) + s2*(p[2][i]-p[0][i]);
  if( mode == 2 ) p[3][zRandI(0,2)] += zRandI(0,1) ? 1 : -1;
  if( mode == 3 ) /* collinear */
    for( i=0; i<3; i++ ) p[2][i] = p[0][i] + s1*(p[1][i]-p[0][i]);
}

/* the exact orientation of four integer points by the determinant of differences. */
long orient3d_int(long p[4][3])
{
  long a[3], b[3], c[3];
  int i;

  for( i=0; i<3; i++ ){
    a[i] = p[1][i] - p[0][i];
    b[i] = p[2][i] - p[0][i];
    c[i] = p[3][i] - p[0][i];
  }
  return a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0]) + a[2]*(b[0]*c[1]-b[1]*c[0]);
}

/* the orientation of four points as an expansion of det[ p1 1; p2 1; p3 1; p0 1 ], which equals to
 * the determinant of differences, expanded along the last column. */
int orient3d_expansion(double p[4][3], double e[])
{
  const int r[4][3] = { { 2, 3, 0 }, { 1, 3, 0 }, { 1, 2, 0 }, { 1, 2, 3 } };
  const double s[4] = { -1, 1, -1, 1 };
  const int perm[6][3] = { { 0, 1, 2 }, { 0, 2, 1 }, { 1, 0, 2 }, { 1, 2, 0 }, { 2, 0, 1 }, { 2, 1, 0 } };
  const double ps[6] = { 1, -1, -1, 1, 1, -1 };
  int i, k, n = 0;

  for( k=0; k<4; k++ )
    for( i=0; i<6; i++ )
      n = zAddExpansionTripleProd( n, e, s[k]*ps[i]*p[r[k][0]][perm[i][0]], p[r[k][1]][perm[i][1]], p[r[k][2]][perm[i][2]] );
  return n;
}

/* arithmetics of expansions, which are checked by the exact signs of determinants of integer points
 * scaled by powers of two and translated by integers, which are computed exactly. */
void assert_expansion(void)
{
  const int testnum = 100000;
  long pi[4][3], offset;
  double p[4][3], e[96+4], h[4*(96+4)], scale, b, p0, p1;
  long det;
  int i, j, k, n, m;
  bool result_grow = true, result_scale = true, result_triple = true, result_sgn = true;

  for( i=0; i<testnum; i++ ){
    create_points_int( pi, i % 4 );
    scale = ldexp( 1.0, zRandI(-60,60) );
    offset = i % 3 == 0 ? (long)zRandI(-1000000,1000000) * 4096 : 0;
    for( j=0; j<4; j++ )
      for( k=0; k<3; k++ ) p[j][k] = (double)( pi[j][k] + offset ) * scale;
    n = orient3d_expansion( p, e );
    det = orient3d_int( pi );
    if( !check_expansion( n, e ) ) result_triple = false;
    if( zSgnExpansion( n, e ) != sgn_long( det ) ) result_sgn = result_triple = false;
    /* the value of the determinant, which is exact if |det| < 2^53 */
    if( labs( det ) < ( 1L << DBL_MANT_DIG ) ){
      m = zGrowExpansion( n, e, -(double)det * scale*scale*scale, h );
      if( !is_zero_expansion( m, h ) ) result_triple = false;
    }
    /* scaling */
    b = zRandF(-1,1) * pow( 10, zRandF(-10,10) );
    m = zScaleExpansion( n, e, b, h );
    if( !check_expansion( m, h ) || zSgnExpansion( m, h ) != zSgnExpansion( n, e ) * zSgn( b ) ) result_scale = false;
    /* cancellation of e b by subtracting each e[j] b, which is exactly represented by two values */
    for( j=0; j<n; j++ ){
      p1 = zMul( e[j], b, &p0 );
      m = zGrowExpansion( m, h, -p1, h );
      if( !check_expansion( m, h ) ) result_grow = false;
      m = zGrowExpansion( m, h, -p0, h );
      if( !check_expansion( m, h ) ) result_grow = false;
    }
    if( !is_zero_expansion( m, h ) ) result_scale = false;
  }
  /* empty and zero expansions, and an expression as an argument */
  e[0] = 0; e[1] = -1; e[2] = 2;
  if( zSgnExpansion( 0, e ) != 0 || zSgnExpansion( 1, e ) != 0 ||
      zSgnExpansion( 1, e+1 ) != -1 || zSgnExpansion( 1, e+2 ) != 1 ) result_sgn = false;
  zAssert( zSgnExpansion, result_sgn );
  zAssert( zGrowExpansion, result_grow );
  zAssert( zScaleExpansion, result_scale );
  zAssert( zAddExpansionTripleProd, result_triple );
}

int main(void)
{
  zRandInit();
  assert_expansion();
  return EXIT_SUCCESS;
}
