#include <zm/zm.h>

void assert_permut(void)
{
  zAssert( zPermut (zFactorial),
    zFactorial( 3 ) == 6 &&
    zFactorial( 4 ) == 24 &&
    zFactorial( 5 ) == 120 &&
    zFactorial( 6 ) == 720 );
}

void assert_combi(void)
{
  zAssert( zCombi,
    zCombi( 2, 0 ) == 1 &&
    zCombi( 2, 1 ) == 2 &&
    zCombi( 2, 2 ) == 1 &&
    zCombi( 3, 0 ) == 1 &&
    zCombi( 3, 1 ) == 3 &&
    zCombi( 3, 2 ) == 3 &&
    zCombi( 3, 3 ) == 1 &&
    zCombi( 4, 0 ) == 1 &&
    zCombi( 4, 1 ) == 4 &&
    zCombi( 4, 2 ) == 6 &&
    zCombi( 4, 3 ) == 4 &&
    zCombi( 4, 4 ) == 1 &&
    zCombi( 5, 0 ) == 1 &&
    zCombi( 5, 1 ) == 5 &&
    zCombi( 5, 2 ) == 10 &&
    zCombi( 5, 3 ) == 10 &&
    zCombi( 5, 4 ) == 5 &&
    zCombi( 5, 5 ) == 1 );
}

void assert_combi_recursive(void)
{
  int n, i, k;
  const int testnum = 100;
  bool result = true;

  for( k=0; k<testnum; k++ ){
    n = zRandI( 1, 20 );
    i = zRandI( 0, n );
    if( !zEqual( zCombi( n, i ), zCombiRecursive( n, i ), zTOL ) ) result = false;
  }
  zAssert( zCombi & zCombiRecursive, result );
}

#define COMBI_MAX_N 100
void assert_combi_series(void)
{
  int n, i;
  double c1[COMBI_MAX_N], c2[COMBI_MAX_N];
  bool result = true;

  n = 20;
  for( i=0; i<=n; i++ )
    c1[i] = zCombi( n, i );
  zCombiSeries( n, COMBI_MAX_N, c2 );
  for( i=0; i<=n; i++ )
    if( !zEqual( c1[i], c2[i], zTOL ) ) result = false;
  zAssert( zCombiSeries, result );
}

int main(void)
{
  double val;

  zRandInit();

  val = zRandF( -10, 10 );
  zAssert( zDeg2Rad, zIsTiny( zRad2Deg( zDeg2Rad(val) ) - val ) );
  zAssert( zRad2Deg, zIsTiny( zDeg2Rad( zRad2Deg(val) ) - val ) );

  zAssert( zIsEven, zIsEven( 0 ) && zIsEven( 2 ) && zIsEven( -2 ) && !zIsEven( 1 ) );
  zAssert( zIsOdd, !zIsOdd( 0 ) && zIsOdd( 1 ) && zIsOdd( -1 ) );

#if 0
  zAssert( zIsSgnOpp, zIsSgnOpp( 1, -1 ) && !zIsSgnOpp( 1, 0 ) && !zIsSgnOpp( -1, 0 ) );
#endif

  zAssert( zPhaseNormalize, zIsTiny( zPhaseNormalize(3*zPI) - zPI ) &&
                            zIsTiny( zPhaseNormalize(-zPI-zTOL) - zPI ) );

  zAssert( zCeil, zCeil(0.5) == 1.0 && zCeil(-0.5) == -1.0 );
  zAssert( zRound, zRound( 1.499999 ) == 1.0 && zRound( 1.50000001 ) == 2.0 &&
                   zRound(-1.499999 ) ==-1.0 && zRound(-1.50000001 ) ==-2.0 );
  zAssert( zFruct, zFruct( zPI, 2 ) == zPI - 2 && zFruct( zPIx2, 2 ) == zPIx2 - 6 );

  val = zRandF( 1, 10 );
  zAssert( zCbrt, zIsTiny( zCube( zCbrt( val ) ) - val ) );
  val = zRandF( zTOL, 1.0e10 );
  zAssert( zLog, zIsTiny( zLog(2,val) - log2(val) ) &&
                 zIsTiny( zLog(10,val) - log10(val) ) );

  assert_permut();
  assert_combi();
  assert_combi_recursive();
  assert_combi_series();
  return EXIT_SUCCESS;
}
