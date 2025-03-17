#include <zm/zm.h>

void assert_permutation(void)
{
  zAssert( zPermutation (zFactorial),
    zFactorial( 3 ) == 6 &&
    zFactorial( 4 ) == 24 &&
    zFactorial( 5 ) == 120 &&
    zFactorial( 6 ) == 720 );
}

void assert_combination(void)
{
  zAssert( zCombination,
    zCombination( 2, 0 ) == 1 &&
    zCombination( 2, 1 ) == 2 &&
    zCombination( 2, 2 ) == 1 &&
    zCombination( 3, 0 ) == 1 &&
    zCombination( 3, 1 ) == 3 &&
    zCombination( 3, 2 ) == 3 &&
    zCombination( 3, 3 ) == 1 &&
    zCombination( 4, 0 ) == 1 &&
    zCombination( 4, 1 ) == 4 &&
    zCombination( 4, 2 ) == 6 &&
    zCombination( 4, 3 ) == 4 &&
    zCombination( 4, 4 ) == 1 &&
    zCombination( 5, 0 ) == 1 &&
    zCombination( 5, 1 ) == 5 &&
    zCombination( 5, 2 ) == 10 &&
    zCombination( 5, 3 ) == 10 &&
    zCombination( 5, 4 ) == 5 &&
    zCombination( 5, 5 ) == 1 );
}

void assert_combination_recursive(void)
{
  int n, i, k;
  const int testnum = 100;
  bool result = true;

  for( k=0; k<testnum; k++ ){
    n = zRandI( 1, 20 );
    i = zRandI( 0, n );
    if( !zEqual( zCombination( n, i ), zCombinationRecursive( n, i ), zTOL ) ) result = false;
  }
  zAssert( zCombination & zCombinationRecursive, result );
}

#define COMBI_MAX_N 100
void assert_combination_pascaltriangle(void)
{
  const int n = 8;
  int i, j;
  double c[COMBI_MAX_N];
  bool result = true;

  for( i=0; i<=n; i++ ){
    zCombinationSeries( i, COMBI_MAX_N, c );
    for( j=0; j<=i; j++ ){
      if( !zEqual( c[j], zCombination(i,j), zTOL ) ) result = false;
      if( !zEqual( c[j], zCombinationRecursive(i,j), zTOL ) ) result = false;
    }
  }
  zAssert( zCombinationSeries (pascal triangle), result );
}

int main(void)
{
  double val;

  zRandInit();

  val = zRandF( -10, 10 );
  zAssert( zDeg2Rad, zEqual( zRad2Deg( zDeg2Rad(val) ), val, zTOL ) );
  zAssert( zRad2Deg, zEqual( zDeg2Rad( zRad2Deg(val) ), val, zTOL ) );

  zAssert( zIsEven, zIsEven( 0 ) && zIsEven( 2 ) && zIsEven( -2 ) && !zIsEven( 1 ) );
  zAssert( zIsOdd, !zIsOdd( 0 ) && zIsOdd( 1 ) && zIsOdd( -1 ) );

#if 0
  zAssert( zIsSgnOpp, zIsSgnOpp( 1, -1 ) && !zIsSgnOpp( 1, 0 ) && !zIsSgnOpp( -1, 0 ) );
#endif

  zAssert( zPhaseNormalize,
    zEqual( zPhaseNormalize(3*zPI), zPI, zTOL ) &&
    zEqual( zPhaseNormalize(-zPI-zTOL), zPI, zTOL ) );

  zAssert( zCeil, zCeil(0.5) == 1.0 && zCeil(-0.5) == -1.0 );
  zAssert( zRound,
    zRound( 1.499999 ) == 1.0 && zRound( 1.50000001 ) == 2.0 &&
    zRound(-1.499999 ) ==-1.0 && zRound(-1.50000001 ) ==-2.0 );
  zAssert( zFruct, zFruct( zPI, 2 ) == zPI - 2 && zFruct( zPIx2, 2 ) == zPIx2 - 6 );

  val = zRandF( 1, 10 );
  zAssert( zCbrt, zIsTiny( zCube( zCbrt( val ) ) - val ) );
  val = zRandF( zTOL, 1.0e10 );
  zAssert( zLog,
    zEqual( zLog(2,val), log2(val), zTOL ) &&
    zEqual( zLog(10,val), log10(val), zTOL ) );

  assert_permutation();
  assert_combination();
  assert_combination_recursive();
  assert_combination_pascaltriangle();
  return EXIT_SUCCESS;
}
