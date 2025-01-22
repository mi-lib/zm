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
  return EXIT_SUCCESS;
}
