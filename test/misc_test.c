#include <zm/zm.h>

int main(void)
{
  double val;

  zRandInit();

  val = zRandF( -10, 10 );
  zAssert( zDeg2Rad, zIsTiny( zRad2Deg( zDeg2Rad(val) ) - val ) );
  zAssert( zRad2Deg, zIsTiny( zDeg2Rad( zRad2Deg(val) ) - val ) );

  zAssert( zIsEven, zIsEven( 0 ) && zIsEven( 2 ) && zIsEven( -2 ) && !zIsEven( 1 ) );
  zAssert( zIsOdd, !zIsOdd( 0 ) && zIsOdd( 1 ) && zIsOdd( -1 ) );
  zAssert( zIsSgnOpp, zIsSgnOpp( 1, -1 ) && !zIsSgnOpp( 1, 0 ) && !zIsSgnOpp( -1, 0 ) );

  zAssert( zPhaseNormalize, zIsTiny( zPhaseNormalize(3*zPI) - zPI ) &&
                            zIsTiny( zPhaseNormalize(-zPI-zTOL) - zPI ) );
  zAssert( zCeil, zCeil(0.5) == 1.0 && zCeil(-0.5) == -1.0 );
  zAssert( zRound, zRound( 1.499999 ) == 1.0 && zRound( 1.50000001 ) == 2.0 &&
                   zRound(-1.499999 ) ==-1.0 && zRound(-1.50000001 ) ==-2.0 );
  zAssert( zFruct, zFruct( zPI, 2 ) == zPI - 2 && zFruct( zPIx2, 2 ) == zPIx2 - 6 );

  val = zRandF( 1, 10 );
  zAssert( zCbrt, zIsTiny( zCube( zCbrt( val ) ) - val ) );
  return EXIT_SUCCESS;
}
