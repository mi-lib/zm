#include <zm/zm_le.h>

#define N 100

void check(zVec a, zVec b, zVec c, zVec d, zVec ans)
{
  register int i;
  double err;

  for( i=0; i<zVecSizeNC(ans); i++ ){
    err = zVecElem(b,i)*zVecElem(ans,i);
    if( i > 0 ) err += zVecElem(a,i)*zVecElem(ans,i-1);
    if( i < zVecSizeNC(a)-1 ) err += zVecElem(c,i)*zVecElem(ans,i+1);
    err -= zVecElem(d,i);
    eprintf( "%d:err=%g\n", i, err );
    if( !zIsTiny(err) )
      eprintf( "might be a bug.\n" );
  }
}

int main(void)
{
  zVec a, b, c, d, ans;

  zRandInit();
  a = zVecAlloc( N );
  b = zVecAlloc( N );
  c = zVecAlloc( N );
  d = zVecAlloc( N );
  ans = zVecAlloc( N );
  zVecRandUniform( a, -10, 10 );
  zVecRandUniform( b, -10, 10 );
  zVecRandUniform( c, -10, 10 );
  zVecRandUniform( d, -10, 10 );
  zTridiagSolve( a, b, c, d, ans );
  check( a, b, c, d, ans );

  zVecFreeAO( 5, a, b, c, d, ans );
  return 0;
}
