#include <zm/zm_le.h>

void set_weight(zVec w, double v)
{
  register int i;

  for( i=0; i<zVecSize(w); i++ ) zVecSetElem( w, i, v );
}

int main(void)
{
  double v;
  zMat a;
  zVec b, b_, w1, w2, x, ref;

  a = zMatCreateList( 2, 3,
    2.0, 1.0, 1.0,
    1.0, 2.0, 1.0 );
  b = zVecCreateList( 2, 4.0, 4.0 );
  b_ = zVecAlloc( 2 );
  w1 = zVecCreateList( 2, 1.0, 1.0 );
  w2 = zVecAlloc( 3 );
  x = zVecAlloc( 3 );
  ref = zVecCreateList( 3, 0.0, 1.0, 2.0 );

  zMatWrite( a );
  zVecWrite( b );
  zVecWrite( w1 );
  for( v=10000; v>0.0000001; v*=0.1 ){
    set_weight( w2, v );
    zLESolveRSR( a, b, w1, w2, ref, x );
    zVecWrite( x );
    zMulMatVec( a, x, b_ );
    zIndent( 2 );
    zVecWrite( b_ );
  }
  return 0;
}
