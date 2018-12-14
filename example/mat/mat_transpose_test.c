#include <zm/zm_mat.h>

int main(void)
{
  double marray[] = {
    1, 2, 3, 4, 5, 6, 7, 8,
    9,10,11,12,13,14,15,16,
    17,18,19,20,21,22,23,24,
    25,26,27,28,29,30,31,32,
  };
  const int r = 4, c = 8;
  zMat m, mt, mp;

  m = zMatCloneArray( marray, r, c );
  mt = zMatAlloc( zMatColSize(m), zMatRowSize(m) );
  printf( "m  : " ); zMatWrite( m );
  printf( "m^T: " ); zMatWrite( zMatT( m, mt ) );
  printf( "m  : " ); zMatWrite( zMatTDST( m ) );
  zMatFree( mt );

  mt = zMatTClone( m );
  printf( "m^T: " ); zMatWrite( mt );

  mp = zMatAlloc( 3, 2 );
  zMatGet( m, 3, 1, mp );
  printf( "(3,1)-(5,2): " ); zMatWrite( mp );
  zMatPut( mt, 0, 0, mp );
  printf( "put into m^T at (0,0): " ); zMatWrite( mt );
  zMatTGet( m, 3, 1, mp );
  printf( "(3,1)-(4,3): " ); zMatWrite( mp );
  zMatTPut( mt, 1, 4, mp );
  printf( "put into m^T at (1,4): " ); zMatWrite( mt );

  zMatFree( mt );
  zMatFree( m );
  return 0;
}
