#include <zm/zm_mat.h>

int main(void)
{
  zMat m, pm;
  double val[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  int sr = 1, sc = 2, dr = 0, dc = 0;

  m = zMatCloneArray( val, 3, 4 );
  pm = zMatAlloc( 2, 2 );

  printf( "get (%d,%d)-(%d,%d)\n", sr, sc, sr+zMatRowSize(pm), sc+zMatColSize(pm) );
  zMatGet( m, sr, sc, pm );
  zMatWrite( m );
  zMatWrite( pm );
  printf( "put (%d,%d)-(%d,%d)\n", dr, dc, dr+zMatRowSize(pm), dc+zMatColSize(pm) );
  zMatPut( m, dr, dc, pm );
  zMatWrite( m );
  return 0;
}
