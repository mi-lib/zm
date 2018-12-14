#include <zm/zm_mat.h>

int main(void)
{
  zMat m;

  zRandInit();
  m = zMatAlloc( 5, 6 );
  zMatRandUniform( m, -10, 10 );
  zMatWrite( m );

  zMatColReg( m, 3 );
  zMatWrite( m );

  zMatFree( m );
  return 0;
}
