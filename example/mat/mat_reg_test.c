#include <zm/zm_mat.h>

int main(void)
{
  zMat m;

  zRandInit();
  m = zMatAlloc( 5, 6 );
  zMatRandUniform( m, -10, 10 );
  zMatPrint( m );

  zMatColReg( m, 3 );
  zMatPrint( m );

  zMatFree( m );
  return 0;
}
