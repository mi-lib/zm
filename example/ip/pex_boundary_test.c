#include <zm/zm_ip.h>

#define T 1000

int main(int argc, char *argv[])
{
  zPexIP pc;
  register int i;
  double x;

  x = ( argc <= 1 ) ? 1 : atof(argv[1]);
  zPexIPCreateFromBoundary( &pc, T, 0, 0, 0, x, 0, 0, NULL );

  for( i=0; i<=T; i++ )
    printf( "%f\n", zPexIPVal( &pc, i ) );
  return 0;
}
