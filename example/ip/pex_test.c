#include <zm/zm_ip.h>

int main(int argc, char *argv[])
{
  zPexIP pc;
  double t;
  int dim, i;

  for( dim=1; dim<=3; dim++ ){
    zPexIPCreate( &pc, 3, dim );
    for( i=0; i<=dim; i++ )
      zPexIPSetCoeff( &pc, i, 1 );
    for( t=-zPexIPTerm(&pc); t<=zPexIPTerm(&pc); t+=0.01 )
      printf( "%f\n", zPexIPVal( &pc, t ) );
    zPexIPDestroy( &pc );
  }
  return 0;
}
