/* interpolation-in-order test */
#include <zm/zm_ip.h>

int main(void)
{
  double p, dt;
  double t;
  zIPIO ip;

  zIPIOCreate( &ip, 0.0 );
  while( 1 ){
    eprintf( "enter dx-value to the next step and y-value.\n" );
    eprintf( "(to terminate the operation, enter 0 for dx)\n" );
    scanf( "%lf %lf", &dt, &p );
    if( dt == 0 ) break;
    zIPIOSetNextValue( &ip, p, dt );
    zIPIOUpdate( &ip );
    for( t=0; t+zTOL<zIPIODT(&ip); t+=0.1 ){
      eprintf( "%f\n", t );
      printf( "%f %f\n", zIPIOValue( &ip, t ), zIPIOLinValue( &ip, t ) );
    }
  }
  return 0;  
}
