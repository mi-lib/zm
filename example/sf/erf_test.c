#include <zm/zm_sf.h>

double erf_ap1(double x)
{
  double a, x2, val;

  x2 = x*x;
  a = ( 8 * ( zPI - 3 ) / ( 3 * zPI * ( zPI - 4 ) ) ) * x2;
  val = sqrt( 1 - exp( -x2 * ( 4/zPI + a ) / ( 1 + a ) ) );
  return x >= 0 ? val : -val;
}

/* http://invar6.blog.fc2.com/blog-entry-16.html */
double erf_ap2(double x)
{
  long double x2, val;

  x2 = x*x;
  val = sqrt( 1 - exp( x2 *  -1.2732395447351626861510701069801 )
          * ( 1 + x2 * x2 * ( 3.8256894018072952549677698975590e-2
                     + x2 * (-3.3924118631479317799557719571888e-3
                     + x2 * ( 6.8005689892920536537036277150071e-4
                     + x2 * (-7.0227051009098722650256537001220e-5
                     + x2 * ( 7.6148050059630352704610740846049e-6
                     + x2 * (-6.6700609990522378462586732796781e-7
                     + x2 * ( 5.3870107438215609185130430012239e-8
                     + x2 * (-3.8693709096400193758628383851182e-9
                     + x2 *   2.5441273130247078605512789263896e-10 ) ) ) ) ) ) ) ) ) );
  return x >= 0 ? val : -val;
}

/* Williamson=Yamashita formula */
double erf_wy(double x)
{
  long double x2, val;

  x2 = x*x;
  val = sqrt( 1 - exp( x2 * -1.2732395447351626861510701069801 )
          * ( 1 + x2 * x2 * ( 0.1101999999998335
                   / ( x2 +   7.2085122317063002 )
                            + 0.0219999999999997 ) ) );
  return x >= 0 ? val : -val;
}


#define W 5.0
#define DIV 1000

int main(int argc, char *argv[])
{
  register int i;
  double x, f, f0, f1, f2, f3;

  for( i=0; i<=DIV; i++ ){
    x = W * ( (double)i/DIV - 0.5 );
    f  = erf(x);
    f0 = zErf(x);
    f1 = erf_ap1(x);
    f2 = erf_ap2(x);
    f3 = erf_wy(x);
    printf( "%.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g\n", x, f, f0, f1, f2, f3, f-f0, f-f1, f-f2, f-f3 );
  }
  return 0;
}
