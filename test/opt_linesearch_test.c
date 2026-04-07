#include <zm/zm_opt.h>

double eval0(double x, void *dummy){ /* x=1 */
  return (x-1)*(x-1)*(x*x+6*x+15);
}
double eval1(double x, void *dummy){ /* x=2 */
  return (x-2)*(x-2)-0.2*cos(zPIx2*(x-2));
}
double eval2(double x, void *dummy){ /* x=3 */
  return fabs( (x-3)*(x-12)*(x+12) );
}
double eval3(double x, void *dummy){ /* x=4 */
  return (x*x+x+3)*(x-4)*(x-4);
}
double eval4(double x, void *dummy){ /* x=2.5 */
  return (x*x+0.3)*(x-1)*(x-4)*((x-5)*(x-5)+0.3);
}
void assert_opt_linesearch(void)
{
  double (* eval_fp[])(double,void*) = {
    eval0,
    eval1,
    eval2,
    eval3,
    eval4,
    NULL,
  };
  double answer[] = {
    1.0,
    2.0,
    3.0,
    4.0,
    2.5,
  };
  double (** fp)(double,void*);
  double x;
  int i;
  bool result1, result2, result3;
  const double tol = 1.0e-6;

  result1 = result2 = result3 = true;
  for( i=0, fp=eval_fp; *fp; fp++, i++ ){
    x = zOptLineGoldenSection( *fp, -10, 10, NULL, 0 );
    if( !zEqual( x, answer[i], tol ) ){
      eprintf( "golden section method:  x* = %.10f ( evaluation = %.10f )\n", x, (*fp)( x, NULL ) );
      result1 = false;
    }
    x = zOptLineTrisection( *fp, -10, 10, NULL, 0 );
    if( !zEqual( x, answer[i], tol ) ){
      eprintf( "trisection method:      x* = %.10f ( evaluation = %.10f )\n", x, (*fp)( x, NULL ) );
      result2 = false;
    }
    x = zOptLineBrent( *fp, -10, 10, NULL, 0 );
    if( !zEqual( x, answer[i], tol ) ){
      eprintf( "Brent's method :        x* = %.10f ( evaluation = %.10f )\n", x, (*fp)( x, NULL ) );
      result3 = false;
    }
  }
  zAssert( zOptLineGoldenSection, result1 );
  zAssert( zOptLineTrisection, result2 );
  zAssert( zOptLineBrent, result3 );
}

int main(void)
{
  assert_opt_linesearch();
  return 0;
}
