#include <zm/zm_intg.h>

#define TEST 1

double f_org(double t)
{
#if TEST == 1
  return sin(t) + t + 1;
#else
  return 0.5 * t * ( t - 1 ) + 1;
#endif
}

double f_dif(double t)
{
#if TEST == 1
  return cos(t) + 1;
#else
  return t - 0.5;
#endif
}

void output(double t, double s_true, zIntg intg[], int n)
{
  int i;

  printf( "%f %.16f", t, s_true );
  for( i=0; i<n; i++ )
    printf( " %.16f %.16f", intg[i].s, s_true-intg[i].s );
  printf( "\n" );
}

#define STEP 1000

#define INTG_N 3

int main(void)
{
  register int i;
  zIntg intg[INTG_N];
  double s, x, t, dt;

  s = f_org(0);
  x = f_dif(0);
  for( i=0; i<INTG_N; i++ )
    zIntgInit( &intg[i], s, x );

  for( t=0; t<10; ){
    dt = zRandF( 0.01, 0.05 );
    t += dt;
    x = f_dif(t);
    zIntgRect( &intg[0], x, dt );
    zIntgTR( &intg[1], x, dt );
    zIntgQuad( &intg[2], x, dt );
    output( t, f_org(t), intg, INTG_N );
  }
  return 0;
}
