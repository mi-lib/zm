#include <zm/zm.h>

#define G  9.8

#define DT 0.002
#define STEP 3000

#define T1 30
#define T2  0

#define M1 1.0
#define M2 1.0
#define L1 1.0
#define L2 1.0

/* equation of motion */
zVec dualpend_eqm(double t, zVec x, void *dummy, zVec dx)
{
  double a1, a2, a3, v1, v2, det;
  double t1, t2, t1_dot, t2_dot;
  double u;

  u = *(double *)dummy;
  t1 = zVecElem( x, 0 );
  t2 = zVecElem( x, 1 );
  t1_dot = zVecElem( x, 2 );
  t2_dot = zVecElem( x, 3 );

  a1 = M1*L1*L1 + M2*(L1*L1+2*L1*L2*cos(t2)+L2*L2);
  a2 = M2*(L1*cos(t2)+L2)*L2;
  a3 = M2*L2*L2;
  det = a1*a3 - a2*a2;
  v1 = M2*L1*L2*(2*t1_dot+t2_dot)*t2_dot*sin(t2)
    - (M1+M2)*G*L1*sin(t1) - M2*G*L2*sin(t1+t2) + u;
  v2 =-2*M2*L1*L2*t1_dot*(t1_dot+t2_dot)*sin(t2) - M2*G*L2*sin(t1+t2);

  zVecSetElem( dx, 0, t1_dot );
  zVecSetElem( dx, 1, t2_dot );
  zVecSetElem( dx, 2, ( a3*v1 - a2*v2 ) / det );
  zVecSetElem( dx, 3, (-a2*v1 + a1*v2 ) / det );
  return dx;
}

void output(int i, zVec x)
{
  FILE *fp;
  static char buf[BUFSIZ], nbuf[5];
  double x1, y1, x2, y2;

  sprintf( buf, "dat%s", itoa_zerofill(i,4,nbuf) );
  fp = fopen( buf, "w" );
  eprintf( "%f %f\n", zVecElem(x,0), zVecElem(x,1) );
  x1 = L1 * sin( zVecElem(x,0) );
  y1 =-L1 * cos( zVecElem(x,0) );
  x2 = x1 + L2 * sin( zVecElem(x,0)+zVecElem(x,1) );
  y2 = y1 - L2 * cos( zVecElem(x,0)+zVecElem(x,1) );
  fprintf( fp, "0 0\n" );
  fprintf( fp, "%f %f\n", x1, y1 );
  fprintf( fp, "%f %f\n", x2, y2 );
  fclose( fp );
}

int main(int argc, char *argv[])
{
  register int i;
  zVec x;
  double u = 0;
  zODE ode;

  zODEAssign( &ode, RK4, NULL, NULL );
  zODEInit( &ode, 4, 0, dualpend_eqm );
  x = zVecAlloc( 4 );
  zVecSetElem( x, 0, zDeg2Rad( argc>1 ? atof(argv[1]) : T1 ) );
  zVecSetElem( x, 1, zDeg2Rad( argc>2 ? atof(argv[2]) : T2 ) );
  for( i=0; i<STEP; i++ ){
    output( i, x );
    zODEUpdate( &ode, 0, x, DT, &u );
  }
  output( STEP, x );
  zVecFree( x );
  zODEDestroy( &ode );
  return 0;
}
