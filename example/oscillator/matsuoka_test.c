#include <zm/zm.h>

typedef struct{
  double y;     /* output */
  double x, dx; /* membrane potential */
  double f, df; /* degree of fatigue */
} zOscMatUnit;

void zOscMatUnitInit(zOscMatUnit *u, double x0)
{
  u->x  = x0;
  u->f  = 0;
  u->dx = u->df = 0;
  u->y  = 0;
}

void zOscMatUnitUpdate(zOscMatUnit *u, double s, double yo, double tr, double ta, double b, double dt)
{
  u->dx = ( s - yo - u->x - b * u->f ) / tr;
  u->df = (          u->y -     u->f ) / ta;
  u->x += u->dx * dt;
  u->f += u->df * dt;
}

double zOscMatUnitOutput(zOscMatUnit *u)
{
  return ( u->y = u->x > 0 ? u->x : 0 );
}

typedef struct{
  double tr; /* rise time constant */
  double ta; /* adaptation time constant */
  double b;  /* steady-state firing rate */
  double a;  /* coupling strength */
  zOscMatUnit u[2];
} zOscMat;

void zOscMatCreate(zOscMat *osc, double tr, double ta, double b, double a)
{
  osc->tr = tr;
  osc->ta = ta;
  osc->b  = b;
  osc->a  = a;
}

void zOscMatInit(zOscMat *osc, double x0, double x1)
{
  zOscMatUnitInit( &osc->u[0], x0 );
  zOscMatUnitInit( &osc->u[1], x1 );
}

double zOscMatOutput(zOscMat *osc)
{
  return osc->u[0].y - osc->u[1].y;
}

double zOscMatUpdate(zOscMat *osc, double u, double dt)
{
  zOscMatUnitUpdate( &osc->u[0], u, osc->u[1].y * osc->a, osc->tr, osc->ta, osc->b, dt );
  zOscMatUnitUpdate( &osc->u[1], u, osc->u[0].y * osc->a, osc->tr, osc->ta, osc->b, dt );
  zOscMatUnitOutput( &osc->u[0] );
  zOscMatUnitOutput( &osc->u[1] );
  return zOscMatOutput( osc );
}


#define DT 0.005
#define T  100.0

int main(int argc, char *argv[])
{
  zOscMat osc;
  double u;
  int i, step;

  u = argc > 1 ? atof( argv[1] ) : 1.0;
  step = T / DT;
  zOscMatCreate( &osc, 1.0, 6.0, 2.5, 1.5 );
  zOscMatInit( &osc, 0.1, 0.0 );
  for( i=0; i<=step; i++ ){
    zOscMatUpdate( &osc, u, DT );
    printf( "%f %f %f %f\n", i*DT, osc.u[0].y, osc.u[1].y, zOscMatOutput(&osc) );
  }
  return 0;
}
