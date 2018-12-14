#include <zm/zm_opt.h>

#define TEST  7
#define SCALE 0
#define DIM   4

#if TEST == 1 /* Rosenbrock function (1,1,1,1) */
double eval(zVec var, void *dummy)
{
  int i;
  double x, y, val;

  for( val=0, i=0; i<zVecSizeNC(var); i+=2 ){
    x = zVecElem(var,i  );
    y = zVecElem(var,i+1);
    val += 100 * zSqr( y - x*x ) + zSqr( 1 - x );
  }
  return val;
}
zVec grad(zVec var, zVec g, void *dummy)
{
  int i;
  double x, y;

  for( i=0; i<zVecSizeNC(var); i+=2 ){
    x = zVecElem(var,i  );
    y = zVecElem(var,i+1);
    zVecSetElem( g, i  ,-400*x*(y-x*x)+2*(x-1) );
    zVecSetElem( g, i+1, 200*(y-x*x) );
  }
  return g;
}
zMat hess(zVec var, zMat h, void *dummy)
{
  int i;
  double x, y;

  for( i=0; i<zVecSizeNC(var); i+=2 ){
    x = zVecElem(var,i  );
    y = zVecElem(var,i+1);
    zMatSetElem( h, i  , i  ,-400*y+1200*x*x+2 );
    zMatSetElem( h, i  , i+1,-400*x );
    zMatSetElem( h, i+1, i  ,-400*x );
    zMatSetElem( h, i+1, i+1, 200 );
  }
  return h;
}
#elif TEST == 2 /* Beale function (3,0.5,3,0.5) */
#define BF_C1 1.5
#define BF_C2 2.25
#define BF_C3 2.625
double eval(zVec var, void *dummy)
{
  int i;
  double x, y, val;

  for( val=0, i=0; i<zVecSizeNC(var); i+=2 ){
    x = zVecElem(var,i  );
    y = zVecElem(var,i+1);
    val += zSqr( BF_C1 - x*(1-y) ) + zSqr( BF_C2 - x*(1-y*y) ) + zSqr( BF_C3 - x*(1-y*y*y) );
  }
  return val;
}
zVec grad(zVec var, zVec g, void *dummy)
{
  int i;
  double x, y;
  double c1, c2, c3;

  for( i=0; i<zVecSizeNC(var); i+=2 ){
    x = zVecElem(var,i  );
    y = zVecElem(var,i+1);
    c1 = 2 * ( BF_C1 - x*(1-y) );
    c2 = 2 * ( BF_C2 - x*(1-y*y) );
    c3 = 2 * ( BF_C3 - x*(1-y*y*y) );
    zVecSetElem( g, i  , -c1*(1-y)-c2*(1-y*y)-c3*(1-y*y*y) );
    zVecSetElem( g, i+1, x*( c1 + c2*2*y + c3 * 3*y*y ) );
  }
  return g;
}
#elif TEST == 3 && DIM == 4 /* Powell function (0,0,0,0) */
double eval(zVec var, void *dummy)
{
  double x1, x2, x3, x4;

  x1 = zVecElem(var,0);
  x2 = zVecElem(var,1);
  x3 = zVecElem(var,2);
  x4 = zVecElem(var,3);
  return zSqr(x1+10*x2) + 5*zSqr(x3-x4) + pow(x2-2*x3,4) + 10*pow(x1-x4,4);
}
zVec grad(zVec var, zVec g, void *dummy)
{
  double x1, x2, x3, x4, a, b, c, d;

  x1 = zVecElem(var,0);
  x2 = zVecElem(var,1);
  x3 = zVecElem(var,2);
  x4 = zVecElem(var,3);
  a = x1+10*x2;
  b = x3-x4;
  c = x2-2*x3;
  d = x1-x4;
  zVecSetElem( g, 0,  2*a+40*pow(d,3) );
  zVecSetElem( g, 1, 20*a+ 4*pow(c,3) );
  zVecSetElem( g, 2, 10*b- 8*pow(c,3) );
  zVecSetElem( g, 3,-10*b-40*pow(d,3) );
  return g;
}
#elif TEST == 4 /* ( 0, -0.01, -0.02, ... ) */
double eval(zVec var, void *dummy)
{
  register int i;
  double val;

  for( val=0, i=0; i<zVecSizeNC(var); i++ )
    val += zSqr( exp(zVecElem(var,i)+0.01*i)-1 );
  return 0.5 * val;
}
zVec grad(zVec var, zVec g, void *dummy)
{
  register int i;
  double e;

  for( i=0; i<zVecSizeNC(var); i++ ){
    e = exp( zVecElem(var,i)+0.01*i );
    zVecSetElem( g, i, (e-1)*e );
  }
  return g;
}
#elif TEST == 5 /* ( 1, 2, 3, ... ) */
double eval(zVec var, void *dummy)
{
  register int i;
  double val;

  for( val=0, i=0; i<zVecSizeNC(var); i++ )
    val += zSqr( log(zVecElem(var,i)-i) );
  return 0.5 * val;
}
zVec grad(zVec var, zVec g, void *dummy)
{
  register int i;

  for( i=0; i<zVecSizeNC(var); i++ )
    zVecSetElem( g, i, 2*log(zVecElem(var,i)-i)/fabs(zVecElem(var,i)-i) );
  return g;
}
#elif TEST == 6 /* (2,-1) */
double eval(zVec var, void *dummy)
{
  register int i;
  double x, y, val;

  for( val=0, i=0; i<zVecSizeNC(var); i+=2 ){
    x = zVecElem(var,i);
    y = zVecElem(var,i+1);
    val += pow(x-2,6) + zSqr(x-2)*y*y+pow(y+1,4);
  }
  return val;
}
zVec grad(zVec var, zVec g, void *dummy)
{
  register int i;
  double x, y, val;

  for( val=0, i=0; i<zVecSizeNC(var); i+=2 ){
    x = zVecElem(var,i);
    y = zVecElem(var,i+1);
    zVecSetElem( g, i,   6*pow(x-2,5)+2*(x-2)*y*y );
    zVecSetElem( g, i+1, 2*zSqr(x-2)*y+4*pow(y+1,3) );
  }
  return g;
}
#elif TEST == 7 /* from Mathematica tutorial (1.37638,1.67868) */
double eval(zVec var, void *dummy)
{
  register int i;
  double x, y, val;

  for( val=0, i=0; i<zVecSizeNC(var); i+=2 ){
    x = zVecElem(var,i);
    y = zVecElem(var,i+1);
    val += cos(x*x-3*y) + sin(x*x+y*y);
  }
  return val;
}
zVec grad(zVec var, zVec g, void *dummy)
{
  register int i;
  double x, y, val;

  for( val=0, i=0; i<zVecSizeNC(var); i+=2 ){
    x = zVecElem(var,i);
    y = zVecElem(var,i+1);
    zVecSetElem( g, i,   6*pow(x-2,5)+2*(x-2)*y*y );
    zVecSetElem( g, i+1, 2*zSqr(x-2)*y+4*pow(y+1,3) );
  }
  return g;
}
#else /* quadratic programming */
double eval(zVec var, void *dummy)
{
  register int i;
  double val;

  for( val=0, i=0; i<zVecSizeNC(var); i++ )
    val += (i+1)*zSqr( zVecElem(var,i)-i );
  return 0.5*val;
}
zVec grad(zVec var, zVec g, void *dummy)
{
  register int i;

  for( i=0; i<zVecSizeNC(var); i++ )
    zVecSetElem( g, i, (i+1)*(zVecElem(var,i)-i) );
  return g;
}
#endif

void test(zOptDM *opt, zVec var, const char* label)
{
  eprintf( "%s test: Hit enter key.", label );
  getchar();
  zVecSetAll( var, 1.0 );
  eprintf( " - initial vec = " );
  zVecFWrite( stderr, var );
  zOptDMSolve( opt, var, NULL, zTOL, 0, NULL );
  eprintf( " + result vec = " );
  zVecFWrite( stderr, var );
}

int main(int argc, char *argv[])
{
  zOptDM opt;
  zVec var;

  var = zVecAlloc( DIM );
  /*
  zOptDMCreate( &opt, DIM, SCALE, eval, grad, hess );
  zOptDMCreate( &opt, DIM, SCALE, eval, NULL, NULL );
  */
  zOptDMCreate( &opt, DIM, SCALE, eval, grad, NULL );

  zOptDMAssignSD( &opt, NULL );
  test( &opt, var, "SD" );
  zOptDMAssignLM( &opt, NULL );
  test( &opt, var, "LM (MT)" );
  zOptDMAssignLM( &opt, "DS" );
  test( &opt, var, "LM (DS)" );
  zOptDMAssignLM( &opt, "Noc" );
  test( &opt, var, "LM (Noc)" );
  zOptDMAssignVM( &opt, NULL, NULL );
  test( &opt, var, "VM (MT:BFGS)" );
  zOptDMAssignVM( &opt, NULL, "DFP" );
  test( &opt, var, "VM (MT:DFP)" );
  zOptDMAssignVM( &opt, "DS", "DFP" );
  test( &opt, var, "VM (DS:DFP)" );
  zOptDMAssignVM( &opt, "Noc", "DFP" );
  test( &opt, var, "VM (Noc:DFP)" );
  zOptDMAssignCG( &opt );
  test( &opt, var, "CG" );

  zOptDMDestroy( &opt );
  zVecFree( var );
  return 0;
}
