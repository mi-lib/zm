#include <zm/zm_nle.h>

#define DIM 3

zVec f1(zVec var, zVec r, void *dummy)
{
  int i;

  for( i=0; i<zVecSizeNC(var); i++ )
    zVecSetElem( r, i, sin(zVecElem(var,i))+zVecElem(var,i)-0.1*i*zPI );
  return r;
}
zMat jac1(zVec var, zMat j, void *dummy)
{
  int i;

  for( i=0; i<zVecSizeNC(var); i++ )
    zMatSetElem( j, i, i, cos(zVecElem(var,i))+1 );
  return j;
}

zVec f2(zVec var, zVec r, void *dummy)
{
  int i;

  for( i=0; i<zVecSizeNC(var); i++ )
    zVecSetElem( r, i, exp(zVecElem(var,i)-0.01*i)-1 );
  return r;
}
zMat jac2(zVec var, zMat j, void *dummy)
{
  int i;

  for( i=0; i<zVecSizeNC(var); i++ )
    zMatSetElem( j, i, i, exp(zVecElem(var,i)-0.01*i) );
  return j;
}

zVec f3(zVec var, zVec r, void *dummy)
{
  int i;

  for( i=0; i<zVecSizeNC(var); i++ )
    zVecSetElem( r, i, log(zVecElem(var,i)+i+1) );
  return r;
}
zMat jac3(zVec var, zMat j, void *dummy)
{
  int i;

  for( i=0; i<zVecSizeNC(var); i++ )
    zMatSetElem( j, i, i, 1.0/(zVecElem(var,i)+i+1) );
  return j;
}

zVec f4(zVec var, zVec r, void *dummy)
{
  int i;
  double x;

  for( i=0; i<zVecSizeNC(var); i++ ){
    x = zVecElem(var,i);
    zVecSetElem( r, i, (x-10.0*i)*(x*x+x+1)*(x*x-x+1) );
  }
  return r;
}
zMat jac4(zVec var, zMat j, void *dummy)
{
  int i;
  double x;

  for( i=0; i<zVecSizeNC(var); i++ ){
    x = zVecElem(var,i);
    zMatSetElem( j, i, i, (x*x+x+1)*(x*x-x+1)+2*x*(x-10.0*i)*(2*x*x+1) );
  }
  return j;
}

zVec f5(zVec var, zVec r, void *dummy)
{
  double x1, x2, x3;

  x1 = zVecElem(var,0);
  x2 = zVecElem(var,1);
  x3 = zVecElem(var,2);
  zVecSetElem( r, 0, x1*x1*x1-2*x2-2 );
  zVecSetElem( r, 1, x1*x1*x1-5*x3*x3-7 );
  zVecSetElem( r, 2, x2*x3*x3-1 );
  return r;
}
zMat jac5(zVec var, zMat j, void *dummy)
{
  double x1, x2, x3;

  x1 = zVecElem(var,0);
  x2 = zVecElem(var,1);
  x3 = zVecElem(var,2);
  zMatSetElem( j, 0, 0, 3*x1*x1 );
  zMatSetElem( j, 0, 1,-2 );
  zMatSetElem( j, 0, 2, 0 );
  zMatSetElem( j, 1, 0, 3*x1*x1 );
  zMatSetElem( j, 1, 1, 0 );
  zMatSetElem( j, 1, 2,-10*x3 );
  zMatSetElem( j, 2, 0, 0 );
  zMatSetElem( j, 2, 1, x3*x3 );
  zMatSetElem( j, 2, 2,2*x2*x3 );
  return j;
}

zVec f6(zVec var, zVec r, void *dummy)
{
  int i;

  for( i=0; i<zVecSizeNC(var); i++ )
    zVecSetElem( r, i, zVecElem(var,i)-3*i );
  return r;
}
zMat jac6(zVec var, zMat j, void *dummy)
{
  int i;

  for( i=0; i<zVecSizeNC(var); i++ )
    zMatSetElem( j, i, i, 1 );
  return j;
}

zVec (* f[])(zVec,zVec,void*) = {
  f1, f2, f3, f4, f5, f6,
};
zMat (* jac[])(zVec,zMat,void*) = {
  jac1, jac2, jac3, jac4, jac5, jac6,
};

void check(const char *label, zNLE *nle, int id, zVec var, zVec r, void *util)
{
  int ret;

  zVecSetAll( var, 1.0 );
  ret = zNLESolve( nle, var, util, zTOL, zTOL, 0, NULL );
  printf( "(%s) returned %d\n", label, ret );
  printf( "   answer: " );
  zVecTouchup( var, zTOL );
  zVecPrint( var );
  f[id]( var, r, util );
  zVecTouchup( r, zTOL );
  printf( " residual: " );
  zVecPrint( r );
}

int main(int argc, char *argv[])
{
  zNLE nle;
  zVec var, r;
  void *util = (void *)"hoge";
  int i = 0;

  var = zVecAlloc( DIM );
  r   = zVecAlloc( DIM );

  if( argc > 1 && atoi(argv[1]) < 6 )
    i = atoi( argv[1] );
  zNLECreate( &nle, DIM, DIM, f[i], argc > 2 ? jac[i] : NULL );

  zNLEAssignLM( &nle, "grad", "MT", 0 );
  check( "LM (grad, MT)", &nle, i, var, r, util );
  zNLEAssignLM( &nle, "dec", "MT", 0 );
  check( "LM (dec, MT)", &nle, i, var, r, util );
  zNLEAssignLM( &nle, "grad", "DS", 0 );
  check( "LM (grad, DS)", &nle, i, var, r, util );
  zNLEAssignNR( &nle );
  check( "NR", &nle, i, var, r, util );
  zNLEAssignBroyden( &nle );
  check( "Broyden", &nle, i, var, r, util );

  zNLEDestroy( &nle );
  zVecFreeAtOnce( 2, var, r );
  return 0;
}
