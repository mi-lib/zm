#include <zm/zm.h>

zVec f1(zVec var, zVec r, void *dummy)
{
  register int i;

  for( i=0; i<zVecSizeNC(var); i++ )
    zVecSetElem( r, i, sin(zVecElem(var,i))+zVecElem(var,i)-0.1*i*zPI );
  return r;
}
zMat jac1(zVec var, zMat j, void *dummy)
{
  register int i;

  for( i=0; i<zVecSizeNC(var); i++ )
    zMatSetElem( j, i, i, cos(zVecElem(var,i))+1 );
  return j;
}


int NLESolveLM(zVec var, int nc, zVec (* f)(zVec,zVec,void*), zMat (* j)(zVec,zMat,void*), double tol, int iternum, void *util)
{
  zVec x, r;
  zMat jac, _j;
  double e;
  register int i;

  x = zVecCopy( var );
  r = zVecAlloc( nc );
  jac = zMatAlloc( nc, zVecSizeNC(var) );
  _j = zMatAllocSqr( nc );
  ZITERINIT( iternum );
  for( i=0; i<iternum; i++ ){
    f( var, r, util );
    e = zVecSqrNorm( r );
    if( zIsTol( sqrt(e), tol ) ) goto TERMINATE;


  }
  ZITERWARN( i );
 TERMINATE:
  zVecFree( x );
  zVecFree( r );
  zMatFree( jac );
  return i;
}



#define DIM 3

int main(int argc, char *argv[])
{
  zVec var, r;

  var = zVecAlloc( DIM );
  r   = zVecAlloc( DIM );

  zVecSetAll( var, 1.0 );

  NLESolveLM( var, f1, jac1, zTOL, 0, util );
  printf( " answer  : " );
  zVecWrite( var );
  f1( var, r, NULL );
  printf( " residual: " );
  zVecWrite( r );

  zNLEDestroy( &nle );
  zVecFree( var );
  return 0;
}
