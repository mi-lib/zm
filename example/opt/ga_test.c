#include <zm/zm_opt.h>
#include <sys/stat.h>

#if 0
double testfunc(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return   -exp(-0.3*(zSqr(x1-5)+zSqr(x2-7)))
         -5*exp(-2.0*(zSqr(x1+4)+zSqr(x2+5)))
         -2*exp(-0.5*(zSqr(x1-5)+zSqr(x2+8)))
         -3*exp(-0.3*(zSqr(x1+6)+zSqr(x2-6)));
}

void create_min_max(zVec *min, zVec *max)
{
  *min = zVecCreateList( 2,-10.0,-10.0 );
  *max = zVecCreateList( 2, 10.0, 10.0 );
}
#else
/* Goldman & Price function */
double testfunc(zVec x, void *dummy)
{
  double x1, x2;

  x1 = zVecElem(x,0);
  x2 = zVecElem(x,1);
  return ( 1 + zSqr(x1+x2+1) * (19-14*x1+3*x1*x1-14*x2+6*x1*x2+3*x2*x2) )
       * ( 30 + zSqr(2*x1-3*x2) * (18-32*x1+12*x1*x1+48*x2-36*x1*x2+27*x2*x2) );
}

void create_min_max(zVec *min, zVec *max)
{
  *min = zVecCreateList( 2,-2.0,-2.0 );
  *max = zVecCreateList( 2, 2.0, 2.0 );
}
#endif


#if 0
#define GENERATION 100
int main(int argc, char *argv[])
{
  zOptGA ga;
  zVec min, max;
  register int i;
  char filename[BUFSIZ];
  FILE *fp;

  create_min_max( &min, &max );
  /* creation of population */
  mkdir( "log", 755 );
  zOptGACreate( &ga, testfunc, NULL, min, max, 1000, 0.3, 0.05 );
    sprintf( filename, "log/000" );
    fp = fopen( filename, "w" );
    zOptGAFPrint( fp, &ga );
    fclose( fp );
  zVecFree( min );
  zVecFree( max );

  for( i=1; i<=GENERATION; i++ ){
    zOptGAReproduce( &ga, NULL );
    zVecPrint( ga.individual[0].gene );
    sprintf( filename, "log/%03d", i );
    fp = fopen( filename, "w" );
    zOptGAFPrint( fp, &ga );
    fclose( fp );
  }
  zVecPrint( ga.individual[0].gene );

  zOptGADestroy( &ga );
  return 0;
}
#else
int main(int argc, char *argv[])
{
  zVec min, max, ans;
  double val;

  create_min_max( &min, &max );
  ans = zVecAlloc( 2 );
  zOptSolveGADefault( testfunc, NULL, min, max, 0, zTOL, ans, &val );
  zVecPrint( ans );
  printf( "%.10g\n", val );
  zVecFree( min );
  zVecFree( max );
  return 0;
}
#endif
