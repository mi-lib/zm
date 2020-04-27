#include <zm/zm_opt.h>
#include <sys/stat.h>

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

#define GENERATION 100
int main(int argc, char *argv[])
{
  zOptGA ga;
  zVec min, max;
  register int i;
  char filename[BUFSIZ];
  FILE *fp;

  create_min_max( &min, &max );
  /* create population */
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
