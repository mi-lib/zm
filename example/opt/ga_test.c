#include <zm/zm_opt.h>

double fitness(zOptGAChromosome *c, void *dummy)
{
  double x, y;

  x = zVecElem(c->gene,0);
  y = zVecElem(c->gene,1);
  return    exp(-0.3*(zSqr(x-5)+zSqr(y-7)))
         +5*exp(-2.0*(zSqr(x+4)+zSqr(y+5)))
         +2*exp(-0.5*(zSqr(x-5)+zSqr(y+8)))
         +3*exp(-0.3*(zSqr(x+6)+zSqr(y-6)));
}

#define GENERATION 100

int main(int argc, char *argv[])
{
  zOptGA ga;
  zVec min, max, ans;
#if 1
  register int i;
  char filename[BUFSIZ];
  FILE *fp;
#endif

  zRandInit();
  min = zVecCreateList( 2,-10.0,-10.0 );
  max = zVecCreateList( 2, 10.0, 10.0 );
  ans = zVecAlloc( 2 );
  /* creation of population */
  zOptGACreate( &ga, 2, 1000, min, max, fitness, 0.3, 0.05, NULL );
#if 1
    sprintf( filename, "log/000" );
    fp = fopen( filename, "w" );
    zOptGAFWrite( fp, &ga );
    fclose( fp );
#endif
  zVecFree( min );
  zVecFree( max );

#if 1
  for( i=1; i<=GENERATION; i++ ){
    zOptGAReproduce( &ga, NULL );
    zVecWrite( ga.individual[0].gene );
#if 1
    sprintf( filename, "log/%03d", i );
    fp = fopen( filename, "w" );
    zOptGAFWrite( fp, &ga );
    fclose( fp );
#endif
  }
#else
  zOptGASolve( &ga, ans, NULL, GENERATION );
  zVecWrite( ga.individual[0].gene );
#endif

  zOptGADestroy( &ga );
  return 0;
}
