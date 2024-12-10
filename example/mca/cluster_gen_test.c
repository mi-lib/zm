#include <zm/zm_mca.h>
#include <zm/zm_rand.h>

double gen_r_uniform(double rmax){ return zRandF( 0, rmax ); }
double gen_r_normald(double rmax){ return zRandND( NULL, 0, 0.5*rmax ); }

bool gen_vec(double (* gen_func)(double), int dim, int np, int nc, double rate_outlier, double min, double max)
{
  zVec vc, v;
  double rmax, r;
  int i, j, k;
  bool ret = true;

  vc = zVecAlloc( dim );
  v  = zVecAlloc( dim );
  if( !vc || !v ){
    ret = false;
    goto TERMINATE;
  }
  rmax = fabs( max - min ) * sqrt(dim) / (2*nc);
  for( i=0; i<nc; i++ ){
    /* center of a cluster */
    for( k=0; k<dim; k++ )
      zVecSetElem( vc, k, zRandF( min, max ) );
    /* sample */
    for( j=0; j<np; j++ ){
      r = gen_func( rmax );
      for( k=0; k<dim; k++ )
        zVecSetElem( v, k, zRandF(min-max,max-min) );
      zVecNormalizeDRC( v );
      zVecMulDRC( v, r );
      zVecAddDRC( v, vc );
      zVecValuePrint( v );
    }
  }
  /* outliers */
  np *= rate_outlier;
  for( j=0; j<np; j++ ){
    for( k=0; k<dim; k++ )
      zVecSetElem( v, k, zRandF(min,max) );
    zVecValuePrint( v );
  }
 TERMINATE:
  zVecFree( vc );
  zVecFree( v );
  return ret;
}

#define DIM           2
#define NUM_POINTS 1000
#define NUM_CLUSTER   6
#define RATE_OUTLIER  0.01
#define MIN           0
#define MAX          10

int main(int argc, char *argv[])
{
  zRandInit();
  if( !gen_vec( argc > 1 && argv[1][0] == 'u' ? gen_r_uniform : gen_r_normald,
                argc > 2 ? atoi( argv[2] ) : DIM,
                argc > 3 ? atoi( argv[3] ) : NUM_POINTS,
                argc > 4 ? atoi( argv[4] ) : NUM_CLUSTER,
                argc > 5 ? atof( argv[5] ) : RATE_OUTLIER,
                argc > 6 ? atof( argv[6] ) : MIN,
                argc > 7 ? atof( argv[7] ) : MAX ) ) return EXIT_FAILURE;
  return EXIT_SUCCESS;
}
