#include <zm/zm_nurbs_fitting.h>

#define VEC_SIZE 3 /* x,y,z */
#define ORDER 3
#define CP_NUM 10
#define DATA_NUM 30
#define SLICE_NUM 200
#define ITER 200
#define TOL 1.0e-6

/* out : ref_vec_array */
void _create_ref_vec_array(zVecArray *ref_vec_array)
{
  double sin, cos, t;
  int i;

  /* fitting data */
  const int data_num = DATA_NUM;
  zVecArrayAlloc( ref_vec_array, VEC_SIZE, data_num ); /* (x, y, z) x data_num */
  for( i=0; i<data_num; i++ ){
    t = (double)(i)/(double)(data_num - 1);
    zSinCos( 4.0*zPI*t, &sin, &cos );
    zVecArrayElem( ref_vec_array, i, 0 ) = cos; /* x */
    zVecArrayElem( ref_vec_array, i, 1 ) = sin; /* y */
    zVecArrayElem( ref_vec_array, i, 2 ) = 1.5 * t; /* z */
  }
}


int main(int argc, char *argv[]){
  zNURBS nurbs;
  zVecArray ref_vec_array;
  int iter_count;
  double start_knot, end_knot, ds, s;
  void *fit;

  _create_ref_vec_array( &ref_vec_array );
  start_knot = 0;
  end_knot = 30;

  zNURBSInit( &nurbs );
  zNURBSFitCreate( &ref_vec_array, ORDER, CP_NUM, start_knot, end_knot, &nurbs, &fit );
  zVecArrayFree( &ref_vec_array );
  zNURBSFitInitialize( fit, NULL, NULL );
  _zNURBSFitPrintFReport( fit, stdout );
  FILE *iteration_logfp;
  iteration_logfp = fopen( "fitting.csv", "w" );
  iter_count = zNURBSFitting( fit, ITER, TOL, iteration_logfp );
  fclose(iteration_logfp);
  printf("end : iteration count = %d\n\n", iter_count);
  _zNURBSFitPrintFReport( fit, stdout );

  FILE* fp_nurbs;
  fp_nurbs = fopen( "fitting_nurbs_curve.csv", "w" );
  zNURBSFitPrintF( fit )->header_label_of_fitted_curve( fp_nurbs );
  fprintf( fp_nurbs, "\n");
  ds = (end_knot - start_knot) / SLICE_NUM;
  for( s=0; s < end_knot + ds; ){
    zNURBSFitPrintF( fit )->fitted_curve( fp_nurbs, s );
    s += ds;
  }
  fclose( fp_nurbs );

  FILE* fp_cp_w;
  fp_cp_w = fopen( "fitting_nurbs_cp_and_weight.csv", "w" );
  zNURBSFitPrintF( fit )->control_point_and_weight( fp_cp_w );
  fclose( fp_cp_w );

  FILE* fp_data;
  fp_data = fopen( "fitting_target_and_fitted_nurbs_point.csv", "w" );
  zNURBSFitPrintF( fit )->target_and_fitted_point( fp_data );
  fclose( fp_data );


  zNURBSFitFree( &fit );
  zNURBSDestroy( &nurbs );

  return EXIT_SUCCESS;
}



