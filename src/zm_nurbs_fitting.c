/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 2025, Daishi Kaneta (kdaic), Tomomichi Sugihara (Zhidao),
 *
 * zm_nurbs_fitting - NURBS curve fitting with Gauss Newton algorighm
 *
 */

#include <zm/zm_nurbs_fitting.h>
#include <zm/zm_mat.h>
#include <zm/zm_le.h>

#define HUGE_VALUE 1.0e9
#define TINY_VALUE 1.0e-9

/* @param[in] nurbs      zNURBS of which the number of control point array is n */
/* @param[in] knot_m     m-th knot */
/* @param[out] j_knot_m  minimum element of Jacobian w.r.t. knot. the vector size is NJ x m  */
zVec _zNURBSFitKnotJacobianElementVec(const zNURBS *nurbs, const double knot_m, zVec j_knot_m)
{
  return zNURBSVecDiff( nurbs, knot_m, 1, j_knot_m );
}

/* @param[in] nurbs     zNURBS of which the number of control point array of each nurbs is n */
/* @param[in] fit_knots knot vector. the vec size is m */
/* @param[out] J_knot   Jacobian w.r.t. knot vector ( NJ m x m ) */
zMat _zNURBSFitKnotJacobian(const zNURBS *nurbs, const zVec fit_knots, zMat J_knot)
{
  int p_size, m, i;
  zVec j_knot_m;

  zMatZero( J_knot );
  p_size = zVecSize( zNURBSCP(nurbs, 0) );
  if( p_size < 1 ){
    ZRUNERROR( ZM_ERR_MAT_SIZEMISMATCH );
    return NULL;
  }
  j_knot_m = zVecAlloc( p_size );
  for( m=0; m < zVecSize(fit_knots); m++ ){
    _zNURBSFitKnotJacobianElementVec( nurbs, zVecElemNC(fit_knots, m), j_knot_m );
    for( i=0; i < p_size; i++ )
      zMatSetElemNC( J_knot, m*p_size+i, m, zVecElemNC(j_knot_m, i) );
  }
  zVecFree( j_knot_m );

  return J_knot;
}

/* @param[in] nurbs    zNURBS of which the number of control point array is n */
/* @param[in] knot_m   m-th knot */
/* @param[in] idx_cp_n the index of n-th control point v0..vn[0..NJ]. this is the X of the [0..X..NJ] */
/* @param[out] j_cp_n  minimum element of Jacobian w.r.t. idx_cp_n of the n-th control point vn */
double _zNURBSFitControlPointJacobianElement(const zNURBS *nurbs, const double knot_m, const int idx_cp_n, double *j_cp_n)
{
  int seg, i;
  double numerater, den; /* den : denominator */

  seg = zBSplineParamSeg( &nurbs->param, knot_m );
  numerater = zNURBSWeight( nurbs, idx_cp_n ) * zBSplineParamBasis( &nurbs->param, knot_m, idx_cp_n, nurbs->param.order, seg );
  if( zIsNan(numerater) || zIsInf(numerater) )
    numerater = numerater/fabs(numerater) * HUGE_VALUE;
  for( den=0, i=seg-nurbs->param.order; i<=seg; i++ )
    den += zNURBSWeight( nurbs, i ) * zBSplineParamBasis( &nurbs->param, knot_m, i, nurbs->param.order,seg );
  *j_cp_n = ( zIsTiny(den) || zIsNan(den) || zIsInf(den) ) ? numerater : numerater/den;
  return *j_cp_n;
}

/* @param[in] nurbs     zNURBS of which the number of control point array of each nurbs is n */
/* @param[in] fit_knots knot vector. the vec size is m */
/* @param[out] J_cp     Jacobian w.r.t. control point ( N_J n x N_J m ) */
zMat _zNURBSFitControlPointJacobian(const zNURBS *nurbs, const zVec fit_knots, zMat J_cp)
{
  int p_size, m, n, i;
  double j_cp_n;

  zMatZero( J_cp );
  p_size = zVecSize( zNURBSCP(nurbs, 0) );
  for( m=0; m < zVecSize(fit_knots); m++ ){
    for( i=0; i < p_size; i++ ){
      for( n=0; n<zNURBSCPNum( nurbs ); n++ ){
        _zNURBSFitControlPointJacobianElement( nurbs, zVecElemNC(fit_knots, m), n, &j_cp_n );
        zMatSetElemNC( J_cp, m*p_size + i, n*p_size + i, j_cp_n);
      }
    }
  }
  return J_cp;
}

/* @param[in] nurbs  zNURBS of which the number of control point array is n */
/* @param[in] knot_m m-th knot */
/* @param[in] idx_w  the index of the weight w0..wn. this is the X of the w0..wX..wn */
/* @param[out] j_wn  minimum vector element of Jacobian w.r.t. idx_wn of the n-th weight */
zVec _zNURBSFitWeightJacobianElementVec(const zNURBS *nurbs, const double knot_m, const int idx_w, zVec j_wn)
{
  int seg, i;
  double b, den; /* den : denominator */

  seg = zBSplineParamSeg( &nurbs->param, knot_m );
  zVecZero( j_wn );
  zNURBSVec( nurbs, knot_m, j_wn );
  zVecMulDRC( j_wn, -1.0 );
  zVecAddDRC( j_wn, zNURBSCP( nurbs, idx_w ) );

  /* weight multipling means w -> exp(w), exp(w) appears after derivative. */
  b = zNURBSWeight(nurbs,idx_w) * zBSplineParamBasis( &nurbs->param, knot_m, idx_w, nurbs->param.order, seg );

  zVecMulDRC( j_wn, b );

  if( zVecIsNan( j_wn ) )
    zVecSetAll( j_wn, HUGE_VALUE );

  for( den=0, i=seg-nurbs->param.order; i<=seg; i++ )
    den += zNURBSWeight(nurbs,i) * zBSplineParamBasis( &nurbs->param, knot_m, i, nurbs->param.order, seg );
  return ( zIsTiny(den) || zIsNan(den) || zIsInf(den) ) ? j_wn : zVecDivDRC( j_wn, den );
}

/* @param[in] nurbs     zNURBS of which the number of control point array of each nurbs is n */
/* @param[in] fit_knots knot vector. the vec size is m */
/* @param[out] J_w      Jacobian w.r.t. weight ( N_J m x n ) */
zMat _zNURBSFitWeightJacobian(const zNURBS *nurbs, const zVec fit_knots, zMat J_w)
{
  int p_size, m, n, i;
  zVec j_wn;
  zMatZero( J_w );
  p_size = zVecSize( zNURBSCP(nurbs, 0) );
  j_wn = zVecAlloc( p_size );
  for( m=0; m < zVecSize(fit_knots); m++ ){
    for( n=0; n<zNURBSCPNum( nurbs ); n++ ){
      _zNURBSFitWeightJacobianElementVec( nurbs, zVecElemNC(fit_knots, m), n,  j_wn );
      for( i=0; i < p_size; i++ )
        zMatSetElemNC( J_w, m*p_size+i, n, zVecElemNC(j_wn, i) );
    }
  }
  zVecFree( j_wn );
  return J_w;
}

/* @param[in] nurbs          zNURBS of which the number of control point array of each nurbs is n */
/* @param[in] fit_knots      knot vecotr. */
/* @param[in] is_fixed_knots the flag whether knot is fixed or not. If fixed, fit_knots will be not updated */
/* @param[out] J_knot        Jacobian of knot */
/* @param[out] J_cp          Jacobian of control point */
/* @param[out] J_w           Jacobian of weight */
bool _zNURBSFitCalcSubJacobian(const zNURBS *nurbs, const zVec fit_knots, const bool is_fixed_knots, zMat J_knot, zMat J_cp, zMat J_w)
{
  if( !is_fixed_knots )
    _zNURBSFitKnotJacobian( nurbs, fit_knots, J_knot );
  _zNURBSFitControlPointJacobian( nurbs, fit_knots, J_cp );
  _zNURBSFitWeightJacobian( nurbs, fit_knots, J_w );
  return true;
}

/* @param[in] nurbs      zNURBS of which the number of control point array of each nurbs is n */
/* @param[in] fint_knots knot vector */
/* @param[out] J_knot    Jacobian of knot */
/* @param[out] J_cp      Jacobian of control point */
/* @param[out] J_w       Jacobian of weight */
bool _zNURBSCreateSubJacobian(const zNURBS *nurbs, const zVec fit_knots, zMat *J_knot, zMat *J_cp, zMat *J_w)
{
  int data_num, cp_num, p_size;
  data_num = zVecSize( fit_knots );
  cp_num = zNURBSCPNum( nurbs );
  p_size = zVecSize( zNURBSCP(nurbs,0) );
  *J_knot = zMatAlloc( p_size * data_num, data_num  );
  *J_cp = zMatAlloc( p_size * data_num, p_size * cp_num );
  *J_w = zMatAlloc( p_size * data_num, cp_num );

  return true;
}

/* linear eqaution parameter to solve by least squares method */
/* See the definition of zLEWorkspace in zm_le_gen.h */
typedef struct{
  zMat a, ident;
  zVec b, wn, we, residual, err_vec;
  zLEWorkspace workspace;
  double eval_eq, eval_dist;
  int iter_count;
} leParam;

/* */
typedef struct{
  int data_num; /* the number of data = fit_knot size (= m). */
  zVec dx; /* difference of main fitting parameters(=x) to solve optimal non-linear equation. */
  zVec x;  /* solved main fitting parameters. x = [ fit_knots t, control point v, weight w ] ((NJ m + NJ n + n) x 1). */
  zVec _cp_v; /* temporary parameter for update a control point (NJ x 1). */
  zNURBS *nurbs; /* fitting (updating) NURBS main data which has fitting (updating) control point (NJ n) & weight (n). */
  zVec fit_knots; /* fitting (updating) knots for fit_vec within NURBS definition. (m x 1). */
  bool is_fixed_knots;
  zVec ref_vec; /* input target referece points data to fit. (NJ m x 1). */
  zVec fit_vec; /* fitting (updating) points on NURBS. (NJ m x 1). */
  zMat J_knot; /* jacobian of knot */
  zMat J_cp; /* jacobian of control point */
  zMat J_w; /* jacobian of weight */
  zMat J; /* jacobian of fitting knot, control point, weight */
  leParam le; /* paramters of linear equation solver */
  zNURBSFitPrintFClass *print; /* instance of print class */
} zNURBSFitData;


/* Calculates the Jacobian for minimizing the squared error between reference points and a NURBS curve with respect to knots, control points, and weights. */
/* @param[in,out] fit zNURBSFitData */
zMat _zNURBSFitCalcJacobian(zNURBSFitData *fit)
{
  _zNURBSFitCalcSubJacobian( fit->nurbs, fit->fit_knots, fit->is_fixed_knots, fit->J_knot, fit->J_cp, fit->J_w );
  if( !fit->is_fixed_knots )
    zMatPutNC( fit->J, 0, 0, fit->J_knot );
  zMatPutNC( fit->J, 0, zMatColSizeNC(fit->J_knot), fit->J_cp );
  zMatPutNC( fit->J, 0, zMatColSizeNC(fit->J_knot) + zMatColSizeNC(fit->J_cp), fit->J_w );

  return fit->J;
}

/* Creates the Jacobian for minimizing the squared error between reference points and a NURBS curve with respect to knots, control points, and weights. */
/* @param[in,out] fit zNURBSFitData */
zMat _zNURBSFitCreateJacobian(zNURBSFitData *fit)
{
  int row_size, col_size;
  _zNURBSCreateSubJacobian( fit->nurbs, fit->fit_knots, &fit->J_knot, &fit->J_cp, &fit->J_w );
  row_size = zMatRowSizeNC(fit->J_knot);
  col_size = zMatColSizeNC(fit->J_knot) + zMatColSizeNC(fit->J_cp) + zMatColSizeNC(fit->J_w);
  fit->J = zMatAlloc( row_size, col_size );
  _zNURBSFitCalcJacobian( fit );

  return fit->J;
}

/* @param[in] nurbs zNURBS of knot with [start knot, end knot]  */
/* @param[out] knot clamped knot within [start knot, end knot] */
void _clamp_knot(const zNURBS* nurbs, double *knot )
{
  if( *knot < zNURBSKnotS(nurbs) )
    *knot = zNURBSKnotS(nurbs);
  if( *knot > zNURBSKnotE(nurbs) )
    *knot = zNURBSKnotE(nurbs);
}

/* Update fitting parameters knots, control points, and weights */
/* @param[in,out] fit zNURBSFitData */
bool _zNURBSFitUpdate(zNURBSFitData *fit)
{
  int i, skip, cp_num, p_size;
  double x, dx, exp_x_dx;
  /* update knots excluding the start(i=0) and end */
  if( !fit->is_fixed_knots ){
    for( i=1; i<fit->data_num-1; i++ ){
      x = zVecElemNC(fit->x, i);
      dx = zVecElemNC(fit->dx, i);
      zVecElemNC(fit->x, i) = x + dx;
      _clamp_knot( fit->nurbs, &zVecElemNC(fit->x, i) );
      zVecSetElemNC( fit->fit_knots, i, zVecElemNC(fit->x, i) );
    }
  }
  /* update control points excluding the start(i=0) and end */
  skip = fit->data_num;
  cp_num = zNURBSCPNum( fit->nurbs );
  p_size = zVecSize( zNURBSCP(fit->nurbs, 0) );
  for( i=1; i < cp_num-1; i++ ){
    zVecGetNC( fit->x,  skip + p_size * i, zNURBSCP( fit->nurbs, i ) );
    zVecGetNC( fit->dx, skip + p_size * i, fit->_cp_v );
    zVecAddDRC( zNURBSCP( fit->nurbs, i ), fit->_cp_v );
    zVecPutNC( fit->x, skip + p_size * i, zNURBSCP( fit->nurbs, i ) );
  }
  /* update weights */
  skip = fit->data_num + p_size * cp_num;
  for( i=0; i < cp_num; i++ ){
    x  = zVecElemNC(fit->x,  skip + i);
    dx = zVecElemNC(fit->dx, skip + i);
    /* negative weights are considered to be 0 (equivalent to being ignored). */
    if( x + dx < 0 || zIsTiny(x + dx) ){
      zVecSetElemNC( fit->x, skip + i, 0.0 );
      zNURBSSetWeight( fit->nurbs, i, 0.0 );
    } else{
      zVecSetElemNC( fit->x, skip + i, x + dx );
      /* multipling means w + dw -> exp(w) * exp(dw) */
      exp_x_dx = exp( x ) * exp( dx );
      if( zIsNan( exp_x_dx ) || zIsInf( exp_x_dx ) || exp_x_dx > HUGE_VALUE )
        zNURBSSetWeight( fit->nurbs, i, fabs(x+dx) );
      else
        zNURBSSetWeight( fit->nurbs, i, exp_x_dx );
    }
  }
  /* update fitting points on NURBS by new fitting parameters */
  for( i=0;i<fit->data_num; i++ ){
    zNURBSVec( fit->nurbs, zVecElemNC(fit->fit_knots, i), fit->_cp_v );
    zVecPutNC( fit->fit_vec, p_size * i, fit->_cp_v );
  }

  return true;
}

leParam* _createLinearEquationParameter(const int size, const int eval_size, leParam *le)
{
  le->a     = zMatAllocSqr( size );
  le->ident = zMatAllocSqr( size );
  zMatIdent( le->ident );
  le->b  = zVecAlloc( size );
  le->wn = zVecAlloc( size );
  le->we = zVecAlloc( eval_size );
  le->residual = zVecAlloc( size );
  le->err_vec = zVecAlloc( eval_size );
  zLEWorkspaceAlloc( &le->workspace, NULL, size );

  return le;
}

void _header_label_of_fitted_curve(FILE* fp);
void _fitted_curve(FILE *fp, const double s);
void _control_point_and_weight(FILE *fp);
void _target_and_fitted_point(FILE *fp);

/* create zNURBS and zNURBSFitData */
bool zNURBSFitCreate(const zVecArray *ref_vec_array, const int order, const int cp_num, const double start_knot, const double end_knot, zNURBS *nurbs, void **fit_data)
{
  int i, cp_size;
  zNURBSFitData *fit;

  *fit_data = (void *)( zAlloc( zNURBSFitData, 1 ) );
  fit = (zNURBSFitData *)(*fit_data);
  cp_size = zVecSize( *zArrayElemNC(ref_vec_array, 0) );

  fit->data_num = zArraySize( ref_vec_array );
  fit->dx       = zVecAlloc( fit->data_num + cp_size * cp_num + cp_num );
  fit->x        = zVecAlloc( zVecSize( fit->dx ) );
  fit->_cp_v    = zVecAlloc( cp_size );

  if( nurbs->cparray.size > 0 )
    zNURBSDestroy( nurbs );
  if( !zBSplineParamAlloc( &nurbs->param, order, cp_num, 0 ) )
    return false;
  zArrayAlloc( &nurbs->cparray, zNURBSCPCell, cp_num );
  if( zNURBSCPNum(nurbs) == 0 ){
    ZALLOCERROR();
    zNURBSDestroy( nurbs );
    return false;
  }
  zBSplineParamKnotInit( &nurbs->param );
  for( i=0; i < cp_num; i++ ){
    zNURBSSetWeight( nurbs, i, ZM_NURBS_DEFAULT_CP_WEIGHT );
    zNURBSCP(nurbs, i) = zVecAlloc( cp_size );
    zVecSetAll( zNURBSCP(nurbs, i), 0.0 ); /* @TODO given by ref_vec_array */
  }
  zNURBSKnotScale(nurbs, start_knot, end_knot);

  fit->nurbs          = nurbs;
  fit->fit_knots      = zVecAlloc( fit->data_num );
  fit->is_fixed_knots = false;

  zVecCreateFromzVecArray( ref_vec_array, &(fit->ref_vec) );
  fit->fit_vec = zVecAlloc( zVecSize( fit->ref_vec ) );

  _zNURBSFitCreateJacobian( fit );

  _createLinearEquationParameter( zVecSize(fit->x), zVecSize(fit->fit_vec), &fit->le );

  fit->print = zAlloc( zNURBSFitPrintFClass, 1 );
  fit->print->header_label_of_fitted_curve = _header_label_of_fitted_curve;
  fit->print->fitted_curve                 = _fitted_curve;
  fit->print->control_point_and_weight     = _control_point_and_weight;
  fit->print->target_and_fitted_point      = _target_and_fitted_point;

  return true;
}

void _initializeLinearEquationParameter(leParam * le)
{
  zVecSetAll( le->wn, TINY_VALUE );
  zVecSetAll( le->we, 1.0 );
}

/* initialize zNURBSFitData */
void zNURBSFitInitialize(void *fit_data, const zVec opt_fixed_knots, const zVecArray *init_cp_q_array)
{
  int i, cp_num, cp_size, skip;
  double s;
  zNURBSFitData *fit = (zNURBSFitData *)(fit_data);
  cp_num = zNURBSCPNum( fit->nurbs );
  cp_size = zVecSize( zNURBSCP(fit->nurbs, 0) );

  /* start and end of the fitting control points are bound on the start and end of the reference data. */
  zVecGetNC( fit->ref_vec, 0, zNURBSCP( fit->nurbs, 0 ) );
  zVecGetNC( fit->ref_vec, (fit->data_num-1) * cp_size, zNURBSCP( fit->nurbs, cp_num-1 ) );
  /* initial fitting knot */
  if( opt_fixed_knots != NULL && zVecSize( opt_fixed_knots ) == fit->data_num ){
    zVecCopyNC( opt_fixed_knots, fit->fit_knots );
    _zNURBSFitKnotJacobian( fit->nurbs, fit->fit_knots, fit->J_knot );
    zMatPutNC( fit->J, 0, 0, fit->J_knot );
    fit->is_fixed_knots = true;
  } else{
    for( i=0; i<fit->data_num; i++ ){
      s = (double)(i) / (double)(fit->data_num - 1);
      zVecElemNC( fit->fit_knots, i ) = (1 - s) * zNURBSKnotS(fit->nurbs) + s * zNURBSKnotE(fit->nurbs);
    }
    fit->is_fixed_knots = false;
  }
  /* start and end of the fitting knots are bound on the start and end of the NURBS knots. */
  zVecSetElemNC( fit->fit_knots, 0, zNURBSKnotS(fit->nurbs) );
  zVecSetElemNC( fit->fit_knots, zVecSize(fit->fit_knots)-1, zNURBSKnotE(fit->nurbs) );

  zVecPutNC( fit->x, 0, fit->fit_knots );
  skip = fit->data_num;
  for( i=0; i < cp_num; i++ ){
    if( init_cp_q_array != NULL &&
        zArraySize( init_cp_q_array ) == cp_num &&
        zVecSize( *zArrayElemNC( init_cp_q_array, i ) ) == cp_size &&
        i > 0 && i < cp_num-1 ){
      zVecCopy ( *zArrayElemNC( init_cp_q_array, i ), zNURBSCP(fit->nurbs, i) );
    }
    zVecPutNC( fit->x, skip + cp_size * i, zNURBSCP(fit->nurbs, i) );
  }
  skip = fit->data_num + cp_size * cp_num;
  for( i=0; i < cp_num; i++ ){
    zVecSetElemNC( fit->x, skip + i, zNURBSWeight(fit->nurbs, i) );
  }
  zVecSetAll( fit->dx, 1.0 );
  _initializeLinearEquationParameter( &fit->le );
  /**/
  _zNURBSFitUpdate( fit );
  zVecSubNC( fit->fit_vec, fit->ref_vec, fit->le.err_vec );
  fit->le.eval_dist = zVecInnerProd( fit->le.err_vec, fit->le.err_vec );
}

void zNURBSForceSetNURBS(void *fit_data, zNURBS *nurbs)
{
  zNURBSFitData *fit = (zNURBSFitData *)(fit_data);
  fit->nurbs = nurbs;
}

void _freeLinearEquationParameter(leParam *le)
{
  zLEWorkspaceFree( &le->workspace );
  zMatFree( le->a );
  zMatFree( le->ident );
  zVecFree( le->b );
  zVecFree( le->wn );
  zVecFree( le->we );
  zVecFree( le->residual );
  zVecFree( le->err_vec );
  le->eval_eq = 0.0;
  le->eval_dist = 0.0;
}

void zNURBSFitFree(void **fit_data)
{
  zNURBSFitData *fit = (zNURBSFitData *)(*fit_data);
  if( fit == NULL )
    return;
  zVecFree( fit->dx );
  zVecFree( fit->x );
  zVecFree( fit->_cp_v );
  zVecFree( fit->fit_vec );
  zVecFree( fit->ref_vec );
  zVecFree( fit->fit_knots );
  zMatFree( fit->J_knot );
  zMatFree( fit->J_cp );
  zMatFree( fit->J_w );
  zMatFree( fit->J );
  _freeLinearEquationParameter( &fit->le );
  zFree( fit->print );
  zFree( *fit_data );
  *fit_data = NULL;
}

zVec _zNURBSFitFunc(zVec dx, zVec err_vec, void* util)
{
  zNURBSFitData *fit = (zNURBSFitData *)util;
  zVecCopyNC( dx, fit->dx );
  _zNURBSFitUpdate( fit );
  zVecSubNC( fit->fit_vec, fit->ref_vec, err_vec );
  return err_vec;
}

zMat _zNURBSFitJacobian(zVec dx, zMat J, void* util)
{
  zNURBSFitData *fit = (zNURBSFitData *)util;
  _zNURBSFitCalcJacobian( fit );
  zMatCopyNC( fit->J, J );
  return J;
}

double _zNURBSFittingOne(zNURBSFitData *fit)
{
  double bias, size, *velm;
  zVec ret;
  leParam *le = &fit->le;
  _zNURBSFitFunc( fit->dx, le->err_vec, (void *)(fit) );
  _zNURBSFitJacobian( fit->dx, fit->J, (void *)(fit) );
  zMulMatTMat( fit->J, fit->J, le->a );
  bias = zTOL;
  bias += le->eval_eq;
  zMatCatDRC( le->a, bias + le->wn->buf[0], le->ident );
  zMulMatTVecNC( fit->J, le->err_vec, le->b );
  zVecMulDRC( le->b, -1.0 );
  ret = zLESolveSRBiasDST( le->a, le->b, le->wn, le->we, bias, fit->dx, &le->workspace );
  if( ret == NULL )
    return -1.0;
  /* dx */
  size = zVecSize( fit->dx );
  velm = zVecBufNC( fit->dx );
  for( ; size-->0; velm++ ){
    if( fabs(*velm) > HUGE_VALUE ) *velm = *velm/fabs(*velm) * HUGE_VALUE;
  }
  /**/
  zLEResidual( le->a, le->b, fit->dx, le->residual );
  /**/
  size = zVecSize( le->residual );
  velm = zVecBufNC( le->residual );
  for( ; size-->0; velm++ ){
    if( fabs(*velm) > HUGE_VALUE ) *velm = *velm/fabs(*velm) * HUGE_VALUE;
  }
  /**/
  le->eval_eq = zVecInnerProd( le->residual, le->residual );
  le->eval_dist = zVecInnerProd( le->err_vec, le->err_vec );

  return le->eval_dist;
}

/* main Fitting loop */
int zNURBSFitting(void *fit_data, const int iter, const double tol, FILE *fp)
{
  int i;
  double rest, converg;
  zNURBSFitData *fit = (zNURBSFitData *)fit_data;
  if( fp != NULL )
    fprintf( fp, "cnt eval_eq eval_dist\n" );
  converg = 0;
  for( rest=HUGE_VALUE, i=0; i<iter; i++ ){
    if( _zNURBSFittingOne( fit ) < 0.0 )
      return -1;
    converg = fit->le.eval_dist - rest;
    if( fp != NULL ){
      fprintf(fp, "%d %g %g\n", i, fit->le.eval_eq, fit->le.eval_dist );
      if( i%10 == 0 ){
        printf("%03d/%03d : ", i, iter);
        printf("eval_eq=%f\t eval_dist=%f\t eval_dist-rest=%f\r", fit->le.eval_eq, fit->le.eval_dist, converg ); fflush(stdout);
      }  /* monitor */
    }
    if( zIsTol( converg, tol ) ){
      break;
    }
    rest = fit->le.eval_dist;
  }
  if( fp != NULL ){
    printf("\n%03d/%03d : ", i, iter);
    printf("eval_eq=%f\t eval_dist=%f\t eval_dist-rest=%f\r", fit->le.eval_eq, fit->le.eval_dist, converg ); fflush(stdout);
  }

  return i;
}

/***********************************************************************************************/

static const zNURBSFitData *_g_fit;
const zNURBSFitPrintFClass *zNURBSFitPrintF(const void* fit_data)
{
  const zNURBSFitData *fit = (const zNURBSFitData *)(fit_data);
  _g_fit = fit;
  return fit->print;
}

void _header_label_of_fitted_curve(FILE *fp)
{
  int i, cp_size;
  cp_size = zVecSize( zNURBSCP(_g_fit->nurbs, 0) );
  fprintf( fp, "knot(s) " );
  for( i=0; i < cp_size; i++ )
    fprintf( fp, "fitted_nurbs_v[%d] ", i );
}

void _fitted_curve(FILE *fp, const double s)
{
  int i, cp_size;

  cp_size = zVecSize( zNURBSCP(_g_fit->nurbs, 0) );
  zNURBSVec( _g_fit->nurbs, s, _g_fit->_cp_v );
  fprintf( fp, "%g ", s );
  for( i=0; i < cp_size; i++ )
    fprintf( fp, "%g ", _g_fit->_cp_v->buf[i]);
  fprintf( fp, "\n" );
  fflush( fp );
}

void _control_point_and_weight(FILE *fp)
{
  int i, j, cp_size, cp_num;
  double weight;
  zVec cp;
  fprintf( fp, "No. " );
  cp_size = zVecSize( zNURBSCP(_g_fit->nurbs, 0) );
  for( i=0; i < cp_size; i++ )
    fprintf( fp, "fitted_nurbs_cp[%d] ", i );
  fprintf( fp, "fitted_weight \n" );
  cp_num = zNURBSCPNum( _g_fit->nurbs );
  for( i=0; i < cp_num; i++ ){
    cp = zNURBSCP( _g_fit->nurbs, i );
    weight = zNURBSWeight( _g_fit->nurbs, i );
    fprintf( fp, "%d ", i );
    for( j=0; j < cp_size; j++ )
      fprintf( fp, "%g ", cp->buf[j] );
    fprintf( fp, "%g \n", weight );
    fflush( fp );
  }
}

void _target_and_fitted_point(FILE *fp)
{
  int cp_size, i, j;
  double s;
  zVec vec;
  cp_size = zVecSize( zNURBSCP(_g_fit->nurbs, 0) );
  vec = zVecAlloc( cp_size );
  fprintf( fp, "No. fitted_knot(s) " );
  for( i=0; i < cp_size; i++ )
    fprintf( fp, "target_point[%d] ", i );
  for( i=0; i < cp_size; i++ )
    fprintf( fp, "fitted_nurbs_point[%d] ", i );
  fprintf( fp, "\n" );
  for( i=0; i < _g_fit->data_num; i++ ){
    s = zVecElemNC( _g_fit->fit_knots, i );
    zNURBSVec( _g_fit->nurbs, s, vec );
    fprintf( fp, "%d %g ", i, s );
    for( j=0; j < cp_size; j++ )
      fprintf( fp, "%g ", zVecElemNC( _g_fit->ref_vec, i * cp_size + j ) );
    for( j=0; j < cp_size; j++ )
      fprintf( fp, "%g ", vec->buf[j] );
    fprintf( fp, "\n" );
    fflush( fp );
  }
  zVecFree( vec );
}

void _zNURBSFitPrintFReport(const void *fit_data, FILE *fp)
{
  fprintf( fp, "\n----------------------------------------------\n");
  zNURBSFitData *fit = (zNURBSFitData *)fit_data;
  fprintf( fp, "J_knot =\n");
  zMatPrint( fit->J_knot );
  fprintf( fp, "J_cp =\n");
  zMatPrint( fit->J_cp );
  fprintf( fp, "J_w =\n");
  zMatPrint( fit->J_w );
  fprintf( fp, "J = \n");
  zMatPrint( fit->J );
  fprintf( fp, "x =\n");
  zVecPrint(fit->x);
  fprintf( fp, "dx =\n");
  zVecPrint(fit->dx);
  fprintf( fp, "fit_knots=\n");
  zVecPrint(fit->fit_knots);
  fprintf( fp, "CP=\n");
  zNURBSCPFPrint(stdout, fit->nurbs);
  fprintf( fp, "fit->le.err_vec =");
  zVecPrint(fit->le.err_vec);
  fprintf( fp, "ref_vec=\n");
  zVecPrint(fit->ref_vec);
  fprintf( fp, "fit_vec=\n");
  zVecPrint(fit->fit_vec);
  fprintf( fp, "\n fit->le.err_vec =");
  zVecPrint(fit->le.err_vec);
  fprintf( fp, "\n");
  fprintf( fp, "eval_eq = %f\n", fit->le.eval_eq );
  fprintf( fp, "eval_dist = %f\n", fit->le.eval_dist );
  fprintf( fp, "\n----------------------------------------------\n");
}

