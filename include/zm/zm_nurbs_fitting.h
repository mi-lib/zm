/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 2025, Daishi Kaneta (kdaic), Tomomichi Sugihara (Zhidao),
 *
 * zm_nurbs_fitting - NURBS curve fitting with Gauss Newton algorighm
 *
 * This algorighm is based with reference to the following papers.
 * Gulliksson, Mårten E. and Nils Carlsson,
 * “Surface Fitting with NURBS - a Gauss Newton with Trust Region Approach",
 * 13th WSEAS International Conference on APPLIED MATHEMATICS Puerto de la Cruz,
 * SPAIN, DEC 15-17, 2008.
 *
 */

#ifndef __ZM_NURBS_FITTING_H__
#define __ZM_NURBS_FITTING_H__

#include <zm/zm_nurbs.h>

/* @param[in] ref_vec_array vector array of reference points */
/* @param[in] order         order of zNURBS */
/* @param[in] cp_num        the number of control point  */
/* @param[in] start_knot    start knot */
/* @param[in] end_knot      end knot */
/* @param[in,out] nurbs     zNURBS for fitting which the pointer is registerd into fit_data */
/* @param[out] fit_data     the main zNURBSFitData to be created, it has pointer of nurbs */
__ZM_EXPORT bool zNURBSFitCreate(const zVecArray *ref_vec_array, const int order, const int cp_num, const double start_knot, const double end_knot, zNURBS *nurbs, void **fit_data);

/* @param[in,out] fit_data        the main data to be initialized */
/* @param[in]     opt_fixed_knots (option) if this fixed_knots is not NULL, the internal fit_knots is not updated and fixed with this input argument, else if fixed_knots is NULL, the internal fit_knots is update and used as fitting(searching) parameter */
__ZM_EXPORT void zNURBSFitInitialize(void *fit_data, const zVec opt_fixed_knots, const zVecArray *init_cp_q_array);

/* @param[in,out] fit_data the main data which NURBS will be set */
/* @param[in]     nurbs    a NURBS pointer to be set into fit_data */
__ZM_EXPORT void zNURBSForceSetNURBS(void *fit_data, zNURBS *nurbs);

/* @param[in,out] fit_data the main data to be fitted the NURBS */
__ZM_EXPORT void zNURBSFitFree(void **fit_data);

/* @param[in,out] fit_data the main data for fitting */
/* @param[in]     iter     iteration number for iterative calculation */
/* @param[in]     tol      tolerance to check precise of error */
/* @param[in,out] fp       pointer of the FILE */
__ZM_EXPORT int zNURBSFitting(void *fit_data, const int iter, const double tol, FILE *fp);

/* print class */
ZDEF_STRUCT( __ZM_CLASS_EXPORT, zNURBSFitPrintFClass ){
  void (*header_label_of_fitted_curve)(FILE*);
  void (*fitted_curve)(FILE*, const double); /* 2nd argument is a knot to pop-out */
  void (*control_point_and_weight)(FILE*);
  void (*target_and_fitted_point)(FILE*);
};
/* @param[in] fit_data the main data to be filed out */
__ZM_EXPORT const zNURBSFitPrintFClass *zNURBSFitPrintF(const void *fit_data);

/* @param[in] fit_data the main data to be print out */
/* @param[in,out] fp   a FILE pointer (i.e. stdout) to print the internal fitting variables */
__ZM_EXPORT void _zNURBSFitPrintFReport(const void *fit_data, FILE *fp);

#endif
