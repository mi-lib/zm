/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_nurbs - NURBS curve.
 */

#ifndef __ZM_NURBS_H__
#define __ZM_NURBS_H__

#include <zm/zm_seq.h>

__BEGIN_DECLS

/* ********************************************************** */
/*! \struct zBSplineParam
 * \brief B-spline parameter
 *
 * zBSplineParam defines the parameter space of B-spline.
 * It consists of the order of the curve and knots.
 *//* ******************************************************* */
ZDEF_STRUCT( __ZM_CLASS_EXPORT, zBSplineParam ){
  int order; /*!< \brief order of a curve */
  zVec knot; /*!< \brief knot vector */
  int slice; /*!< \brief number of slices */
};

#define ZM_BSPLINE_DEFAULT_SLICE_NUM 20

#define zBSplineParamKnotNum(param)       zVecSizeNC((param)->knot)
#define zBSplineParamKnot(param,i)        zVecElemNC((param)->knot,i)
#define zBSplineParamSetKnot(param,i,v)   ( zBSplineParamKnot(param,i) = (v) )
#define zBSplineParamCPNum(param)         ( zBSplineParamKnotNum(param) - (param)->order - 1 )
#define zBSplineParamKnotS(param)         zBSplineParamKnot( param, (param)->order )
#define zBSplineParamKnotE(param)         zBSplineParamKnot( param, zBSplineParamCPNum(param) )
#define zBSplineParamKnotSlice(param,k)   ( ( zBSplineParamKnotE(param) - zBSplineParamKnotS(param) ) * (k) / ((param)->slice) + zBSplineParamKnotS(param) )
#define zBSplineParamSetSlice(param,s)    ( (param)->slice = (s) )

/*! \brief allocate B-spline parameter. */
__ZM_EXPORT zBSplineParam *zBSplineParamAlloc(zBSplineParam *param, int order, int nc, int slice);

/*! \brief free B-spline parameters. */
__ZM_EXPORT void zBSplineParamFree(zBSplineParam *param);

/*! \brief initialize knots of a B-spline parameter. */
__ZM_EXPORT void zBSplineParamKnotInit(zBSplineParam *param);

/*! \brief normalize knot vector of a NURBS curve. */
__ZM_EXPORT void zBSplineParamKnotNormalize(zBSplineParam *param);

/*! \brief find a knot segment that includes the given parameter. */
__ZM_EXPORT int zBSplineParamSeg(zBSplineParam *param, double t);

/*! \brief basis function of B-spline family. */
__ZM_EXPORT double zBSplineParamBasis(zBSplineParam *param, double t, int i, int r, int seg);

/*! \brief derivative of the basis function of B-spline family. */
__ZM_EXPORT double zBSplineParamBasisDiff(zBSplineParam *param, double t, int i, int r, int seg, int diff);

/* ********************************************************** */
/*! \struct zNURBSCPCell
 * \brief cell of NURBS containing a control point and weight.
 *
 * zNURBSCPCell is a cell of NURBS that contains a control point
 * and associated weight.
 *//* ******************************************************* */
ZDEF_STRUCT( __ZM_CLASS_EXPORT, zNURBSCPCell ){
  zVec cp;  /*!< control point */
  double w; /*!< weight */
};

/* ********************************************************** */
/*! \struct zNURBSCPArray
 * \brief an array of control points for NURBS.
 *
 * zNURBSCPArray is an array of control points for NURBS.
 * It is defined with a macro zArrayClass.
 * \sa zArrayClass.
 *//* ******************************************************* */
zArrayClass( zNURBSCPArray, zNURBSCPCell );

/* ********************************************************** */
/*! \struct zNURBS
 * \brief NURBS curve.
 *
 * zNURBS is a NURBS curve made from a sequence of control
 * points in n-dimensional space.
 *//* ******************************************************* */
ZDEF_STRUCT( __ZM_CLASS_EXPORT, zNURBS ){
  zBSplineParam param;   /*!< \brief B-spline parameter */
  zNURBSCPArray cparray; /*!< \brief an array of control points */
};

#define zNURBSKnotNum(nurbs)       zBSplineParamKnotNum( &(nurbs)->param )
#define zNURBSKnot(nurbs,i)        zBSplineParamKnot( &(nurbs)->param, i )
#define zNURBSSetKnot(nurbs,i,v)   zBSplineParamSetKnot( &(nurbs)->param, i, v )
#define zNURBSKnotS(nurbs)         zBSplineParamKnotS( &(nurbs)->param )
#define zNURBSKnotE(nurbs)         zBSplineParamKnotE( &(nurbs)->param )
#define zNURBSKnotSlice(nurbs,k)   zBSplineParamKnotSlice( &(nurbs)->param, k )
#define zNURBSSlice(nurbs)         (nurbs)->param.slice
#define zNURBSSetSlice(nurbs,s)    zBSplineParamSetSlice( &(nurbs)->param, s )

#define zNURBSCPNum(nurbs)         zArraySize( &(nurbs)->cparray )
#define zNURBSWeight(nurbs,i)      ( zArrayElemNC(&(nurbs)->cparray,i)->w )
#define zNURBSSetWeight(nurbs,i,v) ( zNURBSWeight(nurbs,i) = (v) )
#define zNURBSCP(nurbs,i)          ( zArrayElemNC(&(nurbs)->cparray,i)->cp )
#define zNURBSSetCP(nurbs,i,v)     zVecCopy( v, zNURBSCP(nurbs,i) )

#define ZM_NURBS_DEFAULT_CP_WEIGHT 1.0

/*! \brief create a NURBS curve.
 *
 * zNURBSCreate() creates a NURBS curve \a nurbs from a given sequence of control points.
 * \a seq provides the control points.
 * \a order is the order of the curve, which has to be less than the size of \a seq.
 *
 * The knots are initialized as to be a uniform B-spline curve with fixed boundary points.
 * The weights on each control point and the knot vector can be modified later.
 * \return
 * zNURBSCreate() returns the true value if it succeeds to create the NURBS curve. If
 * \a order is larger than the size of \a seq plus one or it fails to allocate internal
 * workspace, the false value is returned.
 */
__ZM_EXPORT bool zNURBSCreate(zNURBS *nurbs, zSeq *seq, int order);

/*! \brief destroy a NURBS curve.
 *
 * zNURBSDestroy() destroys a NURBS curve \a nurbs.
 */
__ZM_EXPORT void zNURBSDestroy(zNURBS *nurbs);

/*! \brief normalize the knot vector of a NURBS curve.
 *
 * zNURBSKnotNormalize() normalizes the knot vector of a
 * NURBS curve \a nurbs so that it starts from 0 and ends at 1.
 */
#define zNURBSKnotNormalize(nurbs) zBSplineParamKnotNormalize( &(nurbs)->param )

/*! \brief compute a vector on NURBS curve.
 *
 * zNURBSVec() computes a vector on a NURBS curve \a nurbs with
 * respect to a given parameter \a t. The computed vector is stored
 * where \a v points.
 * \return
 * zNURBSVec() returns a pointer \a v if \a t is valid. Otherwise,
 * the null vector is returned.
 */
__ZM_EXPORT zVec zNURBSVec(zNURBS *nurbs, double t, zVec v);

/*! \brief find the derivative of a NURBS curve.
 *
 * zNURBSVecDiff() computes the derivative of a NURBS curve \a nurbs
 * at a given parameter \a t and puts it into \a v.
 * \a diff is the order of derivative.
 * \return
 * zNURBSVecDiff() returns a pointer \a v if it succeeds to compute
 * the derivative. If \a diff is invalid or it fails to allocate the
 * internal workspace, the false value is returned.
 */
__ZM_EXPORT zVec zNURBSVecDiff(zNURBS *nurbs, double t, int diff, zVec v);

/*! \brief nearest neighbor on a NURBS curve.
 *
 * zNURBSVecNN() finds the nearest-neighbor vector on a NURBS curve
 * defined by \a nurbs from a vector \a v. The result is put into
 * \a nn.
 * \return
 * zNURBSVecNN() returns the parameter corresponding to the nearest-
 * neighbor vector found by this function.
 */
__ZM_EXPORT double zNURBSVecNN(zNURBS *nurbs, zVec v, zVec nn);

/* for debug */

#define zNURBSKnotFPrint(fp,nurbs) zVecFPrint( fp, (nurbs)->knot )

__ZM_EXPORT void zNURBSCPFPrint(FILE *fp, zNURBS *nurbs);

__END_DECLS

#endif /* __ZM_NURBS_H__ */
