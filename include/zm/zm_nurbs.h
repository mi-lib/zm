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
/*! \struct zNURBSCPCell
 * \brief cell of NURBS containing a control point and weight.
 *
 * zNURBSCPCell is a cell of NURBS that contains a control point
 * and associated weight.
 *//* ******************************************************* */
typedef struct{
  zVec cp;  /*!< control point */
  double w; /*!< weight */
} zNURBSCPCell;

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
typedef struct{
  int dim;   /*!< \brief dimension of a curve */
  /*! \cond */
  zVec knot; /* knot vector */
  zNURBSCPArray cparray; /* an array of control points */
  /*! \endcond */
} zNURBS;

#define zNURBSKnotNum(n)       zVecSizeNC((n)->knot)
#define zNURBSKnot(n,i)        zVecElemNC((n)->knot,i)
#define zNURBSSetKnot(n,i,v)   ( zNURBSKnot(n,i) = (v) )
#define zNURBSKnotS(n)         zNURBSKnot(n,(n)->dim)
#define zNURBSKnotE(n)         zNURBSKnot(n,zNURBSCPNum(n))
#define zNURBSKnotSlice(n,k,s) ( ( zNURBSKnotE(n) - zNURBSKnotS(n) ) * k / s + zNURBSKnotS(n) )

#define zNURBSCPNum(n)         zArraySize( &(n)->cparray )
#define zNURBSWeight(n,i)      ( zArrayElemNC(&(n)->cparray,i)->w )
#define zNURBSSetWeight(n,i,v) ( zNURBSWeight(n,i) = (v) )
#define zNURBSCP(n,i)          ( zArrayElemNC(&(n)->cparray,i)->cp )
#define zNURBSSetCP(n,i,v)     zVecCopy( v, zNURBSCP(n,i) )

/*! \brief create a NURBS curve.
 *
 * zNURBSCreate() creates a NURBS curve \a nurbs from a given
 * sequence of control points. \a seq provides the control points.
 * \a dim is the dimension of the curve, which has to be less than
 * the size of \a seq.
 * It is initialized as a uniform Bezier spline curve with fixed
 * boundary points. The weights on each control point and the knot
 * vector can be modified later.
 * \return
 * zNURBSCreate() returns the true value if it succeeds to create
 * the NURBS curve. If \a dim is larger than the size of \a seq
 * plus one or it fails to allocate internal workspace, the false
 * value is returned.
 */
__EXPORT bool zNURBSCreate(zNURBS *nurbs, zSeq *seq, int dim);

/*! \brief destroy a NURBS curve.
 *
 * zNURBSDestroy() destroys a NURBS curve \a nurbs.
 */
__EXPORT void zNURBSDestroy(zNURBS *nurbs);

/*! \brief normalize the knot vector of a NURBS curve.
 *
 * zNURBSKnotNormalize() normalizes the knot vector of a
 * NURBS curve \a nurbs so that it starts from 0 and ends at 1.
 */
__EXPORT void zNURBSKnotNormalize(zNURBS *nurbs);

/*! \brief compute a vector on NURBS curve.
 *
 * zNURBSVec() computes a vector on a NURBS curve \a nurbs with
 * respect to a given parameter \a t. The computed vector is stored
 * where \a v points.
 * \return
 * zNURBSVec() returns a pointer \a v if \a t is valid. Otherwise,
 * the null vector is returned.
 */
__EXPORT zVec zNURBSVec(zNURBS *nurbs, double t, zVec v);

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
__EXPORT zVec zNURBSVecDiff(zNURBS *nurbs, double t, int diff, zVec v);

/*! \brief nearest neighbor on a NURBS curve.
 *
 * zNURBSVecNN() finds the nearest-neighbor vector on a NURBS curve
 * defined by \a nurbs from a vector \a v. The result is put into
 * \a nn.
 * \return
 * zNURBSVecNN() returns the parameter corresponding to the nearest-
 * neighbor vector found by this function.
 */
__EXPORT double zNURBSVecNN(zNURBS *nurbs, zVec v, zVec nn);

/* for debug */

#define zNURBSKnotFPrint(fp,n) zVecFPrint( fp, (n)->knot )

__EXPORT void zNURBSCPFPrint(FILE *fp, zNURBS *nurbs);

__END_DECLS

#endif /* __ZM_NURBS_H__ */
