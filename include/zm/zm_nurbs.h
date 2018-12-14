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
 * \brief cell of an array of control points for NURBS.
 *
 * zNURBSCPCell is a cell of an array of control points for NURBS.
 * It contains a control point and a weight put on it.
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
 * zNURBS is an NURBS curve made from a sequence of control points.
 *//* ******************************************************* */
typedef struct{
  int dim;   /*!< \brief dimension of a curve */
  /*! \cond */
  zVec knot;             /* knot vector */
  zNURBSCPArray cparray; /* an array of control points */
  /*! \endcond */
} zNURBS;

#define zNURBSKnot(n,i)    zVecElem((n)->knot,i)
#define zNURBSWeight(n,i)  ( zArrayElem(&(n)->cparray,i)->w )
#define zNURBSCP(n,i)      ( zArrayElem(&(n)->cparray,i)->cp )

#define zNURBSKnotNum(n)   ( zVecSizeNC((n)->knot) - 1 )
#define zNURBSKnot0(n)     zNURBSKnot(n,0)
#define zNURBSKnotE(n)     zNURBSKnot(n,zNURBSKnotNum(n)-1)

#define zNURBSCPNum(n)     zArrayNum( &(n)->cparray )

/*! \brief create a NURBS curve.
 *
 * zNURBSCreate() creates a NURBS curve \a nurbs from a given sequence
 * of control points. \a seq provides the control points.
 * \a dim is the dimension of the curve, which has to be less than the
 * size of \a seq.
 * It is initialized as a uniform Bezier spline curve with fixed
 * boundary points. The weights on each control point and the knot
 * vector can be modified later.
 * \return
 * zNURBSCreate() returns the true value if it succeeds to create a NURBS
 * curve. If \a dim is larger than the size of \a seq plus one or it
 * fails to allocate internal workspace, the false value is returned.
 */
__EXPORT bool zNURBSCreate(zNURBS *nurbs, zSeq *seq, int dim);

/*! \brief destroy a NURBS curve.
 *
 * zNURBSDestroy() destroys a NURBS curve \a nurbs.
 */
__EXPORT void zNURBSDestroy(zNURBS *nurbs);

/*! \brief normalize the knot vector of a NURBS curve.
 *
 * zNURBsKnotNormalize() normalizes the knot vector of a NURBS curve
 * \a nurbs so that it starts from 0 and ends at 1.
 */
__EXPORT void zNURBSKnotNormalize(zNURBS *nurbs);

/*! \brief zNURBSVec
 *
 * zNURBSVec() computes a vector on a NURBS curve \a nurbs with respect
 * to the parameter \a t. The computed vector is stored where \a v points.
 * \return
 * zNURBSVec() returns a pointer \a v if \a t is valid. Otherwise, the
 * null vector is returned.
 */
__EXPORT zVec zNURBSVec(zNURBS *nurbs, double t, zVec v);

/*! \brief zNURBSVecDiff
 *
 * zNURBSVecDiff() computes the derivative of a NURBS curve \a nurbs at
 * the parameter \a t and puts it into \a v.
 * \a diff is the number of differential.
 * \return
 * zNURBSVecDiff() returns a pointer \a v if it succeeds to compute the
 * derivative. If \a diff is invalid or it fails to allocate the internal
 * workspace, the false value is returned.
 */
__EXPORT zVec zNURBSVecDiff(zNURBS *nurbs, double t, zVec v, int diff);

/*! \brief nearest neighbor on NURBS
 *
 * zNURBSVecNN() finds the nearest-neighbor vector on a NURBS curve
 * defined by \a nurbs from a vector \a v. The result is put into \a nn.
 * \return
 * zNURBSVecNN() returns the parameter corresponding to the nearest-
 * neighbor vector found by this function.
 */
__EXPORT double zNURBSVecNN(zNURBS *nurbs, zVec v, zVec nn);

/* for debug */

__EXPORT void zNURBSKnotFWrite(FILE *fp, zNURBS *nurbs);
__EXPORT void zNURBSCPFWrite(FILE *fp, zNURBS *nurbs);

__END_DECLS

#endif /* __ZM_NURBS_H__ */
