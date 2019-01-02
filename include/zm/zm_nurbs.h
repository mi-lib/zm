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
/*! \struct zNURBS1
 * \brief NURBS curve.
 *
 * zNURBS1 is a 1-dimensional NURBS curve made from a sequence
 * of control points in n-dimensional space.
 *//* ******************************************************* */
typedef struct{
  int dim;   /*!< \brief dimension of a curve */
  /*! \cond */
  zVec knot;             /* knot vector */
  zNURBSCPArray cparray; /* an array of control points */
  /*! \endcond */
} zNURBS1;

#define zNURBS1Knot(n,i)    zVecElemNC((n)->knot,i)
#define zNURBS1Weight(n,i)  ( zArrayElemNC(&(n)->cparray,i)->w )
#define zNURBS1CP(n,i)      ( zArrayElemNC(&(n)->cparray,i)->cp )

#define zNURBS1KnotNum(n)   ( zVecSizeNC((n)->knot) - 1 )
#define zNURBS1Knot0(n)     zNURBS1Knot(n,0)
#define zNURBS1KnotE(n)     zNURBS1Knot(n,zNURBS1KnotNum(n)-1)

#define zNURBS1CPNum(n)     zArraySize( &(n)->cparray )

/*! \brief create a 1-dimensional NURBS curve.
 *
 * zNURBS1Create() creates a 1-dimensional NURBS curve \a nurbs
 * from a given sequence of control points. \a seq provides the
 * control points.
 * \a dim is the dimension of the curve, which has to be less than
 * the size of \a seq.
 * It is initialized as a uniform Bezier spline curve with fixed
 * boundary points. The weights on each control point and the knot
 * vector can be modified later.
 * \return
 * zNURBS1Create() returns the true value if it succeeds to create
 * the NURBS curve. If \a dim is larger than the size of \a seq
 * plus one or it fails to allocate internal workspace, the false
 * value is returned.
 */
__EXPORT bool zNURBS1Create(zNURBS1 *nurbs, zSeq *seq, int dim);

/*! \brief destroy a 1-dimensional NURBS curve.
 *
 * zNURBS1Destroy() destroys a 1-dimensional NURBS curve \a nurbs.
 */
__EXPORT void zNURBS1Destroy(zNURBS1 *nurbs);

/*! \brief normalize the knot vector of a 1-dimensional NURBS curve.
 *
 * zNURBS1KnotNormalize() normalizes the knot vector of a
 * 1-dimensional NURBS curve \a nurbs so that it starts from 0 and
 * ends at 1.
 */
__EXPORT void zNURBS1KnotNormalize(zNURBS1 *nurbs);

/*! \brief compute a vector on 1-dimensional NURBS curve.
 *
 * zNURBS1Vec() computes a vector on a 1-dimensional NURBS curve
 * \a nurbs with respect to a given parameter \a t. The computed
 * vector is stored where \a v points.
 * \return
 * zNURBS1Vec() returns a pointer \a v if \a t is valid. Otherwise,
 * the null vector is returned.
 */
__EXPORT zVec zNURBS1Vec(zNURBS1 *nurbs, double t, zVec v);

/*! \brief find the derivative of a 1-dimensional NURBS curve.
 *
 * zNURBS1VecDiff() computes the derivative of a 1-dimensional
 * NURBS curve \a nurbs at a given parameter \a t and puts it into
 * \a v.
 * \a diff is the number of differential.
 * \return
 * zNURBS1VecDiff() returns a pointer \a v if it succeeds to compute
 * the derivative. If \a diff is invalid or it fails to allocate the
 * internal workspace, the false value is returned.
 */
__EXPORT zVec zNURBS1VecDiff(zNURBS1 *nurbs, double t, zVec v, int diff);

/*! \brief nearest neighbor on a 1-dimensional NURBS curve.
 *
 * zNURBS1VecNN() finds the nearest-neighbor vector on a 1-dimensional
 * NURBS curve defined by \a nurbs from a vector \a v. The result
 * is put into \a nn.
 * \return
 * zNURBS1VecNN() returns the parameter corresponding to the nearest-
 * neighbor vector found by this function.
 */
__EXPORT double zNURBS1VecNN(zNURBS1 *nurbs, zVec v, zVec nn);

/* for debug */

__EXPORT void zNURBS1KnotFWrite(FILE *fp, zNURBS1 *nurbs);
__EXPORT void zNURBS1CPFWrite(FILE *fp, zNURBS1 *nurbs);

__END_DECLS

#endif /* __ZM_NURBS_H__ */
