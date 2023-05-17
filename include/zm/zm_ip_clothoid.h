/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_ip_clothoid.h
 * \brief interpolation: clothoid curve interpolation.
 * \author Zhidao
 */

#ifndef __ZM_IP_CLOTHOID_H__
#define __ZM_IP_CLOTHOID_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zClothoid
 * clothoid curve interpolator class
 * ********************************************************** */

typedef struct{
  double x0, y0; /*!< \brief initial point */
  double f0, f1; /*!< \brief initial and terminal angle */
  double fc;     /*!< \brief curvature */
  double fs;     /*!< \brief shrinkage */
  double _h;     /*!< \brief length of segment */
} zClothoid;

/*! \brief create a segment of clothoid curve.
 *
 * zClothoidCreateSegment() creates a segment of clothoid curve \a cl that
 * connects (\a x0, \a y0) and (\a x1, \a y1) with the initial and terminal
 * angles \a f0 and \a f1, respectively.
 * \return
 * zClothoidCreateSegment() returns a pointer \a cl.
 */
__ZM_EXPORT zClothoid *zClothoidCreateSegment(zClothoid *cl, double x0, double y0, double f0, double x1, double y1, double f1);

/*! \brief x-y values of a clothoid curve.
 *
 * zClothoidXY() finds a coordinate on a clothoid curve segment \a cl at
 * a parameter \a s. The result is stored into (\a x, \a y).
 * \return
 * zClothoidXY() returns the true value if it succeeds to compute the values.
 * Otherwise, the false value is returned.
 */
__ZM_EXPORT bool zClothoidXY(zClothoid *cl, double s, double *x, double *y);

__END_DECLS

#endif /* __ZM_IP_CLOTHOID_H__ */
