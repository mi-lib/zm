/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_intg.h
 * \brief numerical integrator.
 * \author Zhidao
 */

#ifndef __ZM_INTG_H__
#define __ZM_INTG_H__

#include <zm/zm_misc.h>

__BEGIN_DECLS

/* ********************************************************** */
/*! \defgroup intg numerical integrator.
 * \{ *//* ************************************************** */

/* ********************************************************** */
/*! \brief numerical integrator class.
 *//* ******************************************************* */
typedef struct{
  double s;  /*!< integral value */
  double x0; /*!< previous value */

  /*! \cond */
  double _dt, _x; /* for quadratic approximation */
  /*! \endcond */
} zIntg;

/*! \brief initialize numerical integrator.
 *
 * zIntgInit() initializes a numerical integrator \a intg.
 * \a s0 is the initial integral value, and \a x0 is the
 * initial function value to be integrated.
 */
__EXPORT void zIntgInit(zIntg *intg, double s0, double x0);

/*! \brief rectangular integration.
 *
 * zIntgRect() integrates sequential values based on
 * rectangular formula. \a x is the latest function value
 * and \a dt is the discrete time step.
 */
__EXPORT double zIntgRect(zIntg *intg, double x, double dt);

/*! \brief trapezoidal integration.
 *
 * zIntgTR() integrates sequential values based on
 * trapezoidal formula. \a x is the latest function value
 * and \a dt is the discrete time step.
 */
__EXPORT double zIntgTR(zIntg *intg, double x, double dt);

/*! \brief quadratic approximation integration.
 *
 * zIntgQuad() integrates sequential values based on
 * quadratic approximation of the function. \a x is the latest
 * function value and \a dt is the discrete time step.
 */
__EXPORT double zIntgQuad(zIntg *intg, double x, double dt);

/*! \} */

__END_DECLS

#endif /* __ZM_INTG_H__ */
