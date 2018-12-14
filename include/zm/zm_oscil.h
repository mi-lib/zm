/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_oscil.h
 * \brief nonlinear oscillator.
 * \author Zhidao
 */

#ifndef __ZM_OSCIL_H__
#define __ZM_OSCIL_H__

#include <zm/zm_misc.h>

__BEGIN_DECLS

/* ********************************************************** */
/*! \defgroup oscil nonlinear oscillator.
 * \{ *//* ************************************************** */

typedef struct{
  void (* init)(void*,double,double);     /*!< initializer */
  double (* update)(void*,double,double); /*!< update */
  double (* output)(void*);               /*!< oscillator output */
} zOscCom;

/* ********************************************************** */
/*! \brief a generalized oscillator class.
 *//* ******************************************************* */
typedef struct{
  void *prp;    /*!< oscillator properties */
  zOscCom *com; /*!< command set */
} zOsc;

/*! \brief destroy oscillator.
 *
 * zOscDestroy() destroys an oscillator instance pointed by \a osc.
 */
#define zOscDestroy(osc) do{\
  zFree( (osc)->prp );\
  (osc)->com = NULL;\
} while(0)

/*! \brief initialize oscillator.
 *
 * zOscInit() initializes an oscillator instance pointed by \a osc
 * such that the initial state is ( \a x0, \a v0 ).
 */
#define zOscInit(osc,x0,v0)  (osc)->com->init( (osc)->prp, x0, v0 )

/*! \brief update the internal state of an oscillator.
 *
 * zOscUpdate() updates the internal state of an oscillator instance
 * pointed by \a osc. \a u is an external input and \a dt is a
 * quantized update time.
 */
#define zOscUpdate(osc,u,dt) (osc)->com->update( (osc)->prp, u, dt )

/*! \brief output an oscillator value.
 *
 * zOscOutput() outputs a value from an oscillator instance
 * pointed by \a osc.
 */
#define zOscOutput(osc)      (osc)->com->output( (osc)->prp )

/* ********************************************************** */
/*! \brief van der Pol's oscillator.
 *//* ******************************************************* */

/*! \brief create van der Pol's oscillator.
 *
 * zOscCreateVDP() creates a van der Pol's oscillator pointed
 * by \a osc. \a t is the natural oscillation term,
 * \a damp is the damping factor, and \a amp is the output amplitude.
 */
__EXPORT zOsc *zOscCreateVDP(zOsc *osc, double t, double damp, double amp);

/* ********************************************************** */
/*! \brief Kuramoto's phase oscillator.
 *//* ******************************************************* */

/*! \brief create Kuramoto's phase oscillator.
 *
 * zOscCreateKura() creates a Kuramotof's phase oscillator
 * pointed by \a osc. \a wn is the natural oscillation frequency,
 * \a ke is the entrainment gain, and \a po is the phase offset.
 */
__EXPORT zOsc *zOscCreateKura(zOsc *osc, double wn, double ke, double po);

/*! \} */

__END_DECLS

#endif /* __ZM_OSCIL_H__ */
