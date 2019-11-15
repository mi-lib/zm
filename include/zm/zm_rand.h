/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_rand.h
 * \brief random number generator.
 *
 * The implementation of Mersenne twister proposed by
 * M. Matsumoto and T. Nishimura (1995) is a simple
 * rearrangement of the code written by Mr. Isaku Wada.
 * The original sources zmtrand.{h,c} are available at:
 *  http://www001.upp.so-net.ne.jp/isaku/index.html
 * \author Zhidao
 */

#ifndef __ZM_RAND_H__
#define __ZM_RAND_H__

#include <zm/zm_misc.h>

__BEGIN_DECLS

/* ********************************************************** */
/*! \defgroup rand random number generator with Mersenne twister.
 * \{ *//* ************************************************** */

/* ********************************************************** */
/* a variety of random distributions
 * ********************************************************** */

/* gamma distribution family */

/*! \brief a random number yielded from gamma distribution. */
__EXPORT double zRandGamma(zRandMT *mt, double a);
/*! \brief a random number yielded from Chi squared distribution. */
__EXPORT double zRandChiSqr(zRandMT *mt, double nu);
/*! \brief a random number yielded from beta distribution. */
__EXPORT double zRandBeta(zRandMT *mt, double a, double b);
/*! \brief a random number yielded from F distribution. */
__EXPORT double zRandFD(zRandMT *mt, double a, double b);

/* normal distribution family */

/*! \brief a random number yielded from normal distribution based on Box-Muller's method. */
__EXPORT double zRandND0(zRandMT *mt);
/*! \brief a random number yielded from normal distribution based on Box-Muller's method.
 *
 * zRandND() returns a randomly generated value in accordance with
 * the normal distribution with the mean \a mu and the variance \a sigma
 * based on Box-Muller's method.
 */
__EXPORT double zRandND(zRandMT *mt, double mu, double sigma);
/*! \brief a random number yielded from t distribution. */
__EXPORT double zRandT(zRandMT *mt, double n);

/* continuous distribution */

/*! \brief a random number yielded from exponential distribution
 * with unit meanvalue. */
__EXPORT double zRandExp(zRandMT *mt);
/*! \brief a random number yielded from logistic distribution. */
__EXPORT double zRandLog(zRandMT *mt);
/*! \brief a random number yielded from power distribution. */
__EXPORT double zRandPower(zRandMT *mt, double n);
/*! \brief a random number yielded from triangle distribution. */
__EXPORT double zRandTri(zRandMT *mt);
/*! \brief a random number yielded from Cauchy's distribution. */
__EXPORT double zRandCauchy(zRandMT *mt);
/*! \brief a random number yielded from Weibull distribution. */
__EXPORT double zRandWeibull(zRandMT *mt, double alpha);

/* discrete distribution */

/*! \brief a random number yielded from binomial distribution. */
__EXPORT int zRandBinom(zRandMT *mt, int n, double p);
/*! \brief a random number yielded from two-variable binomial
 * distribution with relation. */
__EXPORT void zRandBinom2(zRandMT *mt, double r, double *x, double *y);
/*! \brief a random number yielded from Poisson's distribution. */
__EXPORT int zRandPoisson(zRandMT *mt, double lambda);
/*! \brief a random number yielded from geometric distribution. */
__EXPORT int zRandGeo(zRandMT *mt, double p);

/*! \} */

__END_DECLS

#endif /* __ZM_RAND_H__ */
