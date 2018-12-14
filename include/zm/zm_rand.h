/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_rand.h
 * \brief random number generator.
 *
 * This implementation of Mersenne twister proposed by
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

/*! \brief a random value generator based on the normal distribution.
 *
 * zRandFND() returns a randomly generated value in accordance with
 * a Gaussian distribution with the mean \a mu and the variance \a sigma.
 */
__EXPORT double zRandFND(double mu, double sigma);

/* ********************************************************** */
/*! \defgroup rand random number generator with Mersenne twister.
 * \{ *//* ************************************************** */

/*! \cond */
#define Z_RAND_MT_HISTORY 623
/*! \endcond */

/* ********************************************************** */
/*! \brief Mersenne twister class.
 *//********************************************************* */
typedef struct{
  ulong  x[Z_RAND_MT_HISTORY+1]; /*!< state vector */
  int    index;                  /*!< index */
  bool   nd_sw;                  /*!< a switch for normal distribution */
  double nd_last;                /*!< a memory for normal distribution */
} zRandMT;

/*! \brief initialize Mersenne twister.
 *
 * zRandInitMT() initializes Mersenne twister \a mt
 * by seeding the current time for a new history.
 * If the null pointer is given for \a mt, it makes
 * use of the internal default instance.
 */
__EXPORT void zRandInitMT(zRandMT *mt);

/*! \brief a pseudo-random integer between \a min and \a max. */
__EXPORT int zRandMTI(zRandMT *mt, int min, int max);

/*! \brief a pseudo-random double-precision floating-point
 * value between \a min and \a max. */
__EXPORT double zRandMTF(zRandMT *mt, double min, double max);

/*! \brief a pseudo-random double-precision floating-point
 * value in the range of [0,1]. */
__EXPORT double zRandMTN(zRandMT *mt);

/*! \brief a pseudo-random value in the range of [0,1). */
__EXPORT double zRandMTNU(zRandMT *mt);

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

/*! \brief a random number yielded from normal distribution. */
__EXPORT double zRandNormal(zRandMT *mt);
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
