/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_stat.h
 * \brief statistics.
 * \author Zhidao
 */

#ifndef __ZM_STAT_H__
#define __ZM_STAT_H__

#include <zm/zm_misc.h>

__BEGIN_DECLS

/* ********************************************************** */
/*! \defgroup stat statistics.
 * \{ *//* ************************************************** */

/*! \brief permutation.
 *
 * zPermut() computes \a n_P_\a i.
 * The result is returned as a double-precision floating-point value.
 */
__EXPORT double zPermut(int n, int i);

/*! \brief factorial.
 *
 * zFacto() computes \a n!.
 * The result is returned as a double-precision floating-point value.
 */
__EXPORT double zFacto(int n);

/*! \brief combination.
 *
 * zCombi() computes \a n_C_i.
 * The result is returned as a double-precision floating-point value.
 */
__EXPORT double zCombi(int n, int i);

/*! \brief series of combination.
 *
 * zCombiSeries() computes a series of combination \a n_C_i where
 * \a i is from 0 to \a n. The result is stored in the array of
 * double-precision floating-point values \a c. \a size is the size
 * of \a c.
 * If \a n is more than \a size, equal to \a size or less than 0,
 * it fails to compute the series and the null pointer is returned.
 * Otherwise, a pointer to \a c is returned.
 */
__EXPORT double *zCombiSeries(int n, size_t size, double c[]);

/* basic distribution functions */

/*! \brief constant coefficient for normal distribution.
 */
#define zND 0.398942280401433

/*! \brief normal distribution function.
 *
 * zNormalDistrib() calculates normal distribution defined
 * by mean \a mu and variance \a sigma as:
 *  1/(sqrt(2*zPI)*\a sigma)*exp( -0.5*((\a x-\a mu)/\a sigma)^2 ).
 * In this destribution, E(x) = \a mu and V(x) = \a sigma^2.
 */
__EXPORT double zNormalDistrib(double x, double mu, double sigma);

/*! \brief normal cumulative distribution function.
 *
 * zNormalCumDistrib() calculates cumulative distribution of a
 * normal distribution with mean \a mu and variance \a sigma.
 */
__EXPORT double zNormalCumDistrib(double x, double mu, double sigma);

/*! \brief Poisson's distribution function.
 *
 * zPoissonDistrib() calculates Poisson's distribution
 * defined by \a lambda as:
 *  exp(-\a lambda)*\a lambda^\a x/\a x!.
 */
__EXPORT double zPoissonDistrib(int x, double lambda);

/*! \brief binomial distribution function.
 *
 * zBinDistrib() calculates binomial distribution defined by \a n
 * and \a p as:
 *  \a n_C_\a x*\a p^\a x*(1-\a p)^(\a n-\a x).
 */
__EXPORT double zBinDistrib(int x, int n, double p);

/*! \brief chi-squared distribution.
 *
 * zChi2Distrib() calculates chi-squared distribution defined by
 * \a x and \a k, where \a x is the chi-squared value and \a k is
 * the degree.
 */
__EXPORT double zChi2Distrib(double x, int k);

/*! \brief chi-squared cumulative distribution.
 *
 * zChi2CumDistrib() calculates chi-squared cumulative distribution
 * defined by \a x and \a k, where \a x is the chi-squared value and
 * \a k is the degree.
 */
__EXPORT double zChi2CumDistrib(double x, int k);

/* basic statistics computation */

/*! \brief the maximum value in data.
 *
 * \a num is the number of data.
 * The index to the max is stored where \a im points,
 * if it is not the null pointer.
 */
__EXPORT double zDataMax(double *data, int num, int *im);

/*! \brief the minimum value in data.
 *
 * \a num is the number of data.
 * The index to the max is stored where \a im points,
 * if it is not the null pointer.
 */
__EXPORT double zDataMin(double *data, int num, int *im);

/*! \brief the maximum absolute value in data.
 *
 * \a num is the number of data.
 * The index to the max is stored where \a im points,
 * if it is not the null pointer.
 */
__EXPORT double zDataAbsMax(double *data, int num, int *im);

/*! \brief the minimum absolute value in data.
 *
 * \a num is the number of data.
 * The index to the max is stored where \a im points,
 * if it is not the null pointer.
 */
__EXPORT double zDataAbsMin(double *data, int num, int *im);

/*! \brief the summation of data.
 *
 * \a num is the number of data.
 */
__EXPORT double zDataSum(double *data, int num);

/*! \brief the average of data.
 *
 * \a num is the number of data.
 */
__EXPORT double zDataAve(double *data, int num);

/*! \brief the variance of data.
 *
 * \a num is the number of data.
 */
__EXPORT double zDataVar(double *data, int num);

/*! \brief the standard deviation of data.
 *
 * \a num is the number of data.
 */
__EXPORT double zDataSD(double *data, int num);

/*! \} */

__END_DECLS

#endif /* __ZM_STAT_H__ */
