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
__ZM_EXPORT double zNormalDistrib(double x, double mu, double sigma);

/*! \brief normal cumulative distribution function.
 *
 * zNormalCumDistrib() calculates cumulative distribution of a
 * normal distribution with mean \a mu and variance \a sigma.
 */
__ZM_EXPORT double zNormalCumDistrib(double x, double mu, double sigma);

/*! \brief Poisson's distribution function.
 *
 * zPoissonDistrib() calculates Poisson's distribution
 * defined by \a lambda as:
 *  exp(-\a lambda)*\a lambda^\a x/\a x!.
 */
__ZM_EXPORT double zPoissonDistrib(int x, double lambda);

/*! \brief binomial distribution function.
 *
 * zBinDistrib() calculates binomial distribution defined by \a n
 * and \a p as:
 *  \a n_C_\a x*\a p^\a x*(1-\a p)^(\a n-\a x).
 */
__ZM_EXPORT double zBinDistrib(int x, int n, double p);

/*! \brief chi-squared distribution.
 *
 * zChi2Distrib() calculates chi-squared distribution defined by
 * \a x and \a k, where \a x is the chi-squared value and \a k is
 * the degree.
 */
__ZM_EXPORT double zChi2Distrib(double x, int k);

/*! \brief chi-squared cumulative distribution.
 *
 * zChi2CumDistrib() calculates chi-squared cumulative distribution
 * defined by \a x and \a k, where \a x is the chi-squared value and
 * \a k is the degree.
 */
__ZM_EXPORT double zChi2CumDistrib(double x, int k);

/* basic statistics computation */

/*! \brief the maximum value in data.
 *
 * zDataMax() finds the maximum element of data \a data.
 * \a size is the size of \a data.
 * The index of the maximum element is stored where \a im points, unless it is the null pointer.
 * \retval the maximum value of \a data.
 */
__ZM_EXPORT double zDataMax(const double *data, int size, int *im);

/*! \brief the minimum value in data.
 *
 * zDataMin() finds the minimum element of data \a data.
 * \a size is the size of \a data.
 * The index of the minimum element is stored where \a im points, unless it is the null pointer.
 */
__ZM_EXPORT double zDataMin(const double *data, int size, int *im);

/*! \brief minimum and maximum elements of data.
 *
 * zDataMinMax() finds the minimum and maximum elements of data \a data.
 * \a size is the size of the array \a data.
 * The minimum and maximum values are stored where \a min and \a max point, respectively, unless
 * they are the null pointers.
 * The indices of the minimum and maximum elements are stored where \a imin and \a imax point,
 * respectively, unless they are the null pointers.
 * \return
 * zDataMinMax() does not return any values.
 */
__ZM_EXPORT void zDataMinMax(const double *data, int size, double *min, int *imin, double *max, int *imax);

/*! \brief the maximum absolute value in data.
 *
 * zDataAbsMax() finds the maximum absolute element of data \a data.
 * \a size is the size of \a data.
 * The index of the maximum element is stored where \a im points, unless it is the null pointer.
 */
__ZM_EXPORT double zDataAbsMax(const double *data, int size, int *im);

/*! \brief the minimum absolute value in data.
 *
 * zDataAbsMax() finds the minimum absolute element of data \a data.
 * \a size is the size of \a data.
 * The index of the minimum element is stored where \a im points, unless it is the null pointer.
 */
__ZM_EXPORT double zDataAbsMin(const double *data, int size, int *im);

/*! \brief the summation of data.
 *
 * \a size is the size of data.
 */
__ZM_EXPORT double zDataSum(const double *data, int size);

/*! \brief the mean of data.
 *
 * \a size is the size of data.
 */
__ZM_EXPORT double zDataMean(const double *data, int size);

/*! \brief the variance of data.
 *
 * \a size is the size of data.
 */
__ZM_EXPORT double zDataVar(const double *data, int size);

/*! \brief the standard deviation of data.
 *
 * \a size is the size of data.
 */
__ZM_EXPORT double zDataStandardDeviation(const double *data, int size);

/*! \brief check if a value is included in data.
 *
 * zDataIsIncluded() checks if a value \a val is included in
 * a data set \a data. \a size is the size of the data.
 * \a tol is the tolerance to regard two values as the same.
 * \sa zIsEqual
 */
__ZM_EXPORT bool zDataIsIncluded(const double *data, int size, double val, double tol);

/*! \} */

#include <zm/zm_stat_histogram.h>

__END_DECLS

#endif /* __ZM_STAT_H__ */
