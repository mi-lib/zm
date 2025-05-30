/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_vec - vector class.
 */

#ifndef __ZM_VEC_H__
#define __ZM_VEC_H__

#include <zeda/zeda.h>
#include <zm/zm_complex.h>
#include <zm/zm_raw.h>
#include <zm/zm_stat.h>

__BEGIN_DECLS

/* ********************************************************** */
/*! \brief double-precision floating-point value vector class.
 *//* ******************************************************* */
zArrayClass( zVecStruct, double );
typedef zVecStruct * zVec;

/*! \brief size of a vector.
 * \retval the size of a vector if \a v is not null.
 * \retval 0 if \a v is the null pointer.
 */
#define zVecSizeNC(v)          zArraySize(v)
#define zVecSize(v)            ( (v) ? zVecSizeNC(v) : 0 )
/*! \brief set the size of a vector. */
#define zVecSetSizeNC(v,s)     ( zVecSizeNC(v) = (s) )
#define zVecSetSize(v,s)       ( (v) ? zVecSetSizeNC(v,s) : 0 )

/*! \brief check if the sizes of two vectors are equal. */
#define zVecSizeEqual(v1,v2) ( zVecSizeNC(v1) == zVecSizeNC(v2) )

/*! \brief check if the specified position of a vector is valid. */
#define zVecPosIsValid(v,n)    zArrayPosIsValid( v, n )

/*! \brief pointer to the array buffer of of double-precision floating-point values in a vector. */
#define zVecBufNC(v)       zArrayBuf(v)
#define zVecBuf(v)         ( (v) ? zArrayBuf(v) : NULL )

/*! \brief get an element of a vector without checking size. */
#define zVecElemNC(v,n)    zVecBufNC(v)[n]
/*! \brief get an element of a vector. */
#define zVecElem(v,n)      ( (v) && zVecPosIsValid(v,n) ? zVecElemNC(v,n) : 0 )
/*! \brief set an element of a vector without checking size. */
#define zVecSetElemNC(v,n,e) ( zVecElemNC(v,n) = (e) )
/*! \brief set an element of a vector. */
#define zVecSetElem(v,n,e) ( (v) && zVecPosIsValid(v,n) ? zVecSetElemNC(v,n,e) : 0 )

/*! \brief set elements of a vector for values in the argument list.
 *
 * zVecSetElemList() sets all elements of a vector \a v for the
 * values given by the argument list.
 *
 * zVecSetElemVList() is equivalent to zVecSetElemList() except
 * that it accepts va_list \a args instead of variable argument
 * list. It does not call va_end, so that user functions have to
 * call va_end themselves.
 * \note \a args after being called is undefined.
 * \return
 * zVecSetElemList() and zVecSetElemVList() return a pointer \a v.
 */
__ZM_EXPORT zVec zVecSetElemVList(zVec v, va_list args);
__ZM_EXPORT zVec zVecSetElemList(zVec v, ...);

/*! \brief create, destroy, cleanup and copy vector.
 *
 * zVecAlloc() allocates a new vector with \a size elements.
 *
 * zVecCreateList() creates a new vector from the value list
 * given by arguments. \a size is the number of elements.
 *
 * zVecFree() frees the memory allocated for a vector \a v.
 *
 * zVecFreeAtOnce() frees multiple vectors lasted as ... at once.
 * \a n is the number of vectors to be freed.
 *
 * zVecZero() sets all components of a vector \a v for zero.
 *
 * zVecTouchup() replaces all components less than \a tol of \a v
 * for zeros.
 *
 * zVecCopyNC() copies a vector \a src to the other \a dest
 * without checking the size consistency between \a src and
 * \a dest.
 *
 * zVecCopy() copies a vector \a src to the other \a dest.
 *
 * zVecCopyArray() copies \a s elements in the array \a array
 * to a vector \a v.
 *
 * zVecClone() creates a clone of a vector \a src.
 * zVecCloneArray() creates a clone of an array \a array.
 * \a s is the size of the array.
 * \return
 * zVecAlloc(), zVecCreateList(), zVecClone() and zVecCloneArray()
 * return a pointer to the newly created vector.
 *
 * zVecFree() and zVecFreeAtOnce() return no values.
 *
 * zVecZero() and zVecTouchup() return a pointer \a v.
 *
 * zVecCopyNC() returns a pointer \a dest.
 *
 * zVecCopy() returns a pointer \a dest, or the null pointer
 * if the sizes of \a src and \a dest do not coincide.
 *
 * zVecCopyArray() returns a pointer \a v.
 * \notes
 * Since zVecCopyNC() does not check the size consistency,
 * if the size of \a src and \a dest are different from each
 * other, anything might happen.
 * If it is not urgent and you are not hasty, you would better
 * use zVecCopy() for safety.
 *
 * Because of a specification of glibc, the following call
 * does not work as expected.
 *   v = zVecCreateList( 3, 1, 2, 3 );
 * It should be written as follows.
 *   v = zVecCreateList( 3, 1.0, 2.0, 3.0 );
 */
__ZM_EXPORT zVec zVecAlloc(int s);
__ZM_EXPORT zVec zVecCreateList(int size, ...);
__ZM_EXPORT void zVecFree(zVec v);
__ZM_EXPORT void zVecFreeAtOnce(int num, ...);
__ZM_EXPORT zVec zVecZero(zVec v);
__ZM_EXPORT zVec zVecTouchup(zVec v, double tol);
__ZM_EXPORT zVec zVecCopyNC(const zVec src, zVec dest);
__ZM_EXPORT zVec zVecCopy(const zVec src, zVec dest);
__ZM_EXPORT zVec zVecCopyArray(const double array[], int s, zVec v);
__ZM_EXPORT zVec zVecClone(const zVec src);
__ZM_EXPORT zVec zVecCloneArray(const double array[], int s);

/*! \brief get and put a part of a vector.
 *
 * zVecGetNC() and zVecGet() copy \a src to a part of
 * \a dest from the \a pos 'th element. The number of
 * copied elements is the size of \a dest. It is expected
 * that \a src has a larger size than \a pos + the size
 * of \a dest.
 *
 * zVecPutNC() and zVecPut() copy \a src to a part of
 * \a dest from the \a pos 'th element. The number of
 * copied elements is the size of \a src. It is expected
 * that \a dest has a larger size than \a pos + the size
 * of \a src.
 * \return
 * zVecGetNC() and zVecPutNC() always return a pointer
 * \a dest without checking the size validity between \a src
 * and \a dest, while zVecGet() and zVecPut() return the
 * null pointer if the size of \a src and \a dest are
 * inconsistent.
 * \sa
 * zRawVecGet, zRawVecPut
 */
__ZM_EXPORT zVec zVecGetNC(const zVec src, int pos, zVec dest);
__ZM_EXPORT zVec zVecGet(const zVec src, int pos, zVec dest);
__ZM_EXPORT zVec zVecPutNC(zVec dest, int pos, const zVec src);
__ZM_EXPORT zVec zVecPut(zVec dest, int pos, const zVec src);

/*! \brief create a uniform vector, a linear space vector and a random vector.
 *
 * zVecSetAll() sets all components of a vector \a v for \a val.
 *
 * zVecLinSpace() divides the range \a from - \a to into the values (the number of values is the
 * same with the size of \a v), and set them as components of the vector \a v.
 * Each space between the components are the same.
 *
 * zVecRand() sets all the components of \a v randomly within the range from \a min to \a max.
 *
 * zVecShift() shifts all components of a vector \a src by a scalar constant \a shift.
 * The result is put into \a dest.
 * zVecShiftDRC() directly shifts all components of a vector \a v by a scalar constant \a shift.
 * \return
 * zVecSetAll(), zVecLinSpace(), zVecRand(), and zVecShiftDRC() return a pointer \a v.
 * zVecShift() returns a pointer \a dest.
 */
__ZM_EXPORT zVec zVecSetAll(zVec v, double val);
__ZM_EXPORT zVec zVecLinSpace(zVec v, double from, double to);
__ZM_EXPORT zVec zVecRandUniform(zVec v, double min, double max);
__ZM_EXPORT zVec zVecRand(zVec v, zVec min, zVec max);
__ZM_EXPORT zVec zVecShift(const zVec src, double shift, zVec dest);
#define zVecShiftDRC(vec,shift) zVecShift( vec, shift, vec )

/*! \brief - swap vector elements.
 *
 * zVecSwapNC() and zVecSwap() swap \a i1 'th and \a i2 'th
 * elements of \a v', where zVecSwapNC() does not check the
 * size consistency of \a v with \a i1 and \a i2.
 * \return
 * zVecSwapNC() and zVecSwap() return a pointer \a v.
 * \sa
 * zRawVecSwap
 */
__ZM_EXPORT zVec zVecSwapNC(zVec v, int i1, int i2);
__ZM_EXPORT zVec zVecSwap(zVec v, int i1, int i2);

/*! \brief rearrange index so as to sort vector.
 *
 * zVecSort() rearranges the index vector \a idx so as to sort
 * the given vector \a v in descending order. Namely:
 *   v[idx[0]] <= v[idx[1]] <= ... <= v[idx[N]]
 * where N is the size of \a v.
 * \return
 * zVecSort() returns no value.
 * \notes
 * If zVecSort() fails to allocate working memory, or the sizes
 * of \a v and \a idx are inconsistent, it does not anything.
 */
__ZM_EXPORT void zVecSort(zVec v, zIndex idx);

/*! \brief reorder components of a vector along with a given index.
 *
 * zVecReorder() reorders components of a vector \a src along with a given index \a idx, and puts the
 * result into \a dest.
 * zVecReorderDRC() directly reorders a vector \a v along with a given index \a idx.
 * \return
 * zVecReorder() returns the pointer \a dest.
 * zVecReorderDRC() returns the pointer \a v.
 */
__ZM_EXPORT zVec zVecReorder(const zVec src, const zIndex idx, zVec dest);
__ZM_EXPORT zVec zVecReorderDRC(zVec v, const zIndex idx);

/*! \brief maximum, minimum, summation, mean and variance of vector elements.
 *
 * zVecElemMax() and zVecElemMin() find the maximum and minimum component of all components of a vector
 * \a v, respectively.
 * zVecElemAbsMax() and zVecElemAbsMin() find the component of \a v whose absolute value is the maximum
 * and minimum, respectively.
 * For those four functions, the index that gives the maximum/minimum is stored where pointed by \a im,
 * unless it is the null pointer.
 *
 * zVecElemSum(), zVecElemMean() and zVecElemVar() calculate the summation, the mean and the variance of
 * all components of \a v.
 * \return
 * zVecElemMax(), zVecElemMin(), zVecElemAbsMax(), zVecElemAbsMin(), zVecElemSum(), zVecElemMean() and
 * zVecElemVar() return the results.
 */
#define _zVecElemMax(v,im)    zDataMax( zVecBuf(v), zVecSizeNC(v), im )
#define _zVecElemMin(v,im)    zDataMin( zVecBuf(v), zVecSizeNC(v), im )
#define _zVecElemAbsMax(v,im) zDataAbsMax( zVecBuf(v), zVecSizeNC(v), im )
#define _zVecElemAbsMin(v,im) zDataAbsMin( zVecBuf(v), zVecSizeNC(v), im )
#define _zVecElemSum(v)       zDataSum( zVecBuf(v), zVecSizeNC(v) )
#define _zVecElemMean(v)      zDataMean( zVecBuf(v), zVecSizeNC(v) )
#define _zVecElemVar(v)       zDataVar( zVecBuf(v), zVecSizeNC(v) )

__ZM_EXPORT double zVecElemMax(const zVec v, int *im);
__ZM_EXPORT double zVecElemMin(const zVec v, int *im);
__ZM_EXPORT double zVecElemAbsMax(const zVec v, int *im);
__ZM_EXPORT double zVecElemAbsMin(const zVec v, int *im);
__ZM_EXPORT double zVecElemSum(const zVec v);
__ZM_EXPORT double zVecElemMean(const zVec v);
__ZM_EXPORT double zVecElemVar(const zVec v);

/*! \brief check if a value is included in a vector.
 *
 * zVecValIsIncluded() checks if a value \a val is included in a vector \a v.
 */
#define zVecValIsIncluded(v,val,tol) zDataIsIncluded( zVecBuf(v), zVecSizeNC(v), val, tol )

/*! \brief compare two vectors.
 *
 * zVecEqual() checks if the given two vectors \a v1 and \a v2 are equal to each other. \a tol is the
 * tolerance to regard two values as the same.
 *
 * zVecMatch() checks if the given two vectors \a v1 and \a v2 exactly match with each other.
 * \return
 * zVecEqual() returns the true value if \a v1 equals to \a v2, or the false value otherwise.
 * zVecMatch() returns the true value if \a v1 exactly matches with \a v2, or the false value otherwise.
 */
__ZM_EXPORT bool zVecEqual(const zVec v1, const zVec v2, double tol);
__ZM_EXPORT bool zVecMatch(const zVec v1, const zVec v2);

/*! \brief check if a vector is tiny.
 *
 * zVecIsTol() returns the true value if all the components of
 * the vector \a v are less than \a tol, or the false value
 * otherwise.
 * zVecIsTiny() is the same with zVecIsTol() except it compares
 * each component with zTOL(defined in zeda_misc.h) instead of
 * \a tol.
 * \return
 * zVecIsTol() and zVecIsTiny() return a boolean value.
 */
__ZM_EXPORT bool zVecIsTol(const zVec v, double tol);
#define zVecIsTiny(v) zVecIsTol( v, zTOL )

__ZM_EXPORT bool zVecIsNan(const zVec v);

/*! \brief basic arithmetics for vector.
 *
 * zVecAddNC() and zVecAdd() add the two vectors, \a v1 and \a v2, and put the result into \a v.
 *
 * zVecSubNC() and zVecSub() subtract \a v2 from \a v1, and put the result into \a v.
 *
 * zVecRevNC() and zVecRev() reverse \a v1, and put the result into \a v.
 *
 * zVecMulNC() and zVecMul() multiply \a v1 by a scalar value \a k, and put the result into \a v.
 *
 * zVecDivNC() and zVecDiv() divide \a v1 by \a k, and puts the result into \a v.
 *
 * zVecAmpNC() and zVecAmp() amplify each component of \a v1 by the corresponding component of
 * a vector \a amp, and puts the result into \a v.
 *
 * zVecDemNC() and zVecDem() demamgnify each component of \a v1 by the corresponding component
 * of a vector \a dem, and put the result into \a v.
 *
 * zVecCatNC() and zVecCat() concatenate \a v1 by adding multiplied \a v2 by \a k, and puts the
 * result into \a v.
 *
 * zVecAddNCDRC() and zVecAddDRC() directly add \a v2 to \a v1.
 *
 * zVecSubNCDRC() and zVecSubDRC() directly subtract \a v2 from \a v1.
 *
 * zVecRevNCDRC() and zVecRevDRC() directly reverse \a v.
 *
 * zVecMulNCDRC() and zVecMulDRC() directly multiply \a v by \a k.
 *
 * zVecDivNCDRC() and zVecDivDRC() directly divide \a v by \a k.
 *
 * zVecAmpNCDRC() and zVecAmpDRC() directly amplify \a v by \a amp.
 *
 * zVecDemNCDRC() and zVecDemDRC() directly demagnify \a v by \a dem.
 *
 * zVecCatNCDRC() and zVecCatDRC() directly concatenate \a v1 by adding multiplied \a v2 by \a k.
 *
 * zVecCats() concatenates \a n vectors directly with \a v.
 * Arguments follow \a n as:
 *   \a k1, \a v1, \a k2, \a v2, ...
 * where \a kx s are scalar values and \a vx are vectors.
 * Then, the result \a v will be:
 *   \a v + \a k1 *\a v1 + \a k2 * \a v2 + ...
 * zVecLinearSum() computes linear sum of \a n vectors.
 * Arguments follow \a n as:
 *   \a k1, \a v1, \a k2, \a v2, ...
 * The result \a v will be:
 *   \a k1 * \a v1 + \a k2 * \a v2 + ...
 * \return
 * Each of all these functions returns a pointer to the result.
 * \notes
 * The type of NC functions do not check the size consistency. If it is not urgent and you are not
 * hasty, you should not use them.
 */
__ZM_EXPORT zVec zVecAddNC(const zVec v1, const zVec v2, zVec v);
__ZM_EXPORT zVec zVecSubNC(const zVec v1, const zVec v2, zVec v);
__ZM_EXPORT zVec zVecRevNC(const zVec v1, zVec v);
__ZM_EXPORT zVec zVecMulNC(const zVec v1, double k, zVec v);
__ZM_EXPORT zVec zVecDivNC(const zVec v1, double k, zVec v);
__ZM_EXPORT zVec zVecAmpNC(const zVec v1, const zVec amp, zVec v);
__ZM_EXPORT zVec zVecDemNC(const zVec v1, const zVec dem, zVec v);
__ZM_EXPORT zVec zVecCatNC(const zVec v1, double k, const zVec v2, zVec v);

#define zVecAddNCDRC(v1,v2)  zVecAddNC( v1, v2, v1 )
#define zVecSubNCDRC(v1,v2)  zVecSubNC( v1, v2, v1 )
#define zVecRevNCDRC(v)      zVecRevNC( v, v )
#define zVecMulNCDRC(v,k)    zVecMulNC( v, k, v )
#define zVecDivNCDRC(v,k)    zVecDivNC( v, k, v )
#define zVecAmpNCDRC(v,a)    zVecAmpNC( v, a, v )
#define zVecDemNCDRC(v,d)    zVecDemNC( v, d, v )
#define zVecCatNCDRC(v,k,v2) zVecCatNC( v, k, v2, v )

__ZM_EXPORT zVec zVecAdd(const zVec v1, const zVec v2, zVec v);
__ZM_EXPORT zVec zVecSub(const zVec v1, const zVec v2, zVec v);
__ZM_EXPORT zVec zVecRev(const zVec v1, zVec v);
__ZM_EXPORT zVec zVecMul(const zVec v1, double k, zVec v);
__ZM_EXPORT zVec zVecDiv(const zVec v1, double k, zVec v);
__ZM_EXPORT zVec zVecAmp(const zVec v1, const zVec amp, zVec v);
__ZM_EXPORT zVec zVecDem(const zVec v1, const zVec dem, zVec v);
__ZM_EXPORT zVec zVecCat(const zVec v1, double k, const zVec v2, zVec v);

#define zVecAddDRC(v1,v2)   zVecAdd( v1, v2, v1 )
#define zVecSubDRC(v1,v2)   zVecSub( v1, v2, v1 )
#define zVecRevDRC(v)       zVecRev( v, v )
#define zVecMulDRC(v,k)     zVecMul( v, k, v )
#define zVecDivDRC(v,k)     zVecDiv( v, k, v )
#define zVecAmpDRC(v,a)     zVecAmp( v, a, v )
#define zVecDemDRC(v,d)     zVecDem( v, d, v )
#define zVecCatDRC(v1,k,v2) zVecCat( v1, k, v2, v1 )

__ZM_EXPORT zVec zVecCats(zVec v, int n, ...);
__ZM_EXPORT zVec zVecLinearSum(zVec v, int n, ...);

/*! \brief interior division of two vectors.
 *
 * zVecInterDiv() calculates the interior division vector of two vectors \a v1 and
 * \a v2 with a division ratio \a ratio. The result is put into \a v,
 * i.e. \a v = (1-\a ratio)* \a v1 + \a ratio * \a v2.
 *
 * zVecInterDivDRC() replaces a vector \a v with the interior division with another
 * \a v2,
 * i.e. \a v = (1-\a ratio)* \a v + \a ratio * \a v2.
 * \return
 * zVecInterDiv() and zVecInterDivDRC() return a pointer \a v.
 */
__ZM_EXPORT zVec zVecInterDiv(const zVec v1, const zVec v2, double ratio, zVec v);
__ZM_EXPORT zVec zVecInterDivDRC(zVec v, const zVec v2, double ratio);

/*! \brief midpoint of two vectors.
 *
 * zVecMid() calculates the midpoint of two vectors \a v1 and \a v2.
 * The result is put into \a v.
 *
 * i.e. \a v = ( \a v1 + \a v2 ) / 2.
 * \return
 * zVecMid() returns a pointer \a v.
 */
__ZM_EXPORT zVec zVecMid(const zVec v1, const zVec v2, zVec v);

/*! \brief scale a vector with two boundary vectors.
 *
 * zVecScale() scales a vector \a src with the minimum and maximum boundary vectors \a min and \a max,
 * respectively, and puts it into \a dest. Namely, each component of \a dest is the corresponding
 * component of \a src multiplied by the corresponding components of \a max - \a min and offset by
 * that of \a min.
 *
 * zVecScaleUniform() uniformly scales \a src with the minimum and maximum boundary values \a min
 * and \a max, and puts it into \a dest.
 * \return
 * zVecScale() and zVecScaleUniform() return a pointer \a dest.
 */
__ZM_EXPORT zVec zVecScale(const zVec src, zVec min, zVec max, zVec dest);
__ZM_EXPORT zVec zVecScaleUniform(const zVec src, double min, double max, zVec dest);

/*! \brief inner product of vector.
 *
 * zVecInnerProdNC() and zVecInnerProd() calculates the
 * inner products of the two vector \a v1 and \a v2.
 * \return
 * zVecInnerProdNC() and zVecInnerProd() return the inner
 * product.
 * \notes
 * zVecInnerProdNC() does not check the size consistency
 * between \a v1 and \a v2. If it is not urgent and you
 * are not hasty, you should use zVecInnerProd().
 */
__ZM_EXPORT double zVecInnerProdNC(const zVec v1, const zVec v2);
__ZM_EXPORT double zVecInnerProd(const zVec v1, const zVec v2);

/*! \brief normalize a vector.
 *
 * zVecSqrNorm() calculates the squared norm of the vector
 * \a v. zVecNorm() calculates the norm of \a v.
 *
 * zVecNormalize() normalizes the vector \a src, namely,
 * \a src is divided by the norm of itself. The result is
 * put into \a dest.
 * zVecNormalizeDRC() directly normalizes the vector \a v.
 *
 * zVecSqrDist() calculates the squared distance between
 * \a v1 and \a v2. zVecDist() calculates the distance
 * between the two.
 * \return
 * zVecSqrNorm(), zVecNorm(), zVecSqrDist() and zVecDist()
 * return the value calculated.
 *
 * zVecNormalize() and zVecNormalizeDRC() return the
 * pointer to the result.
 */
__ZM_EXPORT double zVecSqrNorm(const zVec v);
#define zVecNorm(v)         sqrt( zVecSqrNorm(v) )
__ZM_EXPORT double zVecWSqrNormNC(const zVec v, const zVec w);
#define zVecWNormNC(v,w)    sqrt( zVecWSqrNormNC(v,w) )
__ZM_EXPORT double zVecWSqrNorm(const zVec v, const zVec w);
#define zVecWNorm(v,w)      sqrt( zVecWSqrNorm(v,w) )
__ZM_EXPORT double zVecInfNorm(const zVec v);
__ZM_EXPORT zVec zVecNormalize(const zVec src, zVec dest);
#define zVecNormalizeDRC(v) zVecNormalize(v,v)
#define zVecSqrDist(v1,v2)  zRawVecSqrDist(zVecBuf(v1),zVecBuf(v2),zVecSizeNC(v1))
#define zVecDist(v1,v2)     sqrt( zVecSqrDist( v1, v2 ) )

/*! \brief projection and orthogonalization of a vector.
 *
 * zVecProj() projects a vector \a v onto the line directed by \a n, and put the result into \a pv;
 * \a pv is parallel to \a n.
 *
 * zVec3DOrthogonalize() orthogonalizes \a v with respect to \a n, and puts it into \a ov. Namely,
 * \a ov is orthogonal to \a n.
 * \return
 * zVec3DProj() returns a pointer \a pv.
 * zVec3DOrthogonalize() returns a pointer \a ov.
 * Both functions return the null pointer if sizes of given vectors are not equal.
 */
__ZM_EXPORT zVec zVecProj(const zVec v, const zVec n, zVec pv);
__ZM_EXPORT zVec zVecOrthogonalize(const zVec v, const zVec n, zVec ov);

/*! \brief distance from a vector to an edge.
 *
 * zVecEdgeDist() computes the distance from a vector \a v to an edge that connects two ends \a v1 and \a v2,
 * namely, the radius of a circle that centers \a v and is tangential to the edge.
 * \return
 * zVecEdgeDist() returns the computed distance from \a v to the edge \a v1 - \a v2.
 */
__ZM_EXPORT double zVecEdgeDist(const zVec v, const zVec v1, const zVec v2);

/* ********************************************************** */
/* I/O
 * ********************************************************** */

/*! \brief read a vector from a ZTK format processor. */
__ZM_EXPORT zVec zVecFromZTK(ZTK *ztk);

/*! \brief scan and print a vector.
 *
 * These functions scan/print a vector from/to a file.
 * The format is as follows:
 *  \a n ( \a x1 \a x2 \a x3 ... \a xn )
 * where \a n is the size of vector.
 *
 * zVecFScan() scans a sequence of double-precision floating-point values
 * from the current position of a file \a fp, and creates a new vector.
 * zVecScan() scans a sequence of double-precision floating-point values
 * from the standard input.
 *
 * zVecFPrint() prints the components of a vector \a v out to the current
 * position of a file \a fp.
 * zVecPrint() prints the components of \a v out to the standard output.
 * \return
 * zVecFScan() and zVecScan() return a pointer to the newly created vector.
 *
 * zVecFPrint() and zVecPrint() return no values.
 */
__ZM_EXPORT zVec zVecFScan(FILE *fp);
#define zVecScan()   zVecFScan( stdin )
__ZM_EXPORT void zVecFPrint(FILE *fp, const zVec v);
#define zVecPrint(v) zVecFPrint( stdout, v )

/*! \brief scan and print only components of a vector.
 *
 * zVecValueFScan() scans a sequence of double-precision floating-point values from the current
 * position of a file \a fp, and creates a new vector. The accepted format is only a sequence
 * of values without any parentheses and braces as
 *  \a x1 \a x2 \a x3 ... \a xn
 * zVecValueScan() scans the standard input in the same format and creates a vector.
 *
 * zVecValueFPrint() prints only components of a vector \a v out to the current position of a
 * file \a fp in the following format:
 *  \a x1 \a x2 \a x3 ... \a xn
 * zVecValuePrint() prints components of \a v out to the standard output.
 * \return
 * zVecValueFScan() and zVecValueScan() return the pointer to the newly created vector.
 *
 * zVecValueFPrint() and zVecValuePrint() return no values.
 */
__ZM_EXPORT zVec zVecValueFScan(FILE *fp);
#define zVecValueScan(v) zVecValueFScan( stdin, v )
__ZM_EXPORT void zVecValueFPrint(FILE *fp, const zVec v);
#define zVecValuePrint(v) zVecValueFPrint( stdout, v )

__END_DECLS

#include <zm/zm_vec_array.h> /* vector array */
#include <zm/zm_vec_list.h>  /* vector list */
#include <zm/zm_vec_tree.h>  /* vector tree */
#include <zm/zm_vec_ring.h>  /* vector ring */

#endif /* __ZM_VEC_H__ */
