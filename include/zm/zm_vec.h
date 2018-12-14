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
/* CLASS: zVec
 * double precision floating point value vector class
 * ********************************************************** */

typedef struct{
  int size;
  double *elem;
} zVecStruct;
typedef zVecStruct * zVec;

#define zVecSizeNC(v)           (v)->size
#define zVecSize(v)            ( (v) ? zVecSizeNC(v) : 0 )
#define zVecSetSize(v,s)       ( zVecSizeNC(v) = (s) )
#define zVecSizeIsEqual(v1,v2) ( zVecSize(v1) == zVecSize(v2) )

/*! \brief convert vector to double array. */
#define zVecBuf(v) (v)->elem

/*! \brief get a vector element */
#define zVecElem(v,n)      zVecBuf(v)[n]
/*! \brief set a vector element */
#define zVecSetElem(v,n,e) ( zVecElem(v,n) = (e) )

/*! \brief abstract and set vector element.
 *
 * zVecSetElemList() sets all the components of \a v
 * according to the value list given by the rest of
 * arguments.
 *
 * zVecSetElemVList() is equivalent to zVecSetElemList()
 * except that it accepts va_list \a args instead of
 * variable arguments. It does not call va_end, so that
 * user functions have to call va_end themselves.
 * Note that \a args after being called is undefined.
 * \return
 * zVecSetElemList() and zVecSetElemVList() return
 * a pointer \a v.
 */
__EXPORT zVec zVecSetElemVList(zVec v, va_list args);
__EXPORT zVec zVecSetElemList(zVec v, ...);

/*! \brief create, destroy, cleanup and copy vector.
 *
 * zVecAlloc() allocates a new vector with \a size elements.
 *
 * zVecCreateList() creates a new vector from the value list
 * given by arguments. \a size is the number of elements.
 *
 * zVecFree() frees the memory allocated for a vector \a v.
 *
 * zVecFreeAO() frees multiple vectors lasted as ... at once.
 * \a n is the number of vectors to be freed.
 *
 * zVecClear() zeros all the elements of a vector \a v.
 *
 * zVecTouchup() zeros all elements less than zTOL.
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
 * zVecFree() and zVecFreeAO() return no values.
 *
 * zVecClear() and zVecTouchup() return a pointer \a v.
 *
 * zVecCopyNC() returns a pointer \a dest.
 *
 * zVecCopy() returns a pointer \a dest, or the null pointer
 * if the sizes of \a src and \a dest do not coincide.
 *
 * zVecCopyArray() returns a pointer \a v.
 * \notes
 * Since zVecCopyNC()' does not check the size consistency,
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
__EXPORT zVec zVecAlloc(int s);
__EXPORT zVec zVecCreateList(int size, ...);
__EXPORT void zVecFree(zVec v);
__EXPORT void zVecFreeAO(int, ...);
__EXPORT zVec zVecClear(zVec v);
__EXPORT zVec zVecTouchup(zVec v);
__EXPORT zVec zVecCopyNC(zVec src, zVec dest);
__EXPORT zVec zVecCopy(zVec src, zVec dest);
__EXPORT zVec zVecCopyArray(double array[], int s, zVec v);
__EXPORT zVec zVecClone(zVec src);
__EXPORT zVec zVecCloneArray(double array[], int s);

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
__EXPORT zVec zVecGetNC(zVec src, int pos, zVec dest);
__EXPORT zVec zVecGet(zVec src, int pos, zVec dest);
__EXPORT zVec zVecPutNC(zVec dest, int pos, zVec src);
__EXPORT zVec zVecPut(zVec dest, int pos, zVec src);

/*! \brief create a uniform vector, a linear space vector and a random vector.
 *
 * zVecSetAll() sets all the elements of a vector \a v for \a val.
 *
 * zVecLinSpace() divides the range \a from - \a to into
 * the values (the number of values is the same with the size
 * of \a v), and set them as elements of the vector \a v.
 * Each space between the elements are the same.
 *
 * zVecRand() sets all the elements of \a v randomly within
 * the range from \a min to \a max.
 *
 * zVecShift() shifts \a v by a scalar constant \a shift, namely,
 * \a shift is added to all the elements of \a v.
 * \return
 * zVecSetAll(), zVecLinSpace(), zVecRand() and zVecShift()
 * return a pointer \a v.
 */
__EXPORT zVec zVecSetAll(zVec v, double val);
__EXPORT zVec zVecLinSpace(zVec v, double from, double to);
__EXPORT zVec zVecRandUniform(zVec v, double min, double max);
__EXPORT zVec zVecRand(zVec v, zVec min, zVec max);
__EXPORT zVec zVecShift(zVec v, double shift);

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
__EXPORT zVec zVecSwapNC(zVec v, int i1, int i2);
__EXPORT zVec zVecSwap(zVec v, int i1, int i2);

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
__EXPORT void zVecSort(zVec v, zIndex idx);

/*! \brief maximum, minimum, average and variance of vector elements.
 *
 * 'zVecMax()' and 'zVecMin()' gets the maximum and minimum
 * component of all components of the vector 'v', respectively.
 * 'zVecAbsMax()' and 'zVecAbsMin()' gets the component
 * of 'v' whose absolute value is the maximum and minimum,
 * respectively.
 * When 'im' is not a null pointer, the index which gives
 * the maximum/minimum is stored where pointed by it.
 *
 * 'zVecSum()', 'zVecAve()' and 'zVecVar()' calculates
 * the summation, the avarage and the variance of all components
 * of 'v'.
 * [RETURN VALUE]
 * 'zVecMax()', 'zVecMin()', 'zVecAbsMax()', 'zVecAbsMin()',
 * 'zVecSum()', 'zVecAve()' and 'zVecVar()' return the results.
 */
#define zVecMax(v,i)    zDataMax( zVecBuf(v), zVecSizeNC(v), i )
#define zVecMin(v,i)    zDataMin( zVecBuf(v), zVecSizeNC(v), i )
#define zVecAbsMax(v,i) zDataAbsMax( zVecBuf(v), zVecSizeNC(v), i )
#define zVecAbsMin(v,i) zDataAbsMin( zVecBuf(v), zVecSizeNC(v), i )
#define zVecSum(v)      zDataSum( zVecBuf(v), zVecSizeNC(v) )
#define zVecAve(v)      zDataAve( zVecBuf(v), zVecSizeNC(v) )
#define zVecVar(v)      zDataVar( zVecBuf(v), zVecSizeNC(v) )

/*! \brief compare two vectors.
 *
 * zVecIsEqual() sees if the given two vector \a v1 and \a v2
 * are equal to each other.
 * \return
 * zVecIsEqual() returns the true value if \a v1 equals to
 * \a v2, or the false value otherwise.
 */
__EXPORT bool zVecIsEqual(zVec v1, zVec v2);

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
__EXPORT bool zVecIsTol(zVec v, double tol);
#define zVecIsTiny(v) zVecIsTol( v, zTOL )

__EXPORT bool zVecIsNan(zVec v);

/* METHOD:
 * zVecAddNC, zVecSubNC, zVecRevNC, zVecMulNC, zVecDivNC,
 * zVecAmpNC, zVecDemNC, zVecCatNC, zVecAddNCDRC, zVecSubNCDRC,
 * zVecRevNCDRC, zVecMulNCDRC, zVecDivNCDRC, zVecAmpNCDRC,
 * zVecDemNCDRC, zVecCatNCDRC, zVecAdd, zVecSub, zVecRev,
 * zVecMul, zVecDiv, zVecAmp, zVecDem, zVecCat, zVecAddDRC,
 * zVecSubDRC, zVecRevDRC, zVecMulDRC, zVecDivDRC, zVecAmpDRC,
 * zVecDemDRC, zVecCatDRC, zVecCats, zVecLS
 * - basic arithmetics for vector.
 *
 * zVecAddNC() and zVecAdd() add the two vectors,
 * 'v1' and 'v2', and put the result into 'v'.
 *
 * zVecSubNC() and zVecSub() subtract 'v2' from
 * 'v1', and put the result into 'v'.
 *
 * zVecRevNC() and zVecRev() reverse 'v1', and put
 * the result into 'v'.
 *
 * zVecMulNC() and zVecMul() multiply 'v1' by a
 * scalar value 'k', and put the result into 'v'.
 *
 * zVecDivNC() and zVecDiv() divide 'v1' by 'k',
 * and put the result into 'v'.
 *
 * zVecAmpNC() and zVecAmp() amplify each component
 * of 'v1' by the corresponding component of a vector
 * 'amp', and put the result into 'v'.
 *
 * zVecDemNC() and zVecDem() demamgnify each component
 * of 'v1' by the corresponding component of a vector
 * 'dem', and put the result into 'v'.
 *
 * zVecCatNC() and zVecCat() concatenate 'v1' by adding multiplied
 * 'v2' by 'k', and put the result into 'v'.
 *
 * zVecAddNCDRC() and zVecAddDRC() directly add 'v2' to 'v1'.
 *
 * zVecSubNCDRC() and zVecSubDRC() directly subtract 'v2' from 'v1'.
 *
 * zVecRevNCDRC() and zVecRevDRC() directly reverse 'v'.
 *
 * zVecMulNCDRC() and zVecMulDRC() directly multiply 'v' by 'k'.
 *
 * zVecDivNCDRC() and zVecDivDRC() directly divide 'v' by 'k'.
 *
 * zVecAmpNCDRC() and zVecAmpDRC() directly amplify 'v' by 'amp'.
 *
 * zVecDemNCDRC() and zVecDemDRC() directly demagnify 'v' by 'dem'.
 *
 * zVecCatNCDRC() and zVecCatDRC() directly concatenate 'v1' by
 * adding multiplied 'v2' by 'k'.
 *
 * zVecCats() concatenates 'n' vectors directly to 'v'.
 * Arguments follow 'n' as:
 *   'k1', 'v1', 'k2', 'v2', ...
 * where 'kx's are scalar values and 'vx' are vectors.
 * Then, the resultant 'v' will be:
 *   'v' + 'k1'*'v1' + 'k2'*'v2' + ...
 * zVecLS() computes linear sum of 'n' vectors.
 * Arguments follow 'n' as:
 *   'k1', 'v1', 'k2', 'v2', ...
 * The resultant 'v' will be:
 *   'k1'*'v1' + 'k2'*'v2' + ...
 * \return
 * Each of all these functions returns a pointer to
 * the result.
 * \notes
 * The type of NC functions do not check the size
 * consistency. If it is not urgent and you are not
 * hasty, you should not use them.
 */
__EXPORT zVec zVecAddNC(zVec v1, zVec v2, zVec v);
__EXPORT zVec zVecSubNC(zVec v1, zVec v2, zVec v);
__EXPORT zVec zVecRevNC(zVec v1, zVec v);
__EXPORT zVec zVecMulNC(zVec v1, double k, zVec v);
__EXPORT zVec zVecDivNC(zVec v1, double k, zVec v);
__EXPORT zVec zVecAmpNC(zVec v1, zVec amp, zVec v);
__EXPORT zVec zVecDemNC(zVec v1, zVec dem, zVec v);
__EXPORT zVec zVecCatNC(zVec v1, double k, zVec v2, zVec v);

#define zVecAddNCDRC(v1,v2)  zVecAddNC( v1, v2, v1 )
#define zVecSubNCDRC(v1,v2)  zVecSubNC( v1, v2, v1 )
#define zVecRevNCDRC(v)      zVecRevNC( v, v )
#define zVecMulNCDRC(v,k)    zVecMulNC( v, k, v )
#define zVecDivNCDRC(v,k)    zVecDivNC( v, k, v )
#define zVecAmpNCDRC(v,a)    zVecAmpNC( v, a, v )
#define zVecDemNCDRC(v,d)    zVecDemNC( v, d, v )
#define zVecCatNCDRC(v,k,v2) zVecCatNC( v, k, v2, v )

__EXPORT zVec zVecAdd(zVec v1, zVec v2, zVec v);
__EXPORT zVec zVecSub(zVec v1, zVec v2, zVec v);
__EXPORT zVec zVecRev(zVec v1, zVec v);
__EXPORT zVec zVecMul(zVec v1, double k, zVec v);
__EXPORT zVec zVecDiv(zVec v1, double k, zVec v);
__EXPORT zVec zVecAmp(zVec v1, zVec amp, zVec v);
__EXPORT zVec zVecDem(zVec v1, zVec dem, zVec v);
__EXPORT zVec zVecCat(zVec v1, double k, zVec v2, zVec v);

#define zVecAddDRC(v1,v2)   zVecAdd( v1, v2, v1 )
#define zVecSubDRC(v1,v2)   zVecSub( v1, v2, v1 )
#define zVecRevDRC(v)       zVecRev( v, v )
#define zVecMulDRC(v,k)     zVecMul( v, k, v )
#define zVecDivDRC(v,k)     zVecDiv( v, k, v )
#define zVecAmpDRC(v,a)     zVecAmp( v, a, v )
#define zVecDemDRC(v,d)     zVecDem( v, d, v )
#define zVecCatDRC(v1,k,v2) zVecCat( v1, k, v2, v1 )

__EXPORT zVec zVecCats(zVec v, int n, ...);
__EXPORT zVec zVecLS(zVec v, int n, ...);

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
__EXPORT double zVecInnerProdNC(zVec v1, zVec v2);
__EXPORT double zVecInnerProd(zVec v1, zVec v2);

/* METHOD:
 * zVecSqrNorm, zVecNorm, zVecNormalize, zVecNormalizeDRC,
 * zVecSqrDist, zVecDist
 * - normalize vector.
 *
 * 'zVecSqrNorm()' calculates the squared norm of the
 * vector 'v'. And 'zVecNorm()' calculates the norm
 * of 'v'.
 *
 * 'zVecNormalize()' normalizes the vector 'src',
 * dividing by the norm of itself, and put the result
 * into 'dest'.
 * 'zVecNormalizeDRC()' directly normalizes the
 * vector 'v'.
 *
 * 'zVecSqrDist()' calculates the squared distance
 * between 'v1' and 'v2', and 'zVecDist()' calculates
 * the distance between the two.
 * [RETURN VALUE]
 * 'zVecSqrNorm()', 'zVecNorm()', 'zVecSqrDist()' and
 * 'zVecDist()' return the value calculated.
 *
 * Each of 'zVecNormalize()' and 'zVecNormalizeDRC()'
 * returns the pointer to the result.
 */
__EXPORT double zVecSqrNorm(zVec v);
#define zVecNorm(v)         sqrt( zVecSqrNorm(v) )
__EXPORT double zVecWSqrNormNC(zVec v, zVec w);
#define zVecWNormNC(v,w)    sqrt( zVecWSqrNormNC(v,w) )
__EXPORT double zVecWSqrNorm(zVec v, zVec w);
#define zVecWNorm(v,w)      sqrt( zVecWSqrNorm(v,w) )
__EXPORT double zVecInfNorm(zVec v);
__EXPORT zVec zVecNormalize(zVec src, zVec dest);
#define zVecNormalizeDRC(v) zVecNormalize(v,v)
#define zVecSqrDist(v1,v2)  zRawVecSqrDist(zVecBuf(v1),zVecBuf(v2),zVecSizeNC(v1))
#define zVecDist(v1,v2)     sqrt( zVecSqrDist( v1, v2 ) )

/* METHOD:
 * zVecReadFile, zVecFRead, zVecRead,
 * zVecFWrite, zVecWrite, zVecDataFWrite, zVecDataWrite
 * - input/output vector.
 *
 * 'zVecFRead()' reads a sequence of double floating-point
 * values from the current position of the file 'fp',
 * and create a new vector.
 * The format is as follows:
 *  n ( x1 x2 x3 ... xn )
 * where 'n' is the size of vector.
 * 'zVecRead()' reads a sequence of double values according
 * to the above same format simply from the standard input.
 *
 * 'zVecReadFile()' reads a vector from file 'filename' or
 * 'filename'.zv.
 *
 * 'zVecFWrite()' writes the components of the given vector
 * 'v' to the current position of the file 'fp' in the following
 * format.
 *  n ( x1 x2 x3 ... xn )
 * 'zVecWrite()' writes the components of 'v' in the same
 * format simply to the standard output.
 *
 * 'zVecDataFWrite()' writes the components of the given vector
 * 'v' to the current position of the file 'fp' in the following
 * format.
 *  x1 x2 x3 ... xn
 * 'zVecDataWrite()' writes the components of 'v' in the
 * same format simply to the standard output.
 * [RETURN VALUE]
 * Each of 'zVecReadFile()', 'zVecFRead()' and
 * 'zVecRead()' returns a pointer to the newly created vector.
 *
 * 'zVecFWrite()', 'zVecWrite()', 'zVecDataFWrite()'
 * and 'zVecDataWrite()' return no values.
 */
#define ZVECTOR_SUFFIX "zv"
__EXPORT zVec zVecReadFile(char filename[]);
__EXPORT zVec zVecFRead(FILE *fp);
#define zVecRead()       zVecFRead( stdin )
__EXPORT void zVecFWrite(FILE *fp, zVec v);
#define zVecWrite(v)     zVecFWrite( stdout, v )
__EXPORT void zVecDataFWrite(FILE *fp, zVec v);
#define zVecDataWrite(v) zVecDataFWrite( stdout, v )

__END_DECLS

#include <zm/zm_vec_array.h> /* vector array */
#include <zm/zm_vec_list.h>  /* vector list */
#include <zm/zm_vec_tree.h>  /* vector tree */
#include <zm/zm_vec_ring.h>  /* vector ring */

#endif /* __ZM_VEC_H__ */
