/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_gen - linear equation: generalized inverse matrix.
 */

#ifndef __ZM_LE_GEN_H__
#define __ZM_LE_GEN_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief workspace for generalized linear equation solvers.
 */
typedef struct{
  zMat m; /*!< internal matrix to be inverted. */
  zMat l; /*!< left hand matrix for decomposition (LU/LQ) */
  zMat r; /*!< right hand matrix for decomposition (LU/LQ) */
  zVec b; /*!< copy of b vector */
  zVec c;
  zVec v1;
  zVec v2;
  zVec s;
  zIndex idx1;
  zIndex idx2;
} zLEWorkspace;

/*! \brief initialize workspace for generalized linear equation solver. */
__ZM_EXPORT void zLEWorkspaceInit(zLEWorkspace *workspace);
/*! \brief allocate workspace for generalized linear equation solvers. */
__ZM_EXPORT bool zLEWorkspaceAlloc(zLEWorkspace *workspace, const zVec b, int size);
/*! \brief clone workspace for generalized linear equation solvers. */
__ZM_EXPORT bool zLEWorkspaceClone(zLEWorkspace *src, zLEWorkspace *cln);
/*! \brief free workspace for generalized linear equation solvers. */
__ZM_EXPORT void zLEWorkspaceFree(zLEWorkspace *workspace);

/*! \brief generalized linear equation solver.
 *
 * zLESolve function family solves generalized linear equations,
 * in which the coefficient matrix is not necessarily square.
 * Suppose 'a x = b' is given, where 'a' is a 'r'x'c' matrix, 'b' is
 * a 'c'x1 vector and 'x', or 'ans', is a 'r'x1 vector.
 * The answer is put into \a ans.
 *
 * zLESolveNormMinDST() and zLESolveNormMin() solve redundant
 * equations, in which \a r < \a c, so as to minimize the norm
 * of \a ans. \a a must be row full rank.
 * \a w is a weighting vector, or diagonal components of a
 * weighting matrix for each component of \a ans. When the null
 * pointer is given instead, a uniform weight is applied.
 * In zLESolveNormMinDST(), \a m and \a v are temporary
 * work space matrix and vector, where \a m is \a r x \a r,
 * and \a v is \a r x1.
 * \a index is a temporary index( \a r x1), which should be
 * ordered before being given to the function by zIndexOrder(),
 * for example.
 *
 * zLESolveErrorMinDST() and zLESolveErrorMin() solve inferior
 * equations, in which \a r > \a c, so as to minimize the norm
 * of the error \a a x - \a b. \a a must be column full rank.
 * \a w is a weighting vector, or diagonal components of a
 * weighting matrix for the norm of the error \a a x - \a b.
 * When the null pointer is given instead, a uniform weight is
 * applied.
 * In zLESolveErrorMinDST(), \a m and \a v are temporary work
 * space matrices and a vector, where \a m is \a c x \a c, and
 * \a v is \a c x1. \a index is a temporary index( \a c x1),
 * which should be ordered before being given to the function
 * by zIndexOrder(), for example.
 *
 * zLESolveRefMinDST() and zLESolveRefMin() solve redundant
 * equations almost in the same style with zLESolveNormMinDST()
 * and zLESolveNormMin(). Different from norm-minimizing solvers,
 * these two functions accept a reference of the answer and
 * minimize the weighted norm of the error between \a ref and
 * \a ans. In zLESolveRefMinDST(), \a m, \a v1 and \a v2 are
 * temporary work space matrix and vectors, where \a m is
 * \a r x \a r, \a v1 is \a r x1 and \a v2 is \a r x1.
 *
 * zLESolveMP() solves the equation based on Moore-Penrose's
 * inverse (MP-inverse) matrix, where LQ decomposition is used
 * internally. One difference from the original definition, it
 * can weigh on both the residual error and the answer norm;
 * \a wn is for answer norm and \a we is for residual error.
 * zLESolveMPLU() also solves the equation based on pseudoinverse
 * matrix, where LU decomposition is used internally.
 *
 * zLESolveMPSVD() is another version of an equation solver
 * by pseudoinverse matrix, while it is based on the singular
 * value decomposition along with the original definition of
 * pseudoinverse. It cannot accept weightings.
 *
 * Probably, there is no situation where zLESolveMPLU() and
 * zLESolveMPSVD() are preferred.
 *
 * zLESolveMPAux() solves the equation based on MP-inverse,
 * biasing a vector \a aux in the null space of \a a. Namely:
 *  \a ans = \a a # \a b + ( 1 - \a a # \a a ) \a aux, where
 * \a a # is MP-inverse of \a a and 1 is the identity matrix.
 *
 * zLESolveSRDST() and zLESolveSR() solve the linear equation
 * using singularity robust inverse (or SR-inverse) matrix,
 * proposed by Y. Nakamura and H. Hanafusa (1991). \a wn and
 * \a we are weighting vectors on the norm and residual error,
 * respectively.
 * SR-inverse technique is one of the variations of Tikhonov
 * reguralization of ill-posed matrix.
 *
 * zLESolveSRAux() solves the equation based on SR-inverse,
 * biasing a vector \a aux in the null space of \a a. Namely:
 *  \a ans = \a a * \a b + ( 1 - \a a * \a a) \a aux, where
 * \a a * is SR-inverse of \a a and 1 is the identity matrix.
 *
 * zLESolveRSRDST() and zLESolveRSR() solve the linear equation
 * using referred singularity robust inverse matrix, proposed
 * by T. Sugihara (2004). In addition to \a wn and \a we, it
 * requires \a ref as a referential answer vector. The weighted
 * error between \a ref and \a ans instead of the norm of \a ans
 * is minimized.
 *
 * In zLESolveSRDST() and zLESolveRSRDST(), \a m and \a v are
 * temporary work space matrices and a vector, where \a m is
 * \a c x \a c, and \a v is \a c x1.
 * \a index is a temporary index( \a c x1), which should be
 * ordered before being given to the function by zIndexOrder(),
 * for example.
 *
 * For all DST functions, \a s is used for column-balancing.
 * See zBalancingDST() or zLESolveGaussDST().
 * \return
 * All these functions return a pointer \a ans.
 * \sa
 * zBalancingDST, zLESolveGaussDST
 */
__ZM_EXPORT zVec zLESolveNormMinDST(zMat a, zVec b, zVec w, zVec ans, zLEWorkspace *workspace);
__ZM_EXPORT zVec zLESolveNormMin(const zMat a, const zVec b, const zVec w, zVec ans);
__ZM_EXPORT zVec zLESolveErrorMinDST(zMat a, zVec b, zVec w, zVec ans, zLEWorkspace *workspace);
__ZM_EXPORT zVec zLESolveErrorMin(const zMat a, const zVec b, const zVec w, zVec ans);
__ZM_EXPORT zVec zLESolveRefMinDST(zMat a, zVec b, zVec w, zVec ref, zVec ans, zLEWorkspace *workspace);
__ZM_EXPORT zVec zLESolveRefMin(const zMat a, const zVec b, const zVec w, const zVec ref, zVec ans);

__ZM_EXPORT zVec zLESolveMP(const zMat a, const zVec b, const zVec wn, const zVec we, zVec ans);
__ZM_EXPORT zVec zLESolveMPLU(const zMat a, const zVec b, const zVec wn, const zVec we, zVec ans);
__ZM_EXPORT zVec zLESolveMPSVD(const zMat a, const zVec b, zVec ans);
__ZM_EXPORT zVec zLESolveMPNull(const zMat a, const zVec b, const zVec wn, const zVec we, zVec ans, zMat mn);
__ZM_EXPORT zVec zLESolveMPAux(const zMat a, const zVec b, const zVec wn, const zVec we, zVec ans, const zVec aux);

__ZM_EXPORT zVec zLESolveSRBiasDST(zMat a, zVec b, zVec wn, zVec we, double bias, zVec ans, zLEWorkspace *workspace);
#define zLESolveSRDST(a,b,wn,we,ans,workspace) zLESolveSRBiasDST( a, b, wn, we, 0, ans, workspace )
__ZM_EXPORT zVec zLESolveSRBias(const zMat a, const zVec b, const zVec wn, const zVec we, double bias, zVec ans);
#define zLESolveSR( a, b, wn, we, ans ) zLESolveSRBias( a, b, wn, we, 0, ans )
__ZM_EXPORT zVec zLESolveSRAuxDST(zMat a, zVec b, zVec wn, zVec we, zVec ans, zVec aux, zLEWorkspace *Workspace, zVec bb);
__ZM_EXPORT zVec zLESolveSRAux(const zMat a, const zVec b, const zVec wn, const zVec we, zVec ans, const zVec aux);
__ZM_EXPORT zVec zLESolveRSRDST(zMat a, zVec b, zVec wn, zVec we, zVec ref, zVec ans, zLEWorkspace *workspace);
__ZM_EXPORT zVec zLESolveRSR(const zMat a, const zVec b, const zVec wn, const zVec we, const zVec ref, zVec ans);

__END_DECLS

#endif /* __ZM_LE_GEN_H__ */
