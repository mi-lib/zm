/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_gen - linear equation: generalized inverse matrix.
 */

#ifndef __ZM_LE_GEN_H__
#define __ZM_LE_GEN_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/*! \brief workspace for generalized linear equation solvers. */
typedef struct{
  /*! \cond */
  /* for regular linear equation solver */
  zMat m_regular;   /* an internal matrix to be inverted. */
  zVec v_scale;     /* a scaling vector for balancing */
  zIndex index_le;  /* an index vector for pivotting */
  zVec b_copy;      /* copy of b vector */
  zVec v_ispace;    /* internal-space vector */
  /* for referenced norm-minimization */
  zVec v_ref;
  /* for MP-inverse */
  zVec v_mp_ispace; /* full-rank internal-space vector */
  /* for LU/LQ factorization */
  zMat l_facto;     /* left hand matrix */
  zMat r_facto;     /* right hand matrix */
  zIndex index_lu;  /* an index vector for LU decomposition */
  /*! \endcond */
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
 * These functions solve a generalized linear equation \a \a ans = \a b, where the coefficient
 * matrix \a a is not necessarily regular.
 *
 * zLESolveNormMinDST() and zLESolveNormMin() solve a redundant linear equation, namely, an equation
 * where the number of variables is larger than the number of equations, as to minimize the squared
 * norm of the solution.
 * The coefficient matrix \a must be row-full-rank.
 * \a w is a weighting vector, namely diagonal components of a weighting matrix on each component
 * of \a ans. If the null pointer is assigned for \a w, a uniform weight is applied.
 *
 * zLESolveErrorMinDST() and zLESolveErrorMin() solve an overconstrained linear equation, namely,
 * an equation where the number of equations is larger than the number of variables, as to minimize
 * the residual of \a a \a ans and \a b.
 * The coefficient matrix \a a must be column-full-rank.
 * \a w is a weighting vector, namely, diagonal components of a weighting matrix on the residual of
 * each equation. 
 * of \a ans. If the null pointer is assigned for \a w, a uniform weight is applied.
 *
 * zLESolveRefNormMinDST() and zLESolveRefNormMin() solve a redundant linear equation as to minimize
 * the squared error between the answer \a ans and a reference vector \a ref.
 * \a a must be row-full-rank as well as zLESolveNormMinDST() and zLESolveNormMin().
 * \a w is a weighting vector on the error.
 *
 * zLESolveNormMinDST(), zLESolveErrorMinDST(), and zLESolveRefNormMinDST() use a workspace \a workspace,
 * which should be initialized with zLEWorkspaceAlloc() beforehand and freed with zLEWorkspaceFree
 * afterward, for internal matrix-vector computations.
 * zLESolveNormMin(), zLESolveErrorMin(), and zLESolveRefNormMin() internally allocate workspace.
 * \return
 * All these functions return a pointer \a ans.
 * \sa
 * zLEWorkspaceAlloc, zLEWorkspaceFree
 */
__ZM_EXPORT zVec zLESolveNormMinDST(zMat a, zVec b, zVec w, zVec ans, zLEWorkspace *workspace);
__ZM_EXPORT zVec zLESolveNormMin(const zMat a, const zVec b, const zVec w, zVec ans);
__ZM_EXPORT zVec zLESolveErrorMinDST(zMat a, zVec b, zVec w, zVec ans, zLEWorkspace *workspace);
__ZM_EXPORT zVec zLESolveErrorMin(const zMat a, const zVec b, const zVec w, zVec ans);
__ZM_EXPORT zVec zLESolveRefNormMinDST(zMat a, zVec b, zVec w, zVec ref, zVec ans, zLEWorkspace *workspace);
__ZM_EXPORT zVec zLESolveRefNormMin(const zMat a, const zVec b, const zVec w, const zVec ref, zVec ans);

/*! \brief destructive version of generalized linear equation solver based on Moore-Penrose's inverse matrix.
 *
 * zLESolveMP_LQ_DST() solves a generalized linear equation a ans = b with internal workspace \a workspace,
 * where the coefficient matrix a has to be decomposed into L and Q matrices, and the vector b has to be
 * cloned in advance.
 *
 * zLESolveMP_LU_DST() solves a generalized linear equation a ans = b with internal workspace \a workspace,
 * where the coefficient matrix a has to be decomposed into L and U matrices, and the vector b has to be
 * cloned in advance.
 *
 * zLESolveMPNullDST() solves a generalized linear equation a ans = b and simultaneously finds the null-space
 * projector with internal workspace \a workspace, where the coefficient matrix a has to be decomposed into
 * L and Q matrices, and the vector b has to be cloned in advance.
 *
 * For these functions, \a wn is a weighting vector, namely, diagonal components of a weighting matrix on
 * the norm of the answer vector \a ans. \a we is a weighting vector, namely, diagonal components of a
 * weighting matrix on the residual of \a a \a ans and \a b.
 * \a rank is the rank of the coefficient matrix, which is the same with the column size of L matrix and
 * the row size of the Q/U matrix.
 * The answer is stored where \a ans points.
 * \a mn for zLESolveMPNullDST() is a pointer where the null-space projector is stored.
 * \return
 * These functions return a pointer \a ans.
 * \sa
 * zMatDecompLQ, zMatDecompLU
 */
__ZM_EXPORT zVec zLESolveMP_LQ_DST(zLEWorkspace *workspace, const zVec we, const zVec wn, int rank, zVec ans);
__ZM_EXPORT zVec zLESolveMP_LU_DST(zLEWorkspace *workspace, const zVec we, const zVec wn, int rank, zVec ans);
__ZM_EXPORT zVec zLESolveMPNullDST(zLEWorkspace *workspace, const zVec we, const zVec wn, int rank, zVec ans, zMat mn);

/*! \brief generalized linear equation solver based on Moore-Penrose's inverse matrix.
 *
 * zLESolveMP function family solves a generalized linear equation \a a \a ans = \a b, where the number
 * of equations and the number of variables are not necessarily equal, and the coefficient matrix is not
 * necessarily full-rank.
 * \a wn is a weighting vector, namely, diagonal components of a weighting matrix on the norm of the
 * answer vector \a ans.
 * \a we is a weighting vector, namely, diagonal components of a weighting matrix on the residual
 * of \a a \a ans and \a b.
 * They find the weighted Moore-Penrose's inverse (MP-inverse) solution \a ans, which minimizes error
 * between \a a \a ans and \b weighted by \a we with minimum weighted norm by \a wn.
 *
 * Three methods to find the MP-inverse solution are available.
 * zLESolveMP_LQ() uses LQ decomposition.
 * zLESolveMP_LU() uses LU decomposition.
 * zLESolveMP_SVD() uses the singular value decomposition, which does not accept weights.
 * zLESolveMP() is an alias to zLESolveMP_LQ().
 *
 * zLESolveMPNull finds the MP-inverse solution \a ans, and the null-space projector \a wn simultaneously.
 * \a wn is a square matrix that maps arbitrary vector to the zero vector.
 *
 * zLESolveMPAux() finds the MP-inverse solution plus a bias in the null-space, which is mapped from
 * a vector \a aux. Namely,
 *  \a ans = \a a # \a b + ( I - \a a # \a a ) \a aux,
 * where \a a # is the MP-inverse matrix of \a a, and I is the identity matrix.
 * \return
 * All these functions return a pointer \a ans.
 */
__ZM_EXPORT zVec zLESolveMP_LQ(const zMat a, const zVec b, const zVec wn, const zVec we, zVec ans);
__ZM_EXPORT zVec zLESolveMP_LU(const zMat a, const zVec b, const zVec wn, const zVec we, zVec ans);
__ZM_EXPORT zVec zLESolveMP_SVD(const zMat a, const zVec b, zVec ans);
__ZM_EXPORT zVec zLESolveMPNull(const zMat a, const zVec b, const zVec wn, const zVec we, zVec ans, zMat mn);
__ZM_EXPORT zVec zLESolveMPAux(const zMat a, const zVec b, const zVec wn, const zVec we, zVec ans, const zVec aux);

/* alias */
#define zLESolveMP zLESolveMP_LQ

/*! \brief singularity-robust generalized linear equation solver.
 *
 * These functions solve a generalized linear equation \a \a ans = \a b, where the coefficient
 * matrix \a a is not necessarily regular, with the singularity-robust inverse matrix proposed
 * by Y. Nakamura and H. Hanafusa (1986):
 *  Y. Nakamura and H. Hanafusa, Inverse Kinematic Solutions with Singularity Robustness for
 *  Robot Manipulator Control, Journal of Dynamic Systems, Measurement and Control, 108:163-171, 1986.
 * This is a simplified version of Tikhonov's regularization.
 *
 * zLESolveSRDST() and zLESolveSR() find the SR-inverse solution \a ans, which minimizes the sum
 * of the error between \a a \a ans and \a b weighted by \a we and the weighted norm by \a wn, namely,
 *  (\a a \a ans - \a b)^T diag \a we (\a a \a ans - \a b) + \a ans^T diag \a wn \a ans.
 * \a wn is a weighting vector, namely, diagonal components of a weighting matrix on the norm of the
 * answer vector \a ans.
 * \a we is a weighting vector, namely, diagonal components of a weighting matrix on the residual
 * of \a a \a ans and \a b.
 * zLESolveSRBiasDST() and zLESolveSRBias() find the SR-inverse solution \a ans, where the weight
 * on the norm of \a ans is biased by a scalar value \a bias. Namely, it minimizes
 *  (\a a \a ans - \a b)^T diag \a we (\a a \a ans - \a b) + \a ans^T ( diag \a wn + \a bias I ) \a ans,
 * where I is the identity matrix.
 *
 * zLESolveSRAuxDST() and zLESolveSRAux() find the SR-inverse solution plus a bias in the null-space,
 * which is mapped from a vector \a aux. Namely,
 *  \a ans = \a a * \a b + ( I - \a a * \a a ) \a aux,
 * where \a a * is the SR-inverse matrix of \a a.
 *
 * zLESolveSRRefDST() and zLESolveSRRef() find a solution \a ans that minimizes the sum of the error
 * between \a a \a ans and \a b weighted yb \a we and the error between \a ans and a reference vector
 * \a ref weighted by \a wn.
 *
 * zLESolveSRDST(), zLESolveSRBiasDST(), zLESolveSRAuxDST(), and zLESolveSRRefDST() use a workspace
 * \a workspace, which should be initialized with zLEWorkspaceAlloc() beforehand and freed with
 * zLEWorkspaceFree afterward, for internal matrix-vector computations.
 * zLESolveSR(), zLESolveSRBias(), zLESolveSRAux(), and zLESolveSRRef() internally allocate workspace.
 * \return
 * All these functions return a pointer \a ans.
 */
__ZM_EXPORT zVec zLESolveSRBiasDST(zMat a, zVec b, zVec wn, zVec we, double bias, zVec ans, zLEWorkspace *workspace);
__ZM_EXPORT zVec zLESolveSRBias(const zMat a, const zVec b, const zVec wn, const zVec we, double bias, zVec ans);
#define zLESolveSRDST(a,b,wn,we,ans,workspace) zLESolveSRBiasDST( a, b, wn, we, 0, ans, workspace )
#define zLESolveSR( a, b, wn, we, ans )        zLESolveSRBias( a, b, wn, we, 0, ans )
__ZM_EXPORT zVec zLESolveSRAuxDST(zMat a, zVec b, zVec wn, zVec we, zVec ans, zVec aux, zLEWorkspace *Workspace, zVec bb);
__ZM_EXPORT zVec zLESolveSRAux(const zMat a, const zVec b, const zVec wn, const zVec we, zVec ans, const zVec aux);
__ZM_EXPORT zVec zLESolveSRRefDST(zMat a, zVec b, zVec wn, zVec we, zVec ref, zVec ans, zLEWorkspace *workspace);
__ZM_EXPORT zVec zLESolveSRRef(const zMat a, const zVec b, const zVec wn, const zVec we, const zVec ref, zVec ans);

__END_DECLS

#endif /* __ZM_LE_GEN_H__ */
