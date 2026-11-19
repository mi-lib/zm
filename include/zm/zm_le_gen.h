/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_le_gen - linear equation: generalized inverse matrix and Moore-Penrose inverse matrix.
 */

#ifndef __ZM_LE_GEN_H__
#define __ZM_LE_GEN_H__

/* NOTE: never include this header file in user programs. */

/*! \brief workspace for generalized linear equation solvers. */
ZDEF_STRUCT( __ZM_CLASS_EXPORT, zLEWorkspace ){
  /*! \cond */
  /* for regular linear equation solver */
  zMat m_regular;   /* an internal matrix to be inverted. */
  zVec v_scale;     /* a scaling vector for balancing */
  zIndex index_le;  /* an index vector for pivotting */
  zVec b_copy;      /* copy of b vector */
  zVec v_ispace;    /* internal-space vector */
  /* for Householder transformations */
  zMat m_house;     /* a copy of the coefficient matrix (or its transpose) to be transformed */
  /* for MP-inverse */
  zVec v_mp_ispace; /* full-rank internal-space vector */
  /* for LU/LQ factorization */
  zMat l_facto;     /* left hand matrix */
  zMat r_facto;     /* right hand matrix */
  zIndex index_row; /* an index vector of rows for LU/LQ decomposition */
  zIndex index_col; /* an index vector of columns for LU decomposition with complete pivoting */
  /*! \endcond */
};

__BEGIN_DECLS

/*! \brief initialize workspace for generalized linear equation solver. */
__ZM_EXPORT void zLEWorkspaceInit(zLEWorkspace *workspace);
/*! \brief allocate workspace for generalized linear equation solvers. */
__ZM_EXPORT bool zLEWorkspaceAlloc(zLEWorkspace *workspace, int num_equation, int num_variable);
/*! \brief clone workspace for generalized linear equation solvers. */
__ZM_EXPORT bool zLEWorkspaceClone(zLEWorkspace *src, zLEWorkspace *cln);
/*! \brief free workspace for generalized linear equation solvers. */
__ZM_EXPORT void zLEWorkspaceFree(zLEWorkspace *workspace);

/*! \brief allocate workspace for LQ/LU decomposition and generalized linear equation solvers based on MP inverse. */
__ZM_EXPORT bool zLEWorkspaceAllocMP(zLEWorkspace *workspace, int num_equation, int num_variable);
/*! \brief free workspace for LQ/LU decomposition and generalized linear equation solvers based on MP inverse. */
__ZM_EXPORT void zLEWorkspaceFreeMP(zLEWorkspace *workspace);
/*! \brief resize matrices and a vector for LQ/LU decomposition and generalized linear equation solvers based on MP inverse. */
__ZM_EXPORT void zLEWorkspaceResizeMP(zLEWorkspace *workspace, int rank);

/*! \brief generalized linear equation solver.
 *
 * These functions solve a generalized linear equation \a a \a ans = \a b, where the coefficient
 * matrix \a a is not necessarily regular.
 *
 * zLESolveNormMinDST() and zLESolveNormMin() solve a redundant linear equation, namely, an equation
 * where the number of variables is larger than the number of equations, as to minimize the weighted
 * squared norm of the solution \a ans^T W^-1 \a ans, where W is a diagonal matrix of which the diagonal
 * components are given by a weighting vector \a w. If the null pointer is assigned for \a w, a uniform
 * weight is applied. The coefficient matrix \a a must be row-full-rank.
 *
 * zLESolveErrorMinDST() and zLESolveErrorMin() solve an overconstrained linear equation, namely,
 * an equation where the number of equations is larger than the number of variables, as to minimize
 * the weighted squared norm of the residual ( \a a \a ans - \a b )^T W ( \a a \a ans - \a b ), where W
 * is a diagonal matrix of which the diagonal components are given by a weighting vector \a w. If the
 * null pointer is assigned for \a w, a uniform weight is applied.
 * The coefficient matrix \a a must be column-full-rank.
 *
 * zLESolveRefNormMinDST() and zLESolveRefNormMin() solve a redundant linear equation as to minimize
 * the weighted squared error between the answer \a ans and a reference vector \a ref, namely,
 * ( \a ans - \a ref )^T W^-1 ( \a ans - \a ref ).
 * \a a must be row-full-rank as well as zLESolveNormMinDST() and zLESolveNormMin().
 *
 * The components of the weighting vectors have to be non-negative.
 *
 * They are based on Householder's method without forming the normal equations, so that the error of
 * the answer is proportional to the condition number of \a a, not to its square. zLESolveNormMinDST()
 * and zLESolveRefNormMinDST() decompose \a a W^1/2 into a lower triangular matrix and a row-orthonormal
 * matrix (LQ decomposition), and zLESolveErrorMinDST() decomposes W^1/2 \a a into a column-orthonormal
 * matrix and an upper triangular matrix (QR decomposition).
 * \a a is regarded as rank-deficient if the norm of the residual of a row (a column for
 * zLESolveErrorMinDST()) of the weighted coefficient matrix orthogonalized to the preceding ones is not
 * greater than zTOL times the maximum norm of the rows (columns) of the weighted coefficient matrix,
 * which is the same tolerance with that of zMatDecompLQ() (zMatDecompQR()). Since the rows (columns) are
 * not pivoted, the rank might be judged differently from zMatDecompLQ() for nearly rank-deficient matrices.
 * zLESolveErrorMinDST() destroys \a b. zLESolveNormMinDST() and zLESolveRefNormMinDST() do not, and
 * \a ans can be the same vector with \a b for them.
 *
 * zLESolveNormMinDST(), zLESolveErrorMinDST(), and zLESolveRefNormMinDST() use a workspace \a workspace,
 * which should be initialized with zLEWorkspaceAlloc() beforehand and freed with zLEWorkspaceFree
 * afterward, for internal matrix-vector computations, where the numbers of equations and variables
 * given to zLEWorkspaceAlloc() have to be those of \a a.
 * zLESolveNormMin(), zLESolveErrorMin(), and zLESolveRefNormMin() internally allocate workspace.
 * \return
 * These functions return a pointer \a ans if they succeed. If \a a is rank-deficient, or the workspace
 * is not large enough, they return the null pointer.
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
 * These functions solve a generalized linear equation a ans = b with internal workspace \a workspace,
 * where the weighted coefficient matrix We^1/2 a Wn^1/2 has to be decomposed in advance into the factors
 * l_facto and r_facto of \a workspace, and b (not weighted) has to be copied to b_copy of \a workspace.
 * We and Wn are diagonal matrices of which the diagonal components are given by \a we and \a wn,
 * respectively. The null pointer for \a we or \a wn means the identity matrix.
 *
 * zLESolveMP_LQ_DST() and zLESolveMPNullDST() accept the factors of the LQ decomposition by
 * zMatDecompLQReorth() with the index vector index_row of \a workspace. Those by zMatDecompLQ() are also
 * accepted, while the residual of the answer might increase in proportion to the condition number of the
 * weighted coefficient matrix due to the loss of orthogonality of the right factor.
 *
 * zLESolveMP_LU_DST() accepts the factors of the LU decomposition with complete pivoting by
 * zMatDecompLUCP() with the index vectors index_row and index_col of \a workspace. Note that the
 * decomposition by zMatDecompLU() based on partial pivoting is not accepted, since the U matrix is
 * supposed to be an upper triangular matrix with columns permuted by index_col.
 *
 * zLESolveMPNullDST() simultaneously finds the null-space projector \a mn (see zLESolveMPNull()).
 *
 * \a wn is a weighting vector, namely, diagonal components of a weighting matrix on the norm of the
 * answer vector \a ans. \a we is a weighting vector, namely, diagonal components of a weighting matrix
 * on the residual of \a a \a ans and \a b.
 * \a rank is the rank of the weighted coefficient matrix, which is the same with the column size of L
 * matrix and the row size of the Q/U matrix.
 * The answer is stored where \a ans points.
 * \return
 * These functions return a pointer \a ans if they succeed. If they fail to solve the equations of the
 * factors, they return the null pointer. zLESolveMPNullDST() also returns the null pointer if a
 * component of \a wn is not positive.
 * \sa
 * zMatDecompLQReorth, zMatDecompLQ, zMatDecompLUCP
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
 * They find the weighted Moore-Penrose's inverse (MP-inverse) solution \a ans, which minimizes
 *  ( \a a \a ans - \a b )^T We ( \a a \a ans - \a b ),
 * and then minimizes \a ans^T Wn^-1 \a ans among the minimizers, where We and Wn are diagonal matrices
 * of which the diagonal components are given by \a we and \a wn, respectively. The null pointer for
 * \a we or \a wn means the identity matrix. The components of the weighting vectors have to be
 * non-negative. An equation with zero weight is ignored, and a variable with zero weight is fixed to zero.
 * Since \a ans = Wn^1/2 y, where y is the MP-inverse solution of We^1/2 \a a Wn^1/2 y = We^1/2 \a b,
 * the weighted coefficient matrix We^1/2 \a a Wn^1/2 is decomposed. Its rank is detected in the
 * decomposition (see zMatDecompLQReorth() and zMatDecompLUCP()), and the directions regarded as
 * degenerated are discarded.
 *
 * Three methods to find the MP-inverse solution are available.
 * zLESolveMP_LQ() uses LQ decomposition.
 * zLESolveMP_LU() uses LU decomposition with complete pivoting (see zMatDecompLUCP()).
 * zLESolveMP_SVD() uses the singular value decomposition, which does not accept weights.
 * zLESolveMP() is an alias to zLESolveMP_LQ().
 *
 * zLESolveMPNull() finds the MP-inverse solution \a ans, and simultaneously the null-space projector
 * \a mn = \a a # \a a - I, where \a a # is the weighted MP-inverse matrix of \a a and I is the identity
 * matrix, namely, \a mn = Wn^1/2 ( Q^T Q - I ) Wn^-1/2 for the LQ decomposition L Q of the weighted
 * coefficient matrix. \a a \a mn = O holds if all the components of \a we are positive, and
 * \a ans - \a mn v is the MP-inverse solution biased by a vector v in the null-space (see zLESolveMPAux()).
 * The components of \a wn have to be positive for zLESolveMPNull().
 *
 * zLESolveMPAux() finds the MP-inverse solution plus a bias in the null-space, which is mapped from
 * a vector \a aux. Namely,
 *  \a ans = \a a # \a b + ( I - \a a # \a a ) \a aux.
 * \return
 * All these functions return a pointer \a ans if they succeed. If they fail to allocate internal
 * workspace or to solve the equations of the factors, they return the null pointer.
 * zLESolveMPNull() also returns the null pointer if a component of \a wn is not positive.
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
 * These functions solve a generalized linear equation \a a \a ans = \a b, where the coefficient
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

/*! \brief Moore-Penrose inverse matrix.
 *
 * zMatMPInv() calculates the Moore-Penrose inverse (MP-inverse) matrix of an arbitrary matrix \a m.
 * The result is put into \a mp. Suppose a matrix B is the MP-inverse matrix of a matrix A. The
 * following relations are satisfied.
 *   A B A = A
 *   B A B = B
 *   ( A B )^T = A B
 *   ( B A )^T = B A
 * It is computed with the LQ decomposition \a m = L Q (see zMatDecompLQReorth()) as \a mp = Q^T L#, where
 * L# is the MP-inverse matrix of the column-full-rank matrix L, which is found by forward substitution
 * if L is square, or with the QR decomposition of L based on Householder's method otherwise. The rank
 * is detected in the LQ decomposition, and the directions regarded as degenerated are discarded.
 *
 * zMatMPInvNull() calculates the MP-inverse matrix \a mp of \a m, and simultaneously the null-space
 * projector \a mn = \a mp \a m - I, where I is the identity matrix.
 *
 * zMulMPInvMatMat() multiplies the MP-inverse matrix of \a m1 by \a m2 from the left side without
 * forming the MP-inverse matrix explicitly, and puts the result into \a m. \a m2 and \a m can be the
 * same matrix.
 *
 * zMatMPInvPenrose() calculates the MP-inverse matrix \a mp of \a m based on Penrose's iterative method.
 * \return
 * zMatMPInv(), zMatMPInvNull() and zMatMPInvPenrose() return the rank of \a m. zMatMPInv() and
 * zMatMPInvNull() return -1 if the sizes of the matrices do not match or they fail to allocate
 * internal workspace.
 * zMulMPInvMatMat() returns a pointer \a m if it succeeds, or the null pointer otherwise.
 */
__ZM_EXPORT int zMatMPInv(const zMat m, zMat mp);
__ZM_EXPORT int zMatMPInvNull(const zMat m, zMat mp, zMat mn);
__ZM_EXPORT zMat zMulMPInvMatMat(const zMat m1, const zMat m2, zMat m);
__ZM_EXPORT int zMatMPInvPenrose(const zMat m, zMat mp);

__END_DECLS

#endif /* __ZM_LE_GEN_H__ */
