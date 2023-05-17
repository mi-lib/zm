/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_ga - optimization tools: genetic algorithm.
 */

#ifndef __ZM_OPT_GA_H__
#define __ZM_OPT_GA_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zOptGAChromosome
 * ********************************************************** */

typedef struct{
  double fitness;
  zVec gene;
} zOptGAChromosome;

/*! \brief initialize a chromosome.
 */
__ZM_EXPORT zOptGAChromosome *zOptGAChromosomeInit(zOptGAChromosome *c);

/*! \brief allocate a chromosome.
 *
 * zOptGAChromosomeAlloc() allocates a chromosome \a c with a \a size
 * gene = double-precision floating-point vector.
 * \return
 * zOptGAChromosomeAlloc() returns a pointer to the newly allocated
 * chromosome.
 * \sa
 * zOptGAChromosomeFree
 */
__ZM_EXPORT zOptGAChromosome *zOptGAChromosomeAlloc(zOptGAChromosome *c, size_t size);

/*! \brief free a chromosome.
 *
 * zOptGAChromosomeFree() frees gene of a chromosome \a c.
 */
__ZM_EXPORT void zOptGAChromosomeFree(zOptGAChromosome *c);

/*! \brief copy a chromosome.
 *
 * zOptGAChromosomeCopy() copies a chromosome \a src to another
 * \a dest.
 * \return
 * zOptGAChromosomeCopy() returns a pointer to \a dest.
 */
__ZM_EXPORT zOptGAChromosome *zOptGAChromosomeCopy(zOptGAChromosome *src, zOptGAChromosome *dest);

/*! \brief randomly generate a chromosome.
 *
 * zOptGAChromosomeRand() generates a chromosome \a c randomly
 * within the range between the minimum pattern \a min and the
 * maximum \a max.
 * \return
 * zOptGAChromosomeRand() returns a pointer \a c.
 */
__ZM_EXPORT zOptGAChromosome *zOptGAChromosomeRand(zOptGAChromosome *c, zVec min, zVec max);

/*! \brief crossover two chromosomes.
 *
 * zOptGAChromosomeXover() blends \a c1 and \a c2 in accordance with
 * \a c2 + ratio x ( \a c1 - \a c2 ), where ratio is f1 / ( f1 + f2 )
 * and f1 and f2 are fitnesses of \a c1 and \c2, respectively.
 * \note
 * \a rate and \a ratio should be from 0 to 1.
 */
__ZM_EXPORT zOptGAChromosome *zOptGAChromosomeXover(zOptGAChromosome *c1, zOptGAChromosome *c2, zOptGAChromosome *c);

/* ********************************************************** */
/* CLASS: zOptGA
 * ********************************************************** */

typedef struct{
  int population;
  zOptGAChromosome *individual;
  zVec min;
  zVec max;
  double rate_survival;
  double rate_mutation;
  double (* f)(zVec,void*);
  double eval;
} zOptGA;

__ZM_EXPORT void zOptGAInit(zOptGA *ga);

/*! \brief create a genetic population.
 *
 * zOptGACreate() creates a new genetic population \a ga.
 * \a num is the number of population.
 * \a min is the minimum pattern of gene.
 * \a max is the maximum pattern of gene.
 * \a fitness is an evaluation function to be minimized.
 * \a rate_survival is the survival rate.
 * \a rate_mutation is the mutation rate.
 * \return
 * zOptGACreate() returns the true value if succeed, or the false
 * value otherwise.
 */
__ZM_EXPORT bool zOptGACreate(zOptGA *ga, double (* f)(zVec,void*), void *util, zVec min, zVec max, int population, double rate_survival, double rate_mutation);

/*! \brief create a genetic population with default parameters. */
__ZM_EXPORT bool zOptGACreateDefault(zOptGA *ga, double (* f)(zVec,void*), void *util, zVec min, zVec max);

/*! \brief destroy a genetic population.
 *
 * zOptGADestroy() destroys a genetic population \a ga.
 */
__ZM_EXPORT void zOptGADestroy(zOptGA *ga);

/*! \brief randomly generate a genetic population.
 *
 * zOptGARand() randomly generate a genetic population \a ga.
 * \a util is a utility pointer in computing fitnesses of each
 * individual.
 * \return
 * zOptGAReproduce() returns the fitness of the best individual
 * in the last generation.
 */
__ZM_EXPORT double zOptGARand(zOptGA *ga, void *util);

/*! \brief reproduce a genetic population.
 *
 * zOptGAReproduce() reproduces the next generation of the genetic
 * population \a ga. \a util is a utility pointer in computing
 * fitnesses of each individual.
 * \return
 * zOptGAReproduce() returns the fitness of the best individual
 * in the last generation.
 */
__ZM_EXPORT double zOptGAReproduce(zOptGA *ga, void *util);

/*! \brief solve an optimization problem by genetic algorithm.
 *
 * zOptGASolve() tries to solve an optimization problem by genetic
 * algorithm.
 * \a ga is an instance of the optimizer, which has to be created
 * by zOptGACreate(). The answer is stored to \a ans. \a util is
 * a utility pointer to be used in computing fitnesses of each
 * individual. \a generation is the number of generation up to
 * which the genetic population is reproduced.
 * \return
 * zOptGASolve() returns the fitness of the best individual produced
 * by the genetic algorithm.
 */
__ZM_EXPORT double zOptGASolve(zOptGA *ga, zVec ans, void *util, int generation);

#define ZOPT_GA_DEFAULT_GENERATION 1000
#define ZOPT_GA_DEFAULT_POPULATION  100
#define ZOPT_GA_RATE_SURVIVAL         0.3
#define ZOPT_GA_RATE_MUTATION         0.05

/*! \brief solve an optimization problem by genetic algorithm. */
__ZM_EXPORT int zOptSolveGA(double (* f)(zVec,void*), void *util, zVec min, zVec max, int iter, double tol, int population, double rate_survival, double rate_mutation, zVec ans, double *eval);

/*! \brief solve an optimization problem by genetic algorithm. */
__ZM_EXPORT int zOptSolveGADefault(double (* f)(zVec,void*), void *util, zVec min, zVec max, int iter, double tol, zVec ans, double *eval);

/*! \brief print a genetic population out to a file.
 */
__ZM_EXPORT void zOptGAFPrint(FILE *fp, zOptGA *ga);

__END_DECLS

#endif /* __ZM_OPT_GA_H__ */
