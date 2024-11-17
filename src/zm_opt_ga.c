/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_ga - optimization tools: genetic algorithm.
 */

#include <zm/zm_opt.h>

/* initialize a chromosome. */
zOptGAChromosome *zOptGAChromosomeInit(zOptGAChromosome *c)
{
  c->fitness = 0;
  c->gene = NULL;
  return c;
}

/* allocate a new chromosome. */
zOptGAChromosome *zOptGAChromosomeAlloc(zOptGAChromosome *c, size_t size)
{
  if( ( c->gene = zVecAlloc( size ) ) == NULL ){
    ZALLOCERROR();
    return NULL;
  }
  c->fitness = 0;
  return c;
}

/* free a chromosome. */
void zOptGAChromosomeFree(zOptGAChromosome *c)
{
  zVecFree( c->gene );
  c->fitness = 0;
}

/* copy a chromosome to another. */
zOptGAChromosome *zOptGAChromosomeCopy(zOptGAChromosome *src, zOptGAChromosome *dest)
{
  zVecCopy( src->gene, dest->gene );
  dest->fitness = src->fitness;
  return dest;
}

/* randomly generate a chromosome. */
zOptGAChromosome *zOptGAChromosomeRand(zOptGAChromosome *c, zVec min, zVec max)
{
  zVecRand( c->gene, min, max );
  return c;
}

/* swap two chromosomes. This is utilized for quicksort of individuals. */
static void _zOptGAChromosomeSwap(zOptGAChromosome *c1, zOptGAChromosome *c2)
{
  zSwap( double, c1->fitness, c2->fitness );
  zSwap( zVec, c1->gene, c2->gene );
}

/* crossover two chromosomes. */
zOptGAChromosome *zOptGAChromosomeXover(zOptGAChromosome *c1, zOptGAChromosome *c2, zOptGAChromosome *c)
{
  int i;
  double ratio;

  ratio = c1->fitness / ( c1->fitness + c2->fitness );
  for( i=0; i<zVecSizeNC(c->gene); i++ )
    zVecElemNC(c->gene,i) = zVecElemNC(c2->gene,i) + ratio * ( zVecElemNC(c1->gene,i) - zVecElemNC(c2->gene,i) );
  return c;
}

/* ********************************************************** */
/* CLASS: zOptGA
 * ********************************************************** */

/* initialize a genetic population. */
void zOptGAInit(zOptGA *ga)
{
  ga->population = 0;
  ga->individual = NULL;
  ga->min = ga->max = NULL;
}

/* create a genetic population. */
bool zOptGACreate(zOptGA *ga, double (* f)(zVec,void*), void *util, zVec min, zVec max, int population, double rate_survival, double rate_mutation)
{
  int i;

  if( !zVecSizeIsEqual( min, max ) ){
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );
    return false;
  }
  zOptGAInit( ga );
  if( ( ga->individual = zAlloc( zOptGAChromosome, population ) ) == NULL ){
    ZALLOCERROR();
    return false;
  }
  for( i=0; i<population; i++ )
    if( zOptGAChromosomeAlloc( &ga->individual[i], zVecSizeNC(min) ) == NULL ) goto ERR;

  ga->population = population;
  ga->min = zVecClone( min );
  ga->max = zVecClone( max );
  if( ga->min == NULL || ga->max == NULL ) goto ERR;

  ga->rate_survival = rate_survival;
  ga->rate_mutation = rate_mutation;
  ga->f = f;
  ga->eval = HUGE_VAL;
  zOptGARand( ga, util );
  return true;

 ERR:
  ZALLOCERROR();
  zOptGADestroy( ga );
  return false;
}

/* create a genetic population with default parameters. */
bool zOptGACreateDefault(zOptGA *ga, double (* f)(zVec,void*), void *util, zVec min, zVec max)
{
  return zOptGACreate( ga, f, util, min, max, ZOPT_GA_DEFAULT_POPULATION, ZOPT_GA_RATE_SURVIVAL, ZOPT_GA_RATE_MUTATION );
}

/* destroy a genetic population. */
void zOptGADestroy(zOptGA *ga)
{
  int i;

  for( i=0; i<ga->population; i++ )
    zOptGAChromosomeFree( &ga->individual[i] );
  zFree( ga->individual );
  zVecFree( ga->min );
  zVecFree( ga->max );
  zOptGAInit( ga );
}

/* quick sort of population based on fitness. */
static void _zOptGAQuickSort(zOptGA *ga, int head, int tail)
{
  int i, j;
  double pivot;

  pivot = ga->individual[(head+tail)/2].fitness;
  for( i=head, j=tail; i<=j; i++, j-- ){
    while( ga->individual[i].fitness > pivot ) i++;
    while( ga->individual[j].fitness < pivot ) j--;
    if( i > j ) break;
    _zOptGAChromosomeSwap( &ga->individual[i], &ga->individual[j] );
  }
  if( j > head ) _zOptGAQuickSort( ga, head, j );
  if( i < tail ) _zOptGAQuickSort( ga, i, tail );
}

/* compute fitnesses of each individual and sort them in accordance with it. */
static void _zOptGASort(zOptGA *ga, void *util)
{
  int i;
  double fitness_min, fitness_sum;

  for( i=0; i<ga->population; i++ )
    ga->individual[i].fitness = -ga->f( ga->individual[i].gene, util );
  _zOptGAQuickSort( ga, 0, ga->population-1 );
  fitness_min = ga->individual[ga->population-1].fitness;
  fitness_sum = 0;
  ga->eval = -ga->individual[0].fitness; /* optimum value of the current population */
  for( i=0; i<ga->population; i++ ){
    ga->individual[i].fitness -= fitness_min;
    fitness_sum += ga->individual[i].fitness;
  }
  if( zIsTiny( fitness_sum ) )
    for( i=0; i<ga->population; i++ )
      ga->individual[i].fitness = 1.0;
  else
    for( i=0; i<ga->population; i++ )
      ga->individual[i].fitness /= fitness_sum;
}

/* randomly generate a genetic population. */
double zOptGARand(zOptGA *ga, void *util)
{
  int i;

  for( i=0; i<ga->population; i++ )
    zOptGAChromosomeRand( &ga->individual[i], ga->min, ga->max );
  _zOptGASort( ga, util );
  return ga->individual[0].fitness;
}

/* select individuals of a genetic population */
static zOptGAChromosome *_zOptGASelect(zOptGA *ga)
{
  double p, rate = 0;
  int i;

  p = zRandF( 0, 1 );
  for( i=0; i<ga->population; i++ )
    if( ( rate += ga->individual[i].fitness ) > p )
      return &ga->individual[i];
  return &ga->individual[ga->population-1];
}

/* reproduce a genetic population. */
double zOptGAReproduce(zOptGA *ga, void *util)
{
  zOptGAChromosome *c1, *c2;
  int ns;
  int i;

  /* survival and selection */
  ns = ga->population * ga->rate_survival;
  /* crossover */
  for( i=ns; i<ga->population; i++ ){
    c1 = _zOptGASelect( ga );
    c2 = _zOptGASelect( ga );
    zOptGAChromosomeXover( c1, c2, &ga->individual[i] );
  }
  /* mutation */
  for( i=0; i<ga->population; i++ )
    if( zRandF(0,1) < ga->rate_mutation )
      zOptGAChromosomeRand( &ga->individual[i], ga->min, ga->max );
  /* resort */
  _zOptGASort( ga, util ); /* compute and normalize fitness */
  return ga->individual[0].fitness;
}

/* solve an optimization problem by genetic algorithm. */
double zOptGASolve(zOptGA *ga, zVec ans, void *util, int generation)
{
  int i;

  if( generation == 0 ) generation = ZOPT_GA_DEFAULT_GENERATION;
  for( i=0; i<generation; i++ )
    zOptGAReproduce( ga, util );
  zVecCopy( ga->individual[0].gene, ans );
  return ga->eval;
}

/* solve an optimization problem by genetic algorithm. */
int zOptSolveGA(double (* f)(zVec,void*), void *util, zVec min, zVec max, int iter, double tol, int population, double rate_survival, double rate_mutation, zVec ans, double *eval)
{
  zOptGA opt;

  zRandInit();
  zOptGACreate( &opt, f, util, min, max, population, rate_survival, rate_mutation );
  *eval = zOptGASolve( &opt, ans, util, iter );
  zOptGADestroy( &opt );
  return iter;
}

/* solve an optimization problem by genetic algorithm with default parameters. */
int zOptSolveGADefault(double (* f)(zVec,void*), void *util, zVec min, zVec max, int iter, double tol, zVec ans, double *eval)
{
  return zOptSolveGA( f, util, min, max, iter, tol, ZOPT_GA_DEFAULT_POPULATION, ZOPT_GA_RATE_SURVIVAL, ZOPT_GA_RATE_MUTATION, ans, eval );
}

/* for debug */

/* print a genetic population out to a file. */
void zOptGAFPrint(FILE *fp, zOptGA *ga)
{
  int i;

  for( i=0; i<ga->population; i++ ){
    fprintf( fp, "%g ",-ga->individual[i].fitness );
    zVecDataFPrint( fp, ga->individual[i].gene );
  }
}
