/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_opt_ga - optimization tools: genetic algorithm.
 */

#include <zm/zm_opt.h>

/* ********************************************************** */
/* CLASS: zOptGAChromosome
 * chromosome class, supporting
 * 1. bit-chromosome
 * 2. integer-chromosome
 * 3. double-chromosome
 * ********************************************************** */

static void _zOptGAChromosomeSwap(zOptGAChromosome *c1, zOptGAChromosome *C2);

/* zOptGAChromosomeInit
 * - initialize a chromosome.
 */
zOptGAChromosome *zOptGAChromosomeInit(zOptGAChromosome *c)
{
  c->fitness = 0;
  c->gene = NULL;
  return c;
}

/* zOptGAChromosomeAlloc
 * - allocate a new chromosome.
 */
zOptGAChromosome *zOptGAChromosomeAlloc(zOptGAChromosome *c, size_t size)
{
  if( ( c->gene = zVecAlloc( size ) ) == NULL ){
    ZALLOCERROR();
    return NULL;
  }
  c->fitness = 0;
  return c;
}

/* zOptGAChromosomeFree
 * - free a chromosome.
 */
void zOptGAChromosomeFree(zOptGAChromosome *c)
{
  zVecFree( c->gene );
  c->fitness = 0;
}

/* zOptGAChromosomeCopy
 * - copy a chromosome.
 */
zOptGAChromosome *zOptGAChromosomeCopy(zOptGAChromosome *src, zOptGAChromosome *dest)
{
  zVecCopy( src->gene, dest->gene );
  dest->fitness = src->fitness;
  return dest;
}

/* zOptGAChromosomeRand
 * - randomly generate a chromosome.
 */
zOptGAChromosome *zOptGAChromosomeRand(zOptGAChromosome *c, zVec min, zVec max)
{
  zVecRand( c->gene, min, max );
  return c;
}

/* (static)
 * _zOptGAChromosomeSwap
 * - swap two chromosomes. This is utilized for quicksort of individuals.
 */
void _zOptGAChromosomeSwap(zOptGAChromosome *c1, zOptGAChromosome *c2)
{
  zSwap( double, c1->fitness, c2->fitness );
  zSwap( zVec, c1->gene, c2->gene );
}

/* zOptGAChromosomeXover
 * - crossover two chromosomes.
 */
zOptGAChromosome *zOptGAChromosomeXover(zOptGAChromosome *c1, zOptGAChromosome *c2, zOptGAChromosome *c)
{
  register int i;
  double ratio;

  ratio = c1->fitness / ( c1->fitness + c2->fitness );
  for( i=0; i<zVecSizeNC(c->gene); i++ )
    zVecElem(c->gene,i) = zVecElem(c2->gene,i) + ratio * ( zVecElem(c1->gene,i) - zVecElem(c2->gene,i) );
  return c;
}

/* ********************************************************** */
/* CLASS: zOptGA
 * ********************************************************** */

static zOptGAChromosome *_zOptGASelect(zOptGA *ga);
static void _zOptGAQuickSort(zOptGA *ga, int head, int tail);
static void _zOptGASort(zOptGA *ga, void *util);

/* zOptGAInit
 * - initialize a genetic population.
 */
void zOptGAInit(zOptGA *ga)
{
  ga->population = 0;
  ga->individual = NULL;
  ga->min = ga->max = NULL;
}

/* zOptGAPopulationCreate
 * - create a genetic population.
 */
bool zOptGACreate(zOptGA *ga, size_t size, int population, zVec min, zVec max, double (* fitness)(zOptGAChromosome*,void*), double rate_survive, double rate_mutate, void *util)
{
  register int i;

  zOptGAInit( ga );
  if( ( ga->individual = zAlloc( zOptGAChromosome, population ) ) == NULL ){
    ZALLOCERROR();
    return false;
  }
  for( i=0; i<population; i++ )
    if( zOptGAChromosomeAlloc( &ga->individual[i], size ) == NULL ) goto ERR;

  ga->population = population;
  ga->min = zVecClone( min );
  ga->max = zVecClone( max );
  if( ga->min == NULL || ga->max == NULL ) goto ERR;

  ga->rate_survive = rate_survive;
  ga->rate_mutate = rate_mutate;
  ga->fitness_fp = fitness;
  zOptGARand( ga, util );
  return true;

 ERR:
  ZALLOCERROR();
  zOptGADestroy( ga );
  return false;
}

/* zOptGADestroy
 * - destroy a genetic population.
 */
void zOptGADestroy(zOptGA *ga)
{
  register int i;

  for( i=0; i<ga->population; i++ )
    zOptGAChromosomeFree( &ga->individual[i] );
  zFree( ga->individual );
  zVecFree( ga->min );
  zVecFree( ga->max );
  zOptGAInit( ga );
}

/* (static)
 * _zOptGAQuickSort
 * - quick sort of population based on fitness.
 */
void _zOptGAQuickSort(zOptGA *ga, int head, int tail)
{
  register int i, j;
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

/* (static)
 * _zOptGASort
 * - compute fitnesses of each individual and sort them in accordance with it.
 */
void _zOptGASort(zOptGA *ga, void *util)
{
  register int i;
  double fitness_min, fitness_sum;

  for( i=0; i<ga->population; i++ )
    ga->individual[i].fitness = ga->fitness_fp( &ga->individual[i], util );
  _zOptGAQuickSort( ga, 0, ga->population-1 );
  fitness_min = ga->individual[ga->population-1].fitness;
  fitness_sum = 0;
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

/* zOptGARand
 * - randomly generate a genetic population.
 */
double zOptGARand(zOptGA *ga, void *util)
{
  register int i;

  for( i=0; i<ga->population; i++ )
    zOptGAChromosomeRand( &ga->individual[i], ga->min, ga->max );
  _zOptGASort( ga, util );
  return ga->individual[0].fitness;
}

/* (static)
 * _zOptGASelect
 * - select individuals of a genetic population
 */
zOptGAChromosome *_zOptGASelect(zOptGA *ga)
{
  double p, rate = 0;
  register int i;

  p = zRandF( 0, 1 );
  for( i=0; i<ga->population; i++ )
    if( ( rate += ga->individual[i].fitness ) > p )
      return &ga->individual[i];
  return &ga->individual[ga->population-1];
}

/* zOptGAReproduce
 * - reproduce a genetic population.
 */
double zOptGAReproduce(zOptGA *ga, void *util)
{
  zOptGAChromosome *c1, *c2;
  int ns;
  register int i;

  /* survival and selection */
  ns = ga->population * ga->rate_survive;
  /* crossover */
  for( i=ns; i<ga->population; i++ ){
    c1 = _zOptGASelect( ga );
    c2 = _zOptGASelect( ga );
    zOptGAChromosomeXover( c1, c2, &ga->individual[i] );
  }
  /* mutation */
  for( i=0; i<ga->population; i++ )
    if( zRandF(0,1) < ga->rate_mutate )
      zOptGAChromosomeRand( &ga->individual[i], ga->min, ga->max );
  /* resort */
  _zOptGASort( ga, util ); /* compute and normalize fitness */
  return ga->individual[0].fitness;
}

/* zOptGASolve
 * - solve an optimization problem by genetic algorithm.
 */
double zOptGASolve(zOptGA *ga, zVec ans, void *util, int generation)
{
  register int i;

  for( i=0; i<generation; i++ )
    zOptGAReproduce( ga, util );
  zVecCopy( ga->individual[0].gene, ans );
  return ga->individual[0].fitness;
}

/* for debug */

/* zOptGAFWrite
 * - output a genetic population to a file.
 */
void zOptGAFWrite(FILE *fp, zOptGA *ga)
{
  register int i;

  for( i=0; i<ga->population; i++ ){
    fprintf( fp, "%g ", ga->individual[i].fitness );
    zVecDataFWrite( fp, ga->individual[i].gene );
  }
}
