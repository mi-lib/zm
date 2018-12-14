/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_graph_search.h
 * \brief shortest-path-finding methods for graph class.
 * \author Zhidao
 */

#ifndef __ZM_GRAPH_SEARCH_H__
#define __ZM_GRAPH_SEARCH_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* ********************************************************** */
/*! \defgroup graph.
 * \{ *//* ************************************************** */

/*! \brief find the shortest path of a graph by A* algorithm.
 *
 * zGraphAStar() finds the shortest path of a graph \a graph
 * from \a start to \a goal, where they are the node identifiers.
 * The resulted path is stored in \a path in the order from the
 * start to the goal.
 * \retval the minimum cost to get to the goal from the start.
 * \notes a heuristic function should be assigned in advance as
 * \a graph->h = h, for instance.
 */
__EXPORT double zGraphAStar(zGraph *graph, void *start, void *goal, void *util, zGraphNodeList *path);

/*! \brief find the shortest path of a graph by Dijkstra's method.
 *
 * zGraphDijkstra() finds the shortest path of a graph \a graph
 * from \a start to \a goal, where they are the node identifiers.
 * The resulted path is stored in \a path in the order from the
 * start to the goal.
 * \retval the minimum cost to get to the goal from the start.
 */
__EXPORT double zGraphDijkstra(zGraph *graph, void *start, void *goal, zGraphNodeList *path);

/*! \} */

__END_DECLS

#endif /* __ZM_GRAPH_SEARCH_H__ */
