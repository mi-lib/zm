/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_graph.h
 * \brief graph class.
 * \author Zhidao
 */

#ifndef __ZM_GRAPH_H__
#define __ZM_GRAPH_H__

#include <zm/zm_misc.h>

__BEGIN_DECLS

/* ********************************************************** */
/*! \defgroup graph.
 * \{ *//* ************************************************** */

/* advance type declaration */
struct _zGraphNode;

/* ********************************************************** */
/*! \brief graph arc class.
 *//* ******************************************************* */
typedef struct _zGraphArc{
  struct _zGraphNode *node; /*!< \brief node to connect */
  double cost; /*!< \brief cost to move */
} zGraphArc;

/* ********************************************************** */
/*! \brief graph arc set class.
 *//* ******************************************************* */
zListClass(zGraphArcList, zGraphArcListCell, zGraphArc);

/* ********************************************************** */
/*! \brief graph node class.
 *//* ******************************************************* */
typedef struct _zGraphNode{
  double val;        /*!< \brief value of node */
  double hval;       /*!< \brief heuristic estimation value of node */
  void *data;        /*!< \brief data assigned to node */
  zGraphArcList arc; /*!< \brief set of arcs to adjacent nodes */
  struct _zGraphNode *to; /*!< \brief node providing the minimum cost */
} zGraphNode;

/*! \brief initialize a graph node unidirectionally
 *
 * zGraphNodeInit() initializes a graph node \a node.
 */
__EXPORT void zGraphNodeInit(zGraphNode *node);

/*! \brief connect a graph node to another.
 *
 * zGraphNodeConnect() connects a graph node \a from to another node
 * \a to with a cost \a cost.
 * It is a unilateral arc and the inverse arc should be independently
 * defined if necessary.
 * \notes zGraphNodeBiconnect() is convenient for the bidirectional
 * connection of two nodes.
 * \retval false if it fails to allocate internal workspace.
 * \retval true otherwise.
 */
__EXPORT bool zGraphNodeConnect(zGraphNode *from, zGraphNode *to, double cost);

/*! \brief connect two graph nodes bidirectionally.
 *
 * zGraphNodeBiconnect() connects two graph nodes \a n1 and \a n2
 * bidirectionally with a cost \a cost.
 * \notes If only a unidirectional connection is desired, use
 * zGraphNodeConnect().
 * \retval false if it fails to allocate internal workspace.
 * \retval true otherwise.
 */
__EXPORT bool zGraphNodeBiconnect(zGraphNode *n1, zGraphNode *n2, double cost);

/* ********************************************************** */
/*! \brief graph node list class.
 *//* ******************************************************* */
zListClass(zGraphList, zGraphCell, zGraphNode);

/* ********************************************************** */
/*! \brief graph class.
 *//* ******************************************************* */
typedef struct{
  zGraphList list;
  /* methods */
  void *(* dup)(void*);
  bool (* equal)(void*,void*);
  void (* fwrite)(FILE*,void*);
  void (* destroy)(void*);
  double (* h)(void*,void*,void*);
} zGraph;

/*! \brief initialize a graph.
 *
 * zGraphInit() initializes a graph \a graph.
 */
__EXPORT void zGraphInit(zGraph *graph);

/*! \brief destroy a graph.
 *
 * zGraphDestroy() destroys a graph \a graph by freeing all
 * internal working memory.
 */
__EXPORT void zGraphDestroy(zGraph *graph);

/*! \brief add a node to a graph.
 *
 * zGraphAddNode() adds a node to a graph \a graph.
 * \a id is the identifier of the newly added node, while \a data
 * is data to be assigned to the node.
 * \retval false if it fails to find the node corresponding to \a id.
 * \retval false if it fails to allocate internal workspace.
 * \retval true otherwise.
 */
__EXPORT bool zGraphAddNode(zGraph *graph, void *data);

/*! \brief find a node from a graph.
 *
 * zGraphFindNode() finds a node from a graph \a graph.
 * \a id is the identifier of the node to be found.
 * \retval the pointer to the found node.
 * \retval the null poiter if it fails to find the node corresponding
 * to \a id.
 */
__EXPORT zGraphNode *zGraphFindNode(zGraph *graph, void *data);

/*! \brief connect a graph node to another.
 *
 * zGraphConnect() connects a graph node with the identifier \a from
 * to another node with the identifier \a to with a cost \a cost.
 * It is a unilateral arc and the inverse arc should be independently
 * defined if necessary.
 * \notes zGraphBiconnect() is convenient for the bidirectional
 * connection of two nodes.
 * \retval false if it fails to allocate internal workspace.
 * \retval true otherwise.
 */
__EXPORT bool zGraphConnect(zGraph *graph, void *from, void *to, double cost);

/*! \brief connect two graph nodes bidirectionally.
 *
 * zGraphBiconnect() connects a graph node with the identifier \a n1
 * and another node with the identifier \a n2 bidirectionally
 * with a cost \a cost.
 * \notes If only a unidirectional connection is desired, use
 * zGraphConnect().
 * \retval false if it fails to allocate internal workspace.
 * \retval true otherwise.
 */
__EXPORT bool zGraphBiconnect(zGraph *graph, void *n1, void *n2, double cost);

/*! \brief output information of a graph.
 *
 * zGraphFWrite() outputs information of \a graph to a stream
 * pointed by \a fp.
 */
__EXPORT void zGraphFWrite(FILE *fp, zGraph *graph);

/* ********************************************************** */
/*! \brief graph node list class.
 *//* ******************************************************* */
zListClass(zGraphNodeList, zGraphNodeListCell, zGraphNode*);

/*! \brief insert a node to a graph.
 *
 * zGraphNodeListAdd() inserts a pointer to a node \a node at the
 * head of a list \a list.
 * \retval true if succeeding to allocate a new list cell
 * \retval false if failing to allocate a new list cell
 */
__EXPORT bool zGraphNodeListAdd(zGraphNodeList *list, zGraphNode *node);

/*! \} */

__END_DECLS

#include <zm/zm_graph_search.h> /* shortest-path-finding methods */

#endif /* __ZM_GRAPH_H__ */
