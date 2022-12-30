//
// Created by carlos on 27/05/19.
//

#ifndef DPARP_GRAPH_H
#define DPARP_GRAPH_H

#include "../headers/Arc.h"
#include "../headers/Include.h"

using namespace std;

class Graph
{
  int n, m, b, n_nodes;
  double m_time;
  graph_t boost_graph;

public:
  vector<vector<Arc *>> arcs;
  vector<set<int>> nodesPerBlock;
  vector<vector<Arc *>> arcsPerBlock;
  vector<pair<int, set<int>>> nodes;

  Graph(string instance, int graph_adapt, int instance_type);

  double timeArc(float distance, float speed);

  double timeBlock(float speed, int block);

  float inseticideBlock(float perMeter, int block);

  int getN() const;

  int getM() const;

  int getB() const;

  int getRoot() const;

  double getMtime();

  void showGraph();

  bool exist_arc(int i, int j);

  void fillMissingArcs();

  void fillCompleteDigraph();

  void load_instance(string instance, int graph_adapt);

  void loadSBRPInstance(string instance, int graph_adapt);

  void mtzAdapt();

  void init_boost_graph();

  double shortest_path(int i, int j, vector<pair<int, int>> &arcs);
};

#endif
