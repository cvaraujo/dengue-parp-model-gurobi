//
// Created by carlos on 27/05/19.
//

#ifndef DPARP_GRAPH_H
#define DPARP_GRAPH_H

#include "../headers/Arc.h"
#include "../headers/Include.h"

using namespace std;

class Graph {
  int n, m, b;

public:
  vector<vector<Arc *>> arcs;
  vector<vector<int>> nodesPerBlock;
  vector<vector<Arc *>> arcsPerBlock;
  vector<pair<int, vector<int>>> nodes;

  Graph(string instance);

  int getN() const;

  int getM() const;

  int getB() const;

  int getRoot() const;

  void showGraph();
};

#endif
