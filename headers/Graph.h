//
// Created by carlos on 27/05/19.
//

#ifndef DPARP_GRAPH_H
#define DPARP_GRAPH_H

#include "../headers/Arc.h"
#include "../headers/Include.h"
#include "../headers/Scenario.h"

using namespace std;

class Graph
{
  int n, m, b, s, n_nodes;
  float default_vel = 40, spraying_vel = 10, insecticide_ml_min = 70;
  double m_time;
  graph_t boost_graph;

public:
  vector<vector<Arc *>> arcs;
  vector<set<int>> nodesPerBlock;
  vector<vector<Arc *>> arcsPerBlock;
  vector<pair<int, set<int>>> nodes;
  vector<int> cases_per_block;
  vector<Scenario> scenarios;

  Graph(string instance, int graph_adapt, int instance_type);

  Graph(string instance, string scenarios, int graph_adapt);

  void load_scenarios_instance(string instance);

  float timeArc(float distance, float speed);

  float timeBlock(int block, float speed);

  float inseticideBlock(int block, float perMeter);

  int getN() const;

  int getM() const;

  int getB() const;

  int getS() const;

  int getRoot() const;

  double getMtime();

  void showGraph();

  void showScenarios();

  bool exist_arc(int i, int j);

  void fillMissingArcs();

  void fillCompleteDigraph();

  void load_instance(string instance, int graph_adapt);

  void load_stochastic_instance(string instance, int graph_adapt);

  void loadSBRPInstance(string instance, int graph_adapt);

  void mtzAdapt();

  void init_boost_graph();

  double shortest_path(int i, int j, vector<pair<int, int>> &arcs);
};

#endif
