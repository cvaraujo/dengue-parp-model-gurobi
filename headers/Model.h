//
// Created by carlos on 06/07/21.
//

#ifndef DPARP_MODEL_H
#define DPARP_MODEL_H

#include "Include.h"
#include "Graph.h"
#include <gurobi_c++.h>
#include <vector>
using namespace std;

class Model {
  Graph *graph;
  GRBEnv env = GRBEnv();
  GRBModel model = GRBModel(env);
  vector<vector<GRBVar>> x, y, t;
  // vector<GRBVar> t;

public:
  void objectiveFunction();

  Model(Graph *graph);

  void createVariables();

  void initModelExp(float maxTime, float maxInsecticide);

  void initModelCompact(float maxTime, float maxInsecticide);

  void artificialNodes();

  void flowConservation();

  void maxAttending();

  void attendingPath();

  void timeConstraint(float maxTime);

  void compactTimeConstraint(float maxTime);

  void inseticideConstraint(float maxInseticide);

  float timeArc(float distance, float speed);

  float timeBlock(float speed, int block);

  float inseticideBlock(float perMeter, int block);

  float profitBlock(int block);

  void solveCompact(string timeLimit);

  void solveExp(string timeLimit);

  void showSolution(string output);

  void writeSolution(string instance, int preprocessingTime);
};


#endif //DPARP_MODEL_H
