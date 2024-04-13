//
// Created by carlos on 06/07/21.
//

#ifndef DPARP_STOCHASTIC_MODEL_H
#define DPARP_STOCHASTIC_MODEL_H

#include "Include.h"
#include "Graph.h"
#include <gurobi_c++.h>
#include <vector>
using namespace std;

class StochasticModel
{
  Graph *graph;
  GRBEnv env = GRBEnv();
  GRBModel model = GRBModel(env);
  vector<vector<vector<GRBVar>>> x, y, t;
  vector<vector<GRBVar>> z;
  int num_lazy_cuts, num_frac_cuts;
  float default_vel, spraying_vel, insecticide_ml_min, alpha = 1.0;

public:
  void objectiveFunction();

  StochasticModel(Graph *graph, float default_vel, float spraying_vel, float insecticide_ml_min);

  void createVariables();

  void initModelExp(float maxTime, float maxInsecticide, bool warm_start);

  void initModelCompact(float maxTime, float maxInsecticide);

  void zValue();

  void artificialNodes();

  void flowConservation();

  void maxAttending();

  void attendingPath();

  void timeConstraint(float maxTime);

  void compactTimeConstraint(float maxTime);

  void inseticideConstraint(float maxInseticide);

  float timeArc(float distance, float speed);

  float timeBlock(int block, float speed);

  float inseticideBlock(int block, float perMeter);

  float profitBlock(int block);

  void solveCompact(string timeLimit);

  void solveExp(string timeLimit, bool frac_cut);

  void writeSolution(string result);

  bool check_solution(float max_time, float max_insecticide);

  void WarmStart(float maxTime, float maxInsecticide);
};

#endif // DPARP_MODEL_H
