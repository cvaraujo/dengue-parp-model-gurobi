#include <sys/stat.h>
#include <string>
#include "headers/Graph.h"
#include "headers/Model.h"
#include <chrono>
#include "headers/StochasticModel.h"

void stochastic(char *argv[])
{
  stringstream convTime(argv[5]), convInsec(argv[6]), convGraph(argv[7]),
      convInstanceType(argv[8]), convWarmStart(argv[9]), convFracCut(argv[10]);

  float time, insecticide;
  int graph_adaptation, instance_type = 0, warm_start;
  bool frac_cut;

  convGraph >> graph_adaptation;
  convTime >> time;
  convInsec >> insecticide;
  convInstanceType >> instance_type;
  convWarmStart >> warm_start;
  convFracCut >> frac_cut;

  float default_vel = 40, spraying_vel = 10, insecticide_ml_min = 70;

  // Create Graph
  cout << "Loading the graph!" << endl;
  Graph *g = new Graph(argv[1], graph_adaptation, instance_type);

  cout << "Instantiating Gurobi model!" << endl;
  StochasticModel *model = new StochasticModel(g, default_vel, spraying_vel, insecticide_ml_min);
  model->createVariables();

  if (strcmp(argv[4], "e") == 0)
  {
    model->initModelExp(time, insecticide, warm_start);
    model->solveExp("3600", frac_cut);
  }
  else
  {
    model->initModelCompact(time, insecticide);
    model->solveCompact("3600");
  }

  if (model->check_solution(time, insecticide))
    cout << "[!!!] Right result!" << endl;
  else
    cout << "WRONG!!!!!" << endl;

  model->writeSolution(argv[3]);
}

void deterministic(const char *argv[])
{
  stringstream convTime(argv[4]), convInsec(argv[5]), convGraph(argv[6]),
      convInstanceType(argv[7]), convWarmStart(argv[8]), convFracCut(argv[9]);

  float time, insecticide;
  int graph_adaptation, instance_type = 0, warm_start;
  bool frac_cut;

  convGraph >> graph_adaptation;
  convTime >> time;
  convInsec >> insecticide;
  convInstanceType >> instance_type;
  convWarmStart >> warm_start;
  convFracCut >> frac_cut;

  float default_vel = 40, spraying_vel = 10, insecticide_ml_min = 70;

  // Create Graph
  cout << "Loading the graph!" << endl;
  Graph *g = new Graph(argv[1], graph_adaptation, instance_type);

  cout << "Instantiating Gurobi model!" << endl;
  Model *model = new Model(g, default_vel, spraying_vel, insecticide_ml_min);
  model->createVariables();

  if (strcmp(argv[3], "e") == 0)
  {
    model->initModelExp(time, insecticide, warm_start);
    model->solveExp("3600", frac_cut);
  }
  else
  {
    model->initModelCompact(time, insecticide, warm_start);
    model->solveCompact("3600");
  }

  if (model->check_solution(time, insecticide))
    cout << "[!!!] Right result!" << endl;
  else
    cout << "WRONG!!!!!" << endl;

  model->writeSolution(argv[2]);
}

int main(int argc, const char *argv[])
{
  deterministic(argv);

  return 0;
}
