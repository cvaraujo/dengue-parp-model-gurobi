#include <sys/stat.h>
#include <string>
#include "headers/Graph.h"
#include "headers/Model.h"
#include <chrono>

int main(int argc, const char *argv[])
{
  stringstream convTime(argv[4]), convInsec(argv[5]), convGraph(argv[6]),
      convInstanceType(argv[7]), convWarmStart(argv[8]), convFracCut(argv[9]);

  float time, insecticide;
  int graph_adaptation, instance_type, warm_start;
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

  if (strcmp(argv[3], "e") == 0) {
    model->initModelExp(time, insecticide, warm_start);
    model->solveExp("3600", frac_cut);
  } else {
    model->initModelCompact(time, insecticide);
    model->solveCompact("3600");
  }

  if (model->check_solution(time, insecticide))
    cout << "[!!!] Right result!" << endl;
  else cout << "WRONG!!!!!" << endl;
  
  model->writeSolution(argv[2]);
  return 0;
}
