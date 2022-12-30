#include <sys/stat.h>
#include <string>
#include "headers/Graph.h"
#include "headers/Model.h"
#include <chrono>

int main(int argc, const char *argv[])
{
  stringstream convTime(argv[4]), convInsec(argv[5]), convGraph(argv[6]),
      convInstanceType(argv[7]), convWarmStart(argv[8]);

  float time, insecticide;
  int graph_adaptation, instance_type, warm_start;

  convGraph >> graph_adaptation;
  convTime >> time;
  convInsec >> insecticide;
  convInstanceType >> instance_type;
  convWarmStart >> warm_start;

  // Create Graph
  cout << "Loading the graph!" << endl;
  Graph *g = new Graph(argv[1], graph_adaptation, instance_type);

  cout << "Creating Gurobi model!" << endl;
  if (strcmp(argv[3], "e") == 0)
  {
    Model *model = new Model(g);
    model->createVariables();
    model->initModelExp(time, insecticide, warm_start);
    model->solveExp("3600");
    bool feasible = model->check_solution();
    if (feasible)
      cout << "Right result!" << endl;
    else
      cout << "WRONG!!!!!" << endl;
    model->showSolution(argv[2]);
  }
  else
  {
    cout << "------------------" << endl;
    Model *model2 = new Model(g);
    model2->createVariables();
    model2->initModelCompact(time, insecticide);
    model2->solveCompact("3600");
    bool feasible = model2->check_solution();
    if (feasible)
      cout << "Right result!" << endl;
    else
      cout << "WRONG!!!!!" << endl;
    model2->showSolution(argv[2]);
  }
  return 0;
}
