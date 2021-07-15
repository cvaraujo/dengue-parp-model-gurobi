#include <sys/stat.h>
#include <string>
#include "headers/Graph.h"
#include "headers/Model.h"
#include <chrono>

int main(int argc, const char *argv[]) {
  Graph *g = new Graph(argv[1]);
  Model *model = new Model(g);
  Model *model2 = new Model(g);
  model->createVariables();
  model->initModelExp();
  model->solveExp("3600");
  model->showSolution();

  cout << "------------------" << endl;
  model2->createVariables();
  model2->initModelCompact();
  model2->solveCompact("3600");
  model2->showSolution();
  return 0;
}
