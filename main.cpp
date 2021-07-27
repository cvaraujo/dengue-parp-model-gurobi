#include <sys/stat.h>
#include <string>
#include "headers/Graph.h"
#include "headers/Model.h"
#include <chrono>

int main(int argc, const char *argv[]) {
  Graph *g = new Graph(argv[1]);
  if(argv[3] == "e") {
  Model *model = new Model(g);
 
  model->createVariables();
  model->initModelExp();
  model->solveExp("3600");
  model->showSolution(argv[2]);
  } else {
  cout << "------------------" << endl;
  Model *model2 = new Model(g);
  model2->createVariables();
  model2->initModelCompact();
  model2->solveCompact("3600");
  model2->showSolution(argv[2]);
  }
  return 0;
}
