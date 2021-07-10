#include <sys/stat.h>
#include <string>
#include "headers/Graph.h"
#include "headers/Model.h"
#include <chrono>

int main(int argc, const char *argv[]) {
  Graph *g = new Graph("tab.txt");
  g->showGraph();
  Model * model = new Model(g);
  model->createVariables();
  model->initModelCompact();
  model->solve("3600");
  model->showSolution();
  return 0;
}
