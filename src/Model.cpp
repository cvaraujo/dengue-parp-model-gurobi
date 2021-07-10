//
// Created by carlos on 06/07/21.
//

#include "../headers/Model.h"

class cyclecallback: public GRBCallback {

public:
  double lastiter, lastnode;
  int numvars;
  vector<vector<GRBVar>> x, y;

  Graph *graph;

  cyclecallback(Graph *xgraph, int xnumvars, vector<vector<GRBVar>> xx, vector<vector<GRBVar>> yy){
    lastiter = lastnode = 0;
    numvars = xnumvars;
    x = xx;
    y = yy;
    graph = xgraph;
  }

protected:
  void callback() {
    try {
      if(where == GRB_CB_MIPSOL) {
        int i, n = graph->getN();
        vector<vector<int>> g = vector<vector<int>>(n+1, vector<int>());

        for(i = 0; i <= n; i++)
          for (auto *arc : graph->arcs[i])
            if(getSolution(x[i][arc->getD()]) > 0.5) {
              cout << i << " -> " << arc->getD() << endl;
              g[i].push_back(arc->getD());
            }

        int s;
        vector<bool> visited(n+1, false);
        vector<vector<int>> conn = vector<vector<int>>(n+1, vector<int>());

        for(i = n; i >= 0; i--) {
          vector<int> stack;
          stack.push_back(i);

          while(!stack.empty()) {
            s = stack.back();
            stack.pop_back();

            if(!visited[s]) {
              conn[i].push_back(s);
              visited[s] = true;
            }

            for(auto k : g[s])
              if(!visited[k]) stack.push_back(k);
          }
        }

        for(i = 0; i < n; i++) {
          if(conn[i].size() > 1) {
            for(int u : conn[i]) {
              for(int b : graph->nodes[u].second) {
                if (getSolution(y[u][b]) > 0.5) {
                  GRBLinExpr in = 0, out = 0;
                  //from root
                  for(int v : conn[n]) {
                    for(auto *arc : graph->arcs[v])
                      if(arc->getD() == u) in += x[v][u];
                  }
                  // to root
                  for(int v : graph->nodesPerBlock[b]){
                    for(auto *arc : graph->arcs[v])
                      for(int k : conn[n])
                        if(arc->getD() == k) out += x[v][k];
                  }
                  addLazy(out >= y[u][b]);
                  addLazy(in >= y[u][b]);
                  break;
                }
              }
            }
          }
        }
        getchar();
      }
  } catch(GRBException e) {
    cout << "Error number: " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch (...) {
    cout << "Error during callback" << endl;
  }
  }
};

Model::Model(Graph *graph) {
  if (graph != nullptr) {
    this->graph = graph;
  } else exit(EXIT_FAILURE);
}

void Model::createVariables() {
  int o, d, k, n = graph->getN(), m = graph->getM(), b = graph->getB();
  try {
    env.set("LogFile", "MS_mip.log");
    env.start();

    x = vector<vector<GRBVar>>(n+2, vector<GRBVar>(n+2));
    y = vector<vector<GRBVar>>(n, vector<GRBVar>(b));
    t = vector<GRBVar>(n+2);

    // variables x
    char name[30];
    for (o = 0; o <= n; o++) {
      for (auto *arc : graph->arcs[o]) {
        d = arc->getD();
        sprintf(name, "x_%d_%d", o, d);
        x[o][d] = model.addVar(0.0, 1.0, 0, GRB_BINARY, name);
      }
    }

    // variables y
    for (int i = 0; i < n; i++) {
      o = graph->nodes[i].first;
      for (auto bl : graph->nodes[i].second) {
        sprintf(name, "y_%d_%d", o, bl);
        y[o][bl] = model.addVar(0.0, 1.0, 0, GRB_BINARY, name);
      }
    }

    // variables t
    for (int i = 0; i <= n+1; i++) {
      sprintf(name, "t_%d", i);
      t[i] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, name);
    }

    model.update();
    cout << "Create variables" << endl;
  } catch (GRBException &ex) {
    cout << ex.getMessage() << endl;
    cout << ex.getErrorCode() << endl;
    exit(EXIT_FAILURE);
  }
}

void Model::initModelExp() {
  cout << "Begin the model creation" << endl;
  objectiveFunction();
  artificialNodes(), flowConservation();
  maxAttending(), attendingPath(), timeConstraint(0.1);
  inseticideConstraint(20);
  cout << "All done!" << endl;
}

void Model::initModelCompact() {
  cout << "Begin the model creation" << endl;
  objectiveFunction();
  artificialNodes(), flowConservation();
  maxAttending(), attendingPath(), compactTimeConstraint(0.5);
  inseticideConstraint(20);
  cout << "All done!" << endl;
}

void Model::objectiveFunction() {
  GRBLinExpr objective;
  int i, j, n = graph->getN();

  for(int i = 0; i < n; i++) {
    j = graph->nodes[i].first;
    for(auto b : graph->nodes[i].second) {
      objective += y[j][b];
    }
  }

  model.setObjective(objective, GRB_MAXIMIZE);
  model.update();
  cout << "Objective Function" << endl;
}

void Model::artificialNodes() {
  int n = graph->getN();
  GRBLinExpr sink, target;

  for(int i = 0; i < n; i++) {
    sink += x[n][i];
    target += x[i][n+1];
  }
  model.addConstr(sink == 1, "sink_constraint");
  model.addConstr(target == 1, "target_constraint");
  cout << "Artificial nodes" << endl;
}

void Model::flowConservation() {
  int i, j, n = graph->getN();

  for (i = 0; i < n; i++) {
    GRBLinExpr flow_out, flow_in;
    for (auto *arc : graph->arcs[i]) {
      if(arc->getD() == n) continue;
      j = arc->getD();
      flow_out += x[i][j];
    }
    for (j = 0; j <= n; j++) {
      for (auto *arc : graph->arcs[j]) {
        if (arc->getD() == i) flow_in += x[j][i];
      }
    }
    model.addConstr(flow_in - flow_out == 0, "flow_conservation_" + to_string(i));
  }
  cout << "Flow conservation" << endl;
}

void Model::maxAttending() {
  int bl, b = graph->getB();

  for (bl = 0; bl < b; bl++) {
    GRBLinExpr maxServ;
    for(auto i : graph->nodesPerBlock[bl]) {
      maxServ += y[i][bl];
    }
    model.addConstr(maxServ <= 1, "max_service_block_" + to_string(bl));
  }
  cout << "Attending blocks" << endl;
}

void Model::attendingPath() {
  int j, bl, n = graph->getN(), b = graph->getB();

  for (bl = 0; bl < b; bl++) {
    for(auto i : graph->nodesPerBlock[bl]) {
      GRBLinExpr served;
      for(j = 0; j < n; j++)
        for(auto *arc : graph->arcs[j])
          if (arc->getD() == i) served += x[j][i];
      model.addConstr(served >= y[i][bl], "att_path_" + to_string(bl) + "_" + to_string(i));
    }
  }
  cout << "Make path to attending block" << endl;
}

void Model::timeConstraint(float maxTime) {
  int i, j, n = graph->getN();
  GRBLinExpr arcTravel, blockTravel;

  for(i = 0; i < n; i++) {
    for (auto *arc : graph->arcs[i]) {
      j = arc->getD();
      arcTravel += x[i][j] * timeArc(arc->getLength(), 30000);
    }
    for (auto b : graph->nodes[i].second) {
      blockTravel += y[i][b] * timeBlock(15000, b);
    }
  }
  model.addConstr(arcTravel + blockTravel <= maxTime, "max_time");
  cout << "Time constraint"<< endl;
}

void Model::compactTimeConstraint(float maxTime) {
  int i, j, n = graph->getN();
  model.addConstr(t[n] == 0);
  for(i = 0; i <= n; i++) {
    if (i == n) continue;
    for(auto *arc : graph->arcs[i]) {
      j = arc->getD();
      cout << j << endl;
      GRBLinExpr expr = 0;
      expr += t[i] + (timeArc(arc->getLength(), 1500) * x[i][j] - 100 * (1 - x[i][j]));
      if(arc->getBlock() != -1)
        expr += timeBlock(1500, arc->getBlock()) * y[i][arc->getBlock()];
      model.addConstr(t[j] >= expr);
    }
  }
  cout << "passei"<< endl;
  model.addConstr(t[n+1] <= maxTime);
  cout << "Time constraint done!" << endl;
}

void Model::inseticideConstraint(float maxInseticide) {
  int i, j, n = graph->getN();
  GRBLinExpr insConsumed;

  for(i = 0; i < n; i++) {
    for(auto b : graph->nodes[i].second) {
      insConsumed += y[i][b] * inseticideBlock(0.01, b);
    }
  }
  model.addConstr(insConsumed <= maxInseticide, "max_inseticide");
  cout << "Inseticide constraint" << endl;
}

float Model::timeArc(float distance, float speed) {
  return distance > 0 ? distance/speed : 0;
}

float Model::timeBlock(float speed, int block) {
  float time = 0;
  for (auto *arc : graph->arcsPerBlock[block]) {
    time += timeArc(arc->getLength(), speed);
  }
  return time;
}

float Model::inseticideBlock(float perMeter, int block) {
  float consumed = 0;
  for (auto *arc : graph->arcsPerBlock[block]) {
   consumed += arc->getLength() * perMeter;
  }
  return consumed;
}

void Model::solve(string timeLimit) {
  try {
    model.set("TimeLimit", timeLimit);
    //model.set(GRB_IntParam_LazyConstraints, 1);
    // model.set(GRB_DoubleParam_Heuristics, 0.0);
    //cyclecallback cb = cyclecallback(graph, graph->getN(), x, y);
    model.update();
    //model.setCallback(&cb);
    // model.computeIIS();
    model.set("OutputFlag", "0");
    model.write("model.lp");
    model.optimize();
  } catch (GRBException &ex) {
    cout << ex.getMessage() << endl;
  }
}

void Model::showSolution(){
  try {
    int n = graph->getN(), b = graph->getB();
    cout << "OF: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

    cout << "X" << endl;
    for(int i = 0; i <= n; i++){
      for(auto *arc : graph->arcs[i]) {
        if(x[i][arc->getD()].get(GRB_DoubleAttr_X) > 0) cout << i << " -- " << arc->getD() << endl;
      }
    }
    cout << "---------------------------------------------" << endl;
    cout << "Y" << endl;
    for (int i = 0; i < n; i++) {
      int o = graph->nodes[i].first;
      for (auto bl : graph->nodes[i].second) {
        if(y[i][bl].get(GRB_DoubleAttr_X) > 0) cout << i << " -> " << bl << endl;
      }
    }
  } catch (GRBException &ex) {
    cout << ex.getMessage() << endl;
  }
}
