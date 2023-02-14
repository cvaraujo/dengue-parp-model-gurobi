//
// Created by carlos on 06/07/21.
//

#include "../headers/Model.h"
#include "../headers/WarmStart.h"

using namespace lemon;

class cyclecallback : public GRBCallback
{

public:
  double lastiter, lastnode;
  int numvars, cuts = 0, num_frac_cuts = 0, num_lazy_cuts = 0;
  bool frac_cut = false;
  vector<vector<GRBVar>> x, y;
  typedef ListDigraph G;
  typedef G::Arc Arc;
  typedef G::ArcIt ArcIt;
  typedef G::Node Node;
  typedef G::ArcMap<double> LengthMap;
  typedef G::NodeMap<bool> BoolNodeMap;
  Graph *graph;

  cyclecallback(Graph *xgraph, int xnumvars, vector<vector<GRBVar>> xx, vector<vector<GRBVar>> yy, bool frac_cut)
  {
    lastiter = lastnode = 0;
    numvars = xnumvars;
    x = xx;
    y = yy;
    graph = xgraph;
    this->frac_cut = frac_cut;
  }

protected:
  void callback()
  {
    try {
      if (where == GRB_CB_MIPSOL) {
        int i, j, s, n = graph->getN();
        vector<vector<int>> g = vector<vector<int>>(n + 2, vector<int>());
        vector<bool> used_node = vector<bool>(n + 1);

        for (i = 0; i <= n; i++) {
          for (auto *arc : graph->arcs[i]) {
            if (getSolution(x[i][arc->getD()]) >= 0.1) {
              g[i].push_back(arc->getD());
              used_node[i] = used_node[arc->getD()] = true;
            }
          }
        }

        vector<bool> visited(n + 1, false);
        vector<int> conn_comp = vector<int>(n + 1, -1);
        vector<vector<int>> conn = vector<vector<int>>(n + 1, vector<int>());

        for (i = n; i >= 0; i--) {
          if (!used_node[i]) continue;

          vector<int> stack;
          stack.push_back(i);

          while (!stack.empty()) {
            s = stack.back();
            stack.pop_back();

            if (!visited[s]) {
              conn[i].push_back(s);
              conn_comp[s] = i;
              visited[s] = true;
            }

            for (auto k : g[s])
              if (!visited[k]) stack.push_back(k);
          }
        }

        int num_comp = 0;
        for (i = n; i >= 0; i--) {
          if (conn[i].size() > 1) num_comp++;
          if (num_comp > 1) break;
        }

        if (num_comp == 1) return;

        for (i = 0; i <= n; i++) {
          if (conn[i].size() <= 1) continue;
          GRBLinExpr expr = 0;
          bool has_constr = false;

          for (auto v : conn[i]) {
            for (auto *arc : graph->arcs[v]) {
              j = arc->getD();
              if (conn_comp[j] != i) {
                expr += x[v][j];
                has_constr = true;
              }
            }
          }

          if (has_constr) {
            for (auto v : conn[i]) {
              for (auto b : graph->nodes[v].second) {
                if (getSolution(y[v][b]) >= 0.1) {
                  addLazy(expr >= y[v][b]);
                  num_lazy_cuts++;
                }
              }
            }
          }
        }
      }
      else if (where == GRB_CB_MIPNODE)
      {
        if(!frac_cut) return;

        int mipStatus = getIntInfo(GRB_CB_MIPNODE_STATUS);
        double nodecnt = getDoubleInfo(GRB_CB_MIP_NODCNT);

        if (mipStatus == GRB_OPTIMAL) {
          int i, j, u, v, n = graph->getN();

          // Basic structures to use Lemon
          G flow_graph;
          LengthMap capacity(flow_graph);
          vector<Node> setNodes = vector<Node>(n + 1);
          vector<bool> usedNode = vector<bool>(n, false);
          vector<Arc> setArcs;

          // Create the node set
          for (i = 0; i <= n; i++) setNodes[i] = flow_graph.addNode();

          // Create the edge set
          for (i = 0; i < n; i++) {
            for (auto *arc : graph->arcs[i]) {
              j = arc->getD();

              if (getNodeRel(x[i][j]) >= 0.1) {
                setArcs.push_back(flow_graph.addArc(setNodes[i], setNodes[j]));
                capacity[setArcs[setArcs.size() - 1]] = double(getNodeRel(x[i][j]));
                usedNode[i] = usedNode[j] = true;
              }
            }
          }

          // Init necessary structures
          double mincut_value;
          bool source_side, has_constr, need_cut;

          // cout << "Add Cut" << endl;
          for (i = 0; i < n; i++)
          {
            // If there is no arc using this node, ignore it
            if (!usedNode[i]) continue;

            // Lemon MaxFlow instance
            Preflow<G, LengthMap> preflow(flow_graph, capacity, setNodes[i], setNodes[n]);
            preflow.runMinCut();
            mincut_value = preflow.flowValue();

            // Create basic variables
            // BoolNodeMap bm(flow_graph);
            source_side = preflow.minCut(setNodes[i]);
            
            // Check the need for cutting
            need_cut = false;
            for (auto b : graph->nodes[i].second) {
              if (mincut_value < getNodeRel(y[i][b])) {
                need_cut = true;
                break;
              }
            }

            if (!need_cut) continue;

            GRBLinExpr expr = 0;
            has_constr = false;

            for (u = 0; u < n; u++)
            {
              if (preflow.minCut(setNodes[u]) != source_side)
                continue;

              for (auto *arc : graph->arcs[u]) {
                v = arc->getD();
                if (preflow.minCut(setNodes[v]) != source_side) {
                  expr += x[u][v];
                  has_constr = true;
                }
              }
            }

            if (has_constr) {
              for (u = 0; u < n; u++) {
                if (usedNode[u] && preflow.minCut(setNodes[u])) {
                  usedNode[u] = false;
                  for (auto b : graph->nodes[u].second) {
                    if (getNodeRel(y[u][b]) > mincut_value) {
                      addCut(expr >= y[u][b]);
                      num_frac_cuts++;
                    }
                  }
                }
              }
            } 
          }
        }
      }
    }
    catch (GRBException e)
    {
      cout << "Error number: " << e.getErrorCode() << endl;
      cout << e.getMessage() << endl;
    }
    catch (...)
    {
      cout << "Error during callback" << endl;
    }
  }
};

Model::Model(Graph *graph, float default_vel, float spraying_vel, float insecticide_ml_min)
{
  if (graph != nullptr) {
    this->graph = graph;
    this->default_vel = default_vel;
    this->spraying_vel = spraying_vel;
    this->insecticide_ml_min = insecticide_ml_min;
  }
  else exit(EXIT_FAILURE);
}

void Model::createVariables()
{
  int o, d, k, n = graph->getN(), m = graph->getM(), b = graph->getB();
  try
  {
    env.set("LogFile", "MS_mip.log");
    env.start();

    x = vector<vector<GRBVar>>(n + 1, vector<GRBVar>(n + 1));
    y = vector<vector<GRBVar>>(n, vector<GRBVar>(b));
    t = vector<vector<GRBVar>>(n + 1, vector<GRBVar>(n + 1));

    // X
    char name[40];
    for (o = 0; o <= n; o++) {
      for (auto *arc : graph->arcs[o]) {
        d = arc->getD();
        sprintf(name, "x_%d_%d", o, d);
        x[o][d] = model.addVar(0.0, 1.0, 0, GRB_BINARY, name);
      }
    }

    // Y
    for (int i = 0; i < n; i++) {
      o = graph->nodes[i].first;
      for (auto bl : graph->nodes[i].second) {
        sprintf(name, "y_%d_%d", o, bl);
        y[o][bl] = model.addVar(0.0, 1.0, 0, GRB_BINARY, name);
      }
    }

    // T
    for (o = 0; o <= n; o++) {
      for (auto *arc : graph->arcs[o]) {
        d = arc->getD();
        sprintf(name, "t_%d_%d", o, d);
        t[o][d] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, name);
      }
    }
    model.update();
    cout << "Create variables" << endl;
  }
  catch (GRBException &ex)
  {
    cout << ex.getMessage() << endl;
    cout << ex.getErrorCode() << endl;
    exit(EXIT_FAILURE);
  }
}

void Model::WarmStart(float maxTime, float maxInsecticide)
{
  int i, j;

  // Warm start
  vector<pair<int, int>> x, y;
  double of = WarmStart::compute_solution(graph, maxTime, maxInsecticide, x, y, default_vel, spraying_vel, insecticide_ml_min);

  cout << "[***] Heuristic value = " << of << endl;

  for (int i = 0; i <= graph->getN(); i++)
  {
    for (auto *arc : graph->arcs[i])
    {
      bool is_in = false;

      for (auto pair : x)
      {
        int k = pair.first, j = pair.second;
        if (k == i && j == arc->getD())
        {
          is_in = true;
          break;
        }
      }
      if (is_in)
        this->x[i][arc->getD()].set(GRB_DoubleAttr_Start, 1.0);
      else
        this->x[i][arc->getD()].set(GRB_DoubleAttr_Start, 0.0);
    }
  }

  model.update();

  for (auto pair : y)
  {
    i = pair.first, j = pair.second;
    this->y[i][j].set(GRB_DoubleAttr_Start, 1.0);
  }

  model.update();
}

void Model::initModelExp(float maxTime, float maxInsecticide, bool warm_start)
{
  cout << "[!!!] Creating the model!" << endl;
  objectiveFunction();
  artificialNodes(), flowConservation();
  maxAttending(), attendingPath(), timeConstraint(maxTime);
  inseticideConstraint(maxInsecticide);

  if (warm_start) {
    cout << "[!!!] Calling Warm-Start function!" << endl;
    this->WarmStart(maxTime, maxInsecticide);
  }

  cout << "[***] done!" << endl;
}

void Model::initModelCompact(float maxTime, float maxInsecticide)
{
  cout << "[!!!] Creating the model!" << endl;  objectiveFunction();
  artificialNodes(), flowConservation();
  maxAttending(), attendingPath(), compactTimeConstraint(maxTime);
  inseticideConstraint(maxInsecticide);
  cout << "[***] done!" << endl;
}

void Model::objectiveFunction()
{
  GRBLinExpr objective;
  int i, j, n = graph->getN();

  for (i = 0; i < n; i++) {
    j = graph->nodes[i].first;
    for (auto b : graph->nodes[i].second) {
      objective += (y[j][b] * profitBlock(b));
    }
  }

  model.setObjective(objective, GRB_MAXIMIZE);
  model.update();
  cout << "[***] Obj. Function: Maximize profit" << endl;
}

void Model::artificialNodes()
{
  int n = graph->getN();
  GRBLinExpr sink, target;

  for (int i = 0; i < n; i++) {
    sink += x[n][i];
    target += x[i][n];
  }

  model.addConstr(sink == 1, "sink_constraint");
  model.addConstr(target == 1, "target_constraint");
  cout << "[***] Contraint: dummy depot" << endl;
}

void Model::flowConservation()
{
  int i, j, n = graph->getN();

  for (i = 0; i < n; i++) {
    GRBLinExpr flow_out, flow_in;
    
    for (auto *arc : graph->arcs[i]) {
      if (arc->getD() >= n) continue;
      flow_out += x[i][arc->getD()];
    }

    for (j = 0; j < n; j++) {
      for (auto *arc : graph->arcs[j]) {
        if (arc->getD() == i)
          flow_in += x[j][i];
      }
    }

    flow_out += x[i][n];
    flow_in += x[n][i];
    model.addConstr(flow_in - flow_out == 0, "flow_conservation_" + to_string(i));
  }

  cout << "[***] Constraint: Flow conservation" << endl;
}

void Model::maxAttending()
{
  int bl, b = graph->getB();

  for (bl = 0; bl < b; bl++) {
    GRBLinExpr maxServ;
    for (auto i : graph->nodesPerBlock[bl]) {
      maxServ += y[i][bl];
    }

    model.addConstr(maxServ <= 1, "max_service_block_" + to_string(bl));
  }

  cout << "[***] Constraint: Serve each block at most once" << endl;
}

void Model::attendingPath()
{
  int j, bl, n = graph->getN(), b = graph->getB();

  for (bl = 0; bl < b; bl++) {
    for (auto i : graph->nodesPerBlock[bl]) {
      GRBLinExpr served;
      for (auto *arc : graph->arcs[i]) served += x[i][arc->getD()];
      
      model.addConstr(served >= y[i][bl], "att_path_" + to_string(i) + "_" + to_string(bl));
    }
  }

  cout << "[***] Constraint: Include node in path" << endl;
}

void Model::timeConstraint(float maxTime)
{
  int i, j, n = graph->getN();
  vector<float> block_time = vector<float>(graph->getB(), -1);
  GRBLinExpr arcTravel, blockTravel;

  for (i = 0; i < n; i++) {
    for (auto *arc : graph->arcs[i]) {
      j = arc->getD();
      arcTravel += x[i][j] * timeArc(arc->getLength(), this->default_vel);
    }

    for (auto b : graph->nodes[i].second) {
      if(block_time[b] == -1) block_time[b] = timeBlock(b, this->spraying_vel);
      blockTravel += y[i][b] * block_time[b];
    }
  }

  model.addConstr(arcTravel + blockTravel <= maxTime, "max_time");
  cout << "[***] Constraint: time limit" << endl;
}

void Model::compactTimeConstraint(float maxTime)
{
  int b, i, j, k, n = graph->getN();
  vector<float> block_time = vector<float>(graph->getB(), -1);   

  for (i = 0; i <= n; i++) {
    if(i < n) model.addConstr(t[n][i] == 0);

    for (auto *arc : graph->arcs[i]) {
      j = arc->getD();
      if (j >= n) continue;

      GRBLinExpr time_ij = 0;
      time_ij += t[i][j] + (timeArc(arc->getLength(), this->default_vel) * x[i][j]);
      
      for(auto b : graph->nodes[j].second) {
        if(block_time[b] == -1) block_time[b] = timeBlock(b, this->spraying_vel);
        time_ij += y[j][b] * block_time[b];
      }

      for (auto *arcl : graph->arcs[j]) {
        k = arcl->getD();
        model.addConstr(time_ij >= t[j][k] - ((2 - x[i][j] - x[j][k]) * maxTime), "t_geq_" + to_string(i) + "_" + to_string(j) + "_" + to_string(k));
        model.addConstr(time_ij <= t[j][k] + ((2 - x[i][j] - x[j][k]) * maxTime), "t_leq_" + to_string(i) + "_" + to_string(j) + "_" + to_string(k));
      }
    }
  }
  for (i = 0; i < n; i++) {
    model.addConstr(t[i][n] <= maxTime, "max_time");
  }
  cout << "[***] Constraint: Time limit" << endl;
}

void Model::inseticideConstraint(float maxInseticide)
{
  int i, j, n = graph->getN();
  GRBLinExpr insConsumed;

  for (i = 0; i < n; i++) {
    for (auto b : graph->nodes[i].second) {
      insConsumed += y[i][b] * inseticideBlock(b, this->insecticide_ml_min);
    }
  }

  model.addConstr(insConsumed <= maxInseticide, "max_inseticide");
  cout << "[***] Constraint: Insecticide limit" << endl;
}

float Model::timeArc(float distance, float speed)
{
  return distance > 0 ? distance / ((speed * 1000)/60) : 0;
}

float Model::timeBlock(int block, float speed)
{
  float time = 0;
  for (auto *arc : graph->arcsPerBlock[block]) {
    time += timeArc(arc->getLength(), speed);
  }
  return time;
}

float Model::inseticideBlock(int block, float perMeter)
{
  float consumed = 0;
  for (auto *arc : graph->arcsPerBlock[block]) {
    consumed += timeArc(arc->getLength(), this->spraying_vel) * perMeter;
  }
  return consumed;
}

float Model::profitBlock(int block)
{
  int i, n = graph->getN(), total = 0;
  for (auto *arc : graph->arcsPerBlock[block]) total += arc->getCases();
  return total;
}

void Model::solveCompact(string timeLimit)
{
  try {
    model.set("TimeLimit", timeLimit);

    model.update();
    model.set("OutputFlag", "0");
    model.write("model.lp");
    model.optimize();
  } catch (GRBException &ex) {
    cout << ex.getMessage() << endl;
  }
}

void Model::solveExp(string timeLimit, bool frac_cut)
{
  try {
    model.set("TimeLimit", timeLimit);

    model.set(GRB_DoubleParam_Heuristics, 1.0);
    model.set(GRB_IntParam_LazyConstraints, 1);
    cyclecallback cb = cyclecallback(graph, graph->getN(), x, y, frac_cut);
    model.setCallback(&cb);

    model.update();
    model.set("OutputFlag", "0");
    model.write("model.lp");
    model.optimize();

    // Save the number of cuts
    num_lazy_cuts = cb.num_lazy_cuts, num_frac_cuts = cb.num_frac_cuts;
  } catch (GRBException &ex) {
    cout << ex.getMessage() << endl;
  }
}

void Model::writeSolution(string result)
{
  try {
    ofstream output;
    output.open(result);

    int j, n = graph->getN(), b = graph->getB();

    output << "Nodes: " << n << endl;
    output << "Arcs: " << graph->getM() << endl;
    output << "Blocks: " << b << endl;
    output << "LB: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
    output << "UB: " << model.get(GRB_DoubleAttr_ObjBound) << endl;
    output << "Gurobi Nodes: " << model.get(GRB_DoubleAttr_NodeCount) << endl;
    output << "LAZY_CUTS: " << this->num_lazy_cuts << endl;
    output << "FRAC_CUTS: " << this->num_frac_cuts << endl;
    output << "Runtime: " << model.get(GRB_DoubleAttr_Runtime) << endl;

    float timeUsed = 0, insecUsed = 0;
    for (int i = 0; i < n; i++) {
      for (auto *arc : graph->arcs[i]) {
        j = arc->getD();
        if (x[i][j].get(GRB_DoubleAttr_X) > 0.5) {
          timeUsed += timeArc(arc->getLength(), this->default_vel);
          output << "X: " << i << " " << j << endl;
        }
      }

      for (auto b : graph->nodes[i].second) {
        if (y[i][b].get(GRB_DoubleAttr_X) > 0.5) {
          timeUsed += timeBlock(b, this->spraying_vel);
          insecUsed += inseticideBlock(b, this->insecticide_ml_min);
          output << "Y: " << i << " " << b << endl;
        }
      }
    }

    output << "Route Time: " << timeUsed << endl;
    output << "Insecticide Used: " << insecUsed << endl;

    cout << "OF: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
  } catch (GRBException &ex) {
    ofstream output;
    output.open(result);

    int n = graph->getN(), b = graph->getB(), j;
    output << "Nodes: " << n << endl;
    output << "Arcs: " << graph->getM() << endl;
    output << "Blocks: " << b << endl;
    output << "LB: 0" << endl;
    output << "UB: " << model.get(GRB_DoubleAttr_ObjBound) << endl;
    output << "N. Nodes: " << model.get(GRB_DoubleAttr_NodeCount) << endl;
    output << "Runtime: " << model.get(GRB_DoubleAttr_Runtime) << endl;
    cout << "OF: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
    cout << ex.getMessage() << endl;
  }
}

bool Model::check_solution(float max_time, float max_insecticide)
{
  int n = graph->getN();

  // Check connectivity
  vector<vector<bool>> used_arc = vector<vector<bool>>(n + 1, vector<bool>(n + 1, false));

  int start_node = n, i, j, s, target;
  bool find_next = true;
  float time = 0, insecticide = 0;

  vector<bool> visited(n + 1, false);
  vector<int> conn_comp = vector<int>(n + 1, -1);
  vector<vector<int>> conn = vector<vector<int>>(n + 1, vector<int>());

  // DFS
  deque<int> stack;
  stack.push_back(n);

  while (!stack.empty()) {
    s = stack.front();
    stack.pop_front();

    for (auto *arc : graph->arcs[s])
    {
      j = arc->getD();
      if (x[s][j].get(GRB_DoubleAttr_X) > 0.5)
      {
        used_arc[s][j] = true;

        if (!visited[j])
        {
          stack.push_back(j);
          visited[j] = true;
        }
      }
    }
  }

  // Check visiting
  for (i = 0; i <= n; i++) {
    for (auto b : graph->nodes[i].second) {
      if (y[i][b].get(GRB_DoubleAttr_X) > 0.5) {
        time += timeBlock(b, this->spraying_vel);
        insecticide += inseticideBlock(b, this->insecticide_ml_min);
        
        if (!visited[i]) {
          cout << "[!!!] Not visited node!" << endl;
          return false;
        }
      }
    }

    for (auto *arc : graph->arcs[i]) {
      if (x[i][arc->getD()].get(GRB_DoubleAttr_X) > 0.5) {
        time += timeArc(arc->getLength(), this->default_vel);
        if(!used_arc[i][arc->getD()]) {
          cout << "[!!!] Not used arc!" << endl;
          return false;
        }
      }
    }
  }

  if(time > max_time || insecticide > max_insecticide) {
    cout << "T: " << time << " <= " << max_time << ", I: " << insecticide << " <= " << max_insecticide << endl;
    cout << "[!!!] Resource limitation error!" << endl;
    return false;
  }

  cout << "[***] Instance ok!!!" << endl;
  return true;
}
