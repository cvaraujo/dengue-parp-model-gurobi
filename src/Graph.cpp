//
// Created by carlos on 06/07/21.
//

#include "../headers/Graph.h"
#include <map>

Graph::Graph(string instance, int graph_adapt, int instance_type)
{
  if (instance_type == 1)
    loadSBRPInstance(instance, graph_adapt);
  else
    load_instance(instance, graph_adapt);
}

Graph::Graph(string instance, string scenarios, int graph_adapt)
{
  load_stochastic_instance(instance, graph_adapt);
  load_scenarios_instance(scenarios);
}

void Graph::loadSBRPInstance(string instance, int graph_adapt)
{
  int block, cases, i, j, k;
  float length;
  string token, src, tgt, aux, x, y;
  ifstream file;
  vector<string> vertices, blcs;

  map<string, int> node_map;
  vector<pair<string, string>> arcsMap = vector<pair<string, string>>();
  vector<float> lenghts = vector<float>();

  file.open(instance, fstream::in);

  while (!file.eof())
  {
    file >> token;
    if (token == "VERTICES")
    {
      file >> token >> Graph::n;
      arcs = vector<vector<Arc *>>(n + 1, vector<Arc *>());
    }
    else if (token == "ARCOS_NOREQ")
    {
      file >> token >> Graph::m;
    }
    else if (token == "BLOQUES")
    {
      file >> token >> Graph::b;
      nodesPerBlock = vector<set<int>>(b, set<int>());
      arcsPerBlock = vector<vector<Arc *>>(b, vector<Arc *>());
    }
    else if (token == "LISTA_ARISTAS_REQ")
    {
      file >> token;
      for (i = 0; i < m; i++)
      {
        file >> src >> tgt >> token >> length;
        src.erase(src.begin());
        src.pop_back();
        tgt.pop_back();
        arcsMap.push_back(make_pair(src, tgt));
        lenghts.push_back(length);
      }
    }
    else if (token == "LISTA_VERTICES")
    {
      file >> token;
      for (i = 0; i < n; i++)
      {
        file >> src >> token >> token;
        node_map[src] = i;
        nodes.push_back(make_pair(i, set<int>()));
      }
    }
    else if (token == "LISTA_BLOQUE")
    {
      // Create arcs
      i = 0;
      for (auto p : arcsMap)
      {
        Arc *arc = new Arc(node_map[p.first], node_map[p.second], lenghts[i++], -1, 0);
        arcs[node_map[p.first]].push_back(arc);
      }

      // Ignoring the end of the line
      file.ignore(numeric_limits<streamsize>::max(), '\n');

      // Check the nodes of each block
      for (i = 0; i < b; i++)
      {

        getline(file, token);
        boost::split(vertices, token, boost::is_any_of(","));
        cases = stoi(vertices.back());
        vertices.pop_back();

        for (auto s : vertices)
        {
          trim(s);
          j = node_map[s];
          nodesPerBlock[i].insert(j);
          nodes[j].second.insert(i);
        }

        for (k = 1; k < vertices.size(); k++)
        {
          trim(vertices[k - 1]);
          trim(vertices[k]);

          int source = node_map[vertices[k - 1]];
          int target = node_map[vertices[k]];

          for (auto *arc : arcs[source])
          {
            if (arc->getD() == target)
            {
              arc->setBlock(i);
              arcsPerBlock[i].push_back(arc);
              break;
            }
          }

          if (k == vertices.size() - 1)
          {
            int first_node = node_map[vertices[0]];
            for (auto *arc : arcs[target])
            {
              if (arc->getD() == first_node)
              {
                arc->setCases(cases);
                arc->setBlock(i);
                arcsPerBlock[i].push_back(arc);
                break;
              }
            }
          }
        }
      }
    }
  }

  if (graph_adapt == 1 || graph_adapt == 3)
  {
    this->fillMissingArcs();
  }
  else if (graph_adapt == 2)
  {
    this->fillCompleteDigraph();
  }

  nodes.push_back(make_pair(n, set<int>()));

  for (i = 0; i < n; i++)
  {
    arcs[n].push_back(new Arc(n, i, 0, -1, 0));
    arcs[i].push_back(new Arc(i, n, 0, -1, 0));
  }

  cout << "Load graph successfully" << endl;
}

// Graph adapt
// 0 - No adapt
// 1 - Create reverse missing arcs
// 2 - Complete Digraph
// 3 - MTZ adaptation
void Graph::load_instance(string instance, int graph_adapt)
{
  int block, cases, i, j, k;
  float length;
  string token, aux, x, y;
  ifstream file;
  vector<string> blcs;

  file.open(instance, fstream::in);

  file >> Graph::n >> Graph::m >> Graph::b;

  if (graph_adapt == 3)
    arcs = vector<vector<Arc *>>(n + m + 1, vector<Arc *>());
  else
    arcs = vector<vector<Arc *>>(n + 1, vector<Arc *>());

  nodesPerBlock = vector<set<int>>(b, set<int>());
  arcsPerBlock = vector<vector<Arc *>>(b, vector<Arc *>());
  cases_per_block = vector<int>(b, 0);
  set<int> blocks_node;

  for (i = 0; i < n; i++)
  {
    file >> token >> j >> x >> y >> aux;
    boost::split(blcs, aux, boost::is_any_of(","));
    blocks_node = set<int>();

    for (auto s : blcs)
    {
      if (s == "-1")
        break;
      k = stoi(s);

      blocks_node.insert(k);
      nodesPerBlock[k].insert(j);
    }
    nodes.push_back(make_pair(j, blocks_node));
  }

  vector<double> bigm_time;
  double new_length;
  int new_num_cases;
  for (k = 0; k < m; k++)
  {
    file >> token >> i >> j >> length >> block;
    file.ignore(numeric_limits<streamsize>::max(), '\n');

    if (graph_adapt == 3)
    {
      blocks_node = set<int>();

      if (block != -1)
      {
        blocks_node.insert(block);
        nodesPerBlock[block].insert(n);
      }

      nodes.push_back(make_pair(n, blocks_node));

      new_length = length / 2;
      new_num_cases = int(cases / 2);

      Arc *arc = new Arc(i, n, new_length, block, 0);
      Arc *sec_arc = new Arc(n, j, new_length, block, 0);

      arcs[i].push_back(arc);
      arcs[n].push_back(sec_arc);

      if (block != -1)
      {
        arcsPerBlock[block].push_back(arc);
        arcsPerBlock[block].push_back(sec_arc);
      }
      n++;
    }
    else
    {
      Arc *arc = new Arc(i, j, length, block, cases);

      arcs[i].push_back(arc);
      if (block != -1)
        arcsPerBlock[block].push_back(arc);
    }
  }

  while (!file.eof())
  {
    file >> token;

    if (token == "B")
    {
      file >> i >> j;
      cases_per_block[i] = j;
    }
  }

  if (graph_adapt == 1 || graph_adapt == 3)
  {
    this->fillMissingArcs();
  }
  else if (graph_adapt == 2)
  {
    this->fillCompleteDigraph();
  }

  nodes.push_back(make_pair(n, set<int>()));

  for (i = 0; i < n; i++)
  {
    arcs[n].push_back(new Arc(n, i, 0, -1, 0));
    arcs[i].push_back(new Arc(i, n, 0, -1, 0));
  }

  cout << "Load graph successfully" << endl;
}

void Graph::load_stochastic_instance(string instance, int graph_adapt)
{
  int block, cases, i, j, k;
  float length;
  string token, aux, x, y;
  ifstream file;
  vector<string> blcs;

  file.open(instance, fstream::in);

  file >> Graph::n >> Graph::m >> Graph::b;

  if (graph_adapt == 3)
  {
    arcs = vector<vector<Arc *>>(n + m + 1, vector<Arc *>());
  }
  else
  {
    arcs = vector<vector<Arc *>>(n + 1, vector<Arc *>());
  }

  nodesPerBlock = vector<set<int>>(b, set<int>());
  arcsPerBlock = vector<vector<Arc *>>(b, vector<Arc *>());
  cases_per_block = vector<int>(b);
  set<int> blocks_node;

  for (i = 0; i < n; i++)
  {
    file >> token >> j >> x >> y >> aux;
    boost::split(blcs, aux, boost::is_any_of(","));
    blocks_node = set<int>();

    for (auto s : blcs)
    {
      if (s == "-1")
        break;
      k = stoi(s);

      blocks_node.insert(k);
      nodesPerBlock[k].insert(j);
    }
    nodes.push_back(make_pair(j, blocks_node));
  }

  vector<double> bigm_time;
  double new_length;
  int new_num_cases;
  for (k = 0; k < m; k++)
  {
    file >> token >> i >> j >> length >> block;
    file.ignore(numeric_limits<streamsize>::max(), '\n');

    if (graph_adapt == 3)
    {
      blocks_node = set<int>();

      if (block != -1)
      {
        blocks_node.insert(block);
        nodesPerBlock[block].insert(n);
      }

      nodes.push_back(make_pair(n, blocks_node));

      new_length = length / 2;
      new_num_cases = int(cases / 2);

      Arc *arc = new Arc(i, n, new_length, block, 0);
      Arc *sec_arc = new Arc(n, j, new_length, block, 0);

      arcs[i].push_back(arc);
      arcs[n].push_back(sec_arc);

      if (block != -1)
      {
        arcsPerBlock[block].push_back(arc);
        arcsPerBlock[block].push_back(sec_arc);
      }
      n++;
    }
    else
    {
      Arc *arc = new Arc(i, j, length, block, cases);

      arcs[i].push_back(arc);
      if (block != -1)
        arcsPerBlock[block].push_back(arc);
    }
  }

  for (k = 0; k < b; k++)
  {
    file >> token >> i >> j;
    file.ignore(numeric_limits<streamsize>::max(), '\n');

    cases_per_block[i] = j;
  }

  if (graph_adapt == 1 || graph_adapt == 3)
  {
    this->fillMissingArcs();
  }
  else if (graph_adapt == 2)
  {
    this->fillCompleteDigraph();
  }

  nodes.push_back(make_pair(n, set<int>()));

  for (i = 0; i < n; i++)
  {
    arcs[n].push_back(new Arc(n, i, 0, -1, 0));
    arcs[i].push_back(new Arc(i, n, 0, -1, 0));
  }

  cout << "Load graph successfully" << endl;
}

void Graph::load_scenarios_instance(string instance)
{
  string token;
  ifstream file;
  int i, block, cases;
  float prob;

  file.open(instance, fstream::in);
  file >> Graph::s;
  Graph::scenarios = vector<Scenario>(Graph::s);

  while (!file.eof())
  {
    file >> token;

    if (token == "P")
    {
      file >> i >> prob;
      vector<int> cases_per_block = vector<int>(Graph::b, 0);
      Scenario scn(prob, cases_per_block);
      Graph::scenarios[i] = scn;
    }
    else if (token == "B")
    {
      file >> i >> block >> cases;
      Graph::scenarios[i].SetCases(block, cases);
    }
  }
}

double Graph::getMtime()
{
  return m_time;
}

float Graph::timeArc(float distance, float speed)
{
  return distance > 0 ? distance / ((speed * 1000) / 60) : 0;
}

float Graph::timeBlock(int block, float speed)
{
  float time = 0;
  for (auto *arc : arcsPerBlock[block])
  {
    time += timeArc(arc->getLength(), speed);
  }
  return time;
}

float Graph::inseticideBlock(int block, float perMeter)
{
  float consumed = 0;
  for (auto *arc : arcsPerBlock[block])
  {
    consumed += timeArc(arc->getLength(), this->spraying_vel) * perMeter;
  }
  return consumed;
}

bool Graph::exist_arc(int i, int j)
{
  if (i > n)
    return false;
  for (auto *arc : arcs[i])
  {
    if (arc->getD() == j)
    {
      return true;
    }
  }
  return false;
}

void Graph::showGraph()
{
  for (int i = 0; i <= n; i++)
    for (auto *arc : arcs[i])
      cout << "[" << i << ", " << arc->getD() << "] - " << arc->getLength() << ", " << arc->getBlock() << ", " << arc->getCases() << endl;
}

void Graph::showScenarios()
{

  for (auto scenario : this->scenarios)
  {
    cout << "Probability: " << scenario.probability << endl;
    for (int bl = 0; bl < this->b; bl++)
    {
      cout << "B" << bl << " -> " << scenario.cases_per_block[bl] << endl;
    }
  }
}

void Graph::fillMissingArcs()
{
  vector<double> distance;
  vector<int> pred;
  graph_t g = graph_t(n);

  int i;
  for (i = 0; i < n; i++)
    for (auto *arc : arcs[i])
      add_edge(i, arc->getD(), arc->getLength(), g);

  for (i = 0; i < n; i++)
  {
    for (auto *arc : arcs[i])
    {
      if (this->exist_arc(arc->getD(), i))
        continue;

      distance = vector<double>(n);
      pred = vector<int>(n);

      dijkstra_shortest_paths(g, arc->getD(), predecessor_map(make_iterator_property_map(pred.begin(), get(vertex_index, g))).distance_map(make_iterator_property_map(distance.begin(), get(vertex_index, g))));
      if (distance[i] < 1.79769e+308)
        this->arcs[arc->getD()].push_back(new Arc(arc->getD(), i, distance[i], -1, 0));
    }
  }
}

void Graph::init_boost_graph()
{
  boost_graph = graph_t(n);

  for (int i = 0; i < n; i++)
  {
    for (auto *arc : arcs[i])
    {
      if (arc->getD() >= n)
        continue;
      add_edge(i, arc->getD(), arc->getLength(), boost_graph);
    }
  }
}

double Graph::shortest_path(int i, int j, vector<pair<int, int>> &arcs)
{
  vector<double> distance = vector<double>(n);
  vector<int> pred = vector<int>(n);
  arcs = vector<pair<int, int>>();

  dijkstra_shortest_paths(boost_graph, i,
                          predecessor_map(make_iterator_property_map(pred.begin(),
                                                                     get(vertex_index, boost_graph)))
                              .distance_map(make_iterator_property_map(distance.begin(),
                                                                       get(vertex_index, boost_graph))));
  if (distance[j] < 1.79769e+308)
  {
    int k = j;
    while (k != i)
    {
      arcs.push_back(make_pair(pred[k], k));
      k = pred[k];
    }
    reverse(arcs.begin(), arcs.end());
    return distance[j];
  }
  return -1;
}

void Graph::fillCompleteDigraph()
{
  vector<double> distance;
  vector<int> pred;
  graph_t g = graph_t(n);

  int i, j;
  for (i = 0; i < n; i++)
    for (auto *arc : arcs[i])
      add_edge(i, arc->getD(), arc->getLength(), g);

  for (i = 0; i < n; i++)
  {
    distance = vector<double>(n);
    pred = vector<int>(n);
    dijkstra_shortest_paths(g, i, predecessor_map(make_iterator_property_map(pred.begin(), get(vertex_index, g))).distance_map(make_iterator_property_map(distance.begin(), get(vertex_index, g))));

    for (j = 0; j < n; j++)
    {
      if (!this->exist_arc(i, j))
      {
        if (distance[j] < 1.79769e+308)
          this->arcs[i].push_back(new Arc(i, j, distance[j], -1, 0));
      }
    }
  }
}

int Graph::getN() const
{
  return n;
}

int Graph::getM() const
{
  return m;
}

int Graph::getB() const
{
  return b;
}

int Graph::getS() const
{
  return s;
}

int Graph::getRoot() const
{
  return n + 1;
}
