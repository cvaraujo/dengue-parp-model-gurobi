//
// Created by carlos on 06/07/21.
//

#include "../headers/Graph.h"

Graph::Graph(string instance)
{
  int block, cases, i, j, k;
  float length;
  string token, aux, x, y;
  ifstream file;
  vector<string> blcs;

  file.open(instance, fstream::in);

  file >> Graph::n >> Graph::m >> Graph::b;

  arcs = vector<vector<Arc *>>(n + 1, vector<Arc *>());
  nodesPerBlock = vector<vector<int>>(b, vector<int>());
  arcsPerBlock = vector<vector<Arc *>>(b, vector<Arc *>());

  for (i = 0; i < n; i++)
  {
    file >> token >> j >> x >> y >> aux;
    boost::split(blcs, aux, boost::is_any_of(","));
    vector<int> blocks_node = vector<int>();

    for (auto s : blcs)
    {
      if (s == "-1")
        break;
      k = stoi(s);

      blocks_node.push_back(k);
      nodesPerBlock[k].push_back(j);
    }
    nodes.push_back(make_pair(j, blocks_node));
  }

  nodes.push_back(make_pair(n, vector<int>()));

  vector<double> bigm_time;

  for (k = 0; k < m; k++)
  {
    file >> token >> i >> j >> length >> block >> cases;
    file.ignore(numeric_limits<streamsize>::max(), '\n');

    Arc *arc = new Arc(i, j, length, block, cases);

    bigm_time.push_back(timeArc(length, 350));

    arcs[i].push_back(arc);
    if (block != -1)
      arcsPerBlock[block].push_back(arc);
  }

  for (i = 0; i < n; i++)
  {
    arcs[n].push_back(new Arc(n, i, 0, -1, 0));
    arcs[i].push_back(new Arc(i, n, 0, -1, 0));
  }

  m_time = 0;
  for (i = 0; i < b; i++)
    m_time += timeBlock(250, i);

  sort(bigm_time.begin(), bigm_time.end(), greater<double>());
  for (i = 0; i < n - 1; i++)
    m_time += bigm_time[i];

  this->fillMissingArcs();

  cout << "M Time " << m_time << endl;

  cout << "Load graph successfully" << endl;
}

double Graph::getMtime()
{
  return m_time;
}

double Graph::timeArc(float distance, float speed)
{
  return distance > 0 ? distance / speed : 0;
}

double Graph::timeBlock(float speed, int block)
{
  float time = 0;
  for (auto *arc : arcsPerBlock[block])
  {
    time += timeArc(arc->getLength(), speed);
  }
  return time;
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
  for (int i = 0; i < n; i++)
    for (auto *arc : arcs[i])
      cout << "[" << i << ", " << arc->getD() << "] - " << arc->getLength() << ", " << arc->getBlock() << ", " << arc->getCases() << endl;
}

void Graph::fillMissingArcs()
{
  vector<double> distance;
  vector<int> pred;
  graph_t g = graph_t(n);

  int i;
  for (i = 0; i < n; i++)
    for (auto *arc : arcs[i])
      if (arc->getD() < n)
        add_edge(i, arc->getD(), arc->getLength(), g);

  for (i = 0; i < n; i++)
  {
    for (auto *arc : arcs[i])
    {
      if (arc->getD() >= n || this->exist_arc(arc->getD(), i))
        continue;
      distance = vector<double>(n);
      pred = vector<int>(n);
      dijkstra_shortest_paths(g, arc->getD(), predecessor_map(make_iterator_property_map(pred.begin(), get(vertex_index, g))).distance_map(make_iterator_property_map(distance.begin(), get(vertex_index, g))));
      // cout << "Distance: " << arc->getD() << ", " << i << ": " << distance[i] << " - REV: " << arc->getLength() << endl;
      if (distance[i] < 1.79769e+308)
        this->arcs[arc->getD()].push_back(new Arc(arc->getD(), i, distance[i], -1, 0));
      else
        this->arcs[arc->getD()].push_back(new Arc(arc->getD(), i, this->getMtime(), -1, 0));
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

int Graph::getRoot() const
{
  return n + 1;
}
