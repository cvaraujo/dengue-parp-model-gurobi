//
// Created by carlos on 06/07/21.
//

#include "../headers/Graph.h"

Graph::Graph(string instance) {
  int block, cases, i, j, k;
  float length;
  string token, aux;
  ifstream file;
  vector<string> blcs;

  file.open(instance, fstream::in);

  file >> Graph::n >> Graph::m >> Graph::b;

  arcs = vector<vector<Arc *>>(n+1, vector<Arc *>());
  nodesPerBlock = vector<vector<int>>(b, vector<int>());
  arcsPerBlock = vector<vector<Arc *>>(b, vector<Arc *>());

  for(i = 0; i < n; i++) {
    file >> token >> j >> aux;
    boost::split(blcs, aux, boost::is_any_of(","));
    vector<int> blocks_node = vector<int>();

    for(auto s : blcs) {
      if(s == "-1") break;
      k = stoi(s);

      blocks_node.push_back(k);
      nodesPerBlock[k].push_back(j);
    }
    nodes.push_back(make_pair(j, blocks_node));
  }

  nodes.push_back(make_pair(n, vector<int>()));

  for(k = 0; k < m; k++) {
    file >> token >> i >> j >> length >> block >> cases;
    Arc *arc = new Arc(i, j, length, block, cases);

    arcs[i].push_back(arc);
    if(block != -1) arcsPerBlock[block].push_back(arc);
  }

  for(i = 0; i < n; i++) {
    arcs[n].push_back(new Arc(n, i, 0, -1, 0));
    arcs[i].push_back(new Arc(i, n+1, 0, -1, 0));
  }

  cout << "Load graph successfully" << endl;
}

void Graph::showGraph() {
  for(int i = 0; i < n; i++)
    for(auto *arc : arcs[i])
      cout << "[" << i << ", " << arc->getD() << "] - " << arc->getLength() << ", " << arc->getBlock() << ", " << arc->getCases() << endl;
}

int Graph::getN() const {
  return n;
}

int Graph::getM() const {
  return m;
}

int Graph::getB() const {
  return b;
}

int Graph::getRoot() const {
  return n+1;
}
