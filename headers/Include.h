//
// Created by carlos on 06/07/21.
//

#ifndef DPARP_INCLUDE_H
#define DPARP_INCLUDE_H

#include <iostream>
#include <vector>
#include "string"
#include <iomanip>
#include <bits/ios_base.h>
#include <algorithm>
#include <fstream>
#include <deque>
#include <gurobi_c++.h>
#include "boost/algorithm/string.hpp"
#include "boost/algorithm/string/trim.hpp"
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <lemon/dijkstra.h>
#include <lemon/euler.h>
#include <lemon/gomory_hu.h>
#include <lemon/hao_orlin.h>
#include <lemon/list_graph.h>
#include <lemon/preflow.h>

using namespace lemon;
using namespace std;
using namespace boost;

typedef adjacency_list<listS, vecS, directedS, no_property,
                       property<edge_weight_t, double>>
    graph_t;

#define EPSILON 0.001

#endif // DPARP_INCLUDE_H
