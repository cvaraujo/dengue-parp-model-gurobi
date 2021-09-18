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
#include <gurobi_c++.h>
#include "boost/algorithm/string.hpp"

#include <lemon/list_graph.h>
#include <lemon/gomory_hu.h>
#include <lemon/euler.h>
#include <lemon/dijkstra.h>
#include <lemon/preflow.h>
#include <lemon/hao_orlin.h>

using namespace lemon;
using namespace std;
using namespace boost;


#define EPSILON 0.001

#endif //DPARP_INCLUDE_H
