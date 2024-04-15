#ifndef DPARP_WARMSTART_H
#define DPARP_WARMSTART_H

#include "../headers/Include.h"
#include "../headers/Graph.h"

class WarmStart
{

private:
    int of;
    float default_vel = 40, spraying_vel = 10, insecticide_ml_min = 70;
    vector<pair<int, int>> solution;

    static void blocks_profit(Graph *graph, vector<pair<int, float>> &profits, double alpha);

    static bool attend_max_blocks(Graph *graph, int j, float &available_time,
                                  float max_time, float &available_insecticide, float max_insecticide,
                                  vector<pair<int, float>> profits, float &of,
                                  vector<float> &node_profit, vector<bool> &serviced, vector<pair<int, int>> &y,
                                  float spraying_vel, float insecticide_ml_min);

    static int bfs_first_profit(Graph *graph, int i, float &available_time,
                                vector<vector<bool>> &used_arcs, vector<float> &node_profit,
                                vector<pair<int, int>> &dfs_arcs, float default_vel);

public:
    static double compute_solution(Graph *graph, int max_time, float max_insecticide, vector<pair<int, int>> &x, vector<pair<int, int>> &y,
                                   float default_vel, float spraying_vel, float insecticide_ml_min, double alpha);
};

#endif // MRP_WARMSTART_H