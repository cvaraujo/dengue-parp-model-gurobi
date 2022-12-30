//
// Created by carlos on 28/12/22.
//

#include "../headers/WarmStart.h"
#include "set"

void WarmStart::blocks_profit(Graph *graph, vector<pair<int, float>> &profits)
{
    int k, b = graph->getB();
    float total;
    profits = vector<pair<int, float>>();

    for (k = 0; k < b; k++)
    {
        total = 0.0;
        for (auto *arc : graph->arcsPerBlock[k])
        {
            total += arc->getCases();
        }
        profits.push_back(make_pair(k, total));
    }
}

bool WarmStart::attend_max_blocks(Graph *graph, int j, float &available_time,
                                  float max_time, float &available_insecticide, float max_insecticide,
                                  vector<pair<int, float>> profits, float &of,
                                  vector<float> &node_profit, vector<bool> &serviced,
                                  vector<pair<int, int>> &y)
{
    bool attend = false;
    for (auto block : graph->nodes[j].second)
    {
        if (serviced[block])
            continue;

        float time = available_time - graph->timeBlock(166.7, block);
        float insecticide = available_insecticide - graph->inseticideBlock(70, block);

        if (time >= 0 && insecticide >= 0)
        {
            // Mark the block as served
            y.push_back(make_pair(j, block));
            serviced[block] = true;

            of += profits[block].second;

            // Change the available amount of resources
            available_insecticide = insecticide, available_time = time;
            attend = true;

            // Remove profit from other nodes of the block
            for (auto k : graph->nodesPerBlock[block])
            {
                node_profit[k] -= profits[block].second;
            }
        }
        else
            break;
    }

    return attend;
}

int WarmStart::bfs_first_profit(Graph *graph, int i, float &available_time,
                                vector<vector<bool>> &used_arcs, vector<float> &node_profit,
                                vector<pair<int, int>> &dfs_arcs)
{
    int s, j, n = graph->getN();
    vector<bool> visited = vector<bool>(n, false);
    vector<int> stack = vector<int>(), pred = vector<int>(n, -1);
    vector<float> dist = vector<float>(n, 0);
    dfs_arcs = vector<pair<int, int>>();

    stack.push_back(i);
    visited[i] = true, pred[i] = i;
    int next_node = -1;

    while (!stack.empty())
    {
        s = stack.back();
        stack.pop_back();

        if (!visited[s])
            visited[s] = true;

        for (auto *arc : graph->arcs[s])
        {
            j = arc->getD();
            if (j >= n)
                continue;

            if (!used_arcs[s][j] && !visited[j])
            {
                stack.push_back(j);
                pred[j] = s, dist[j] = dist[s] + arc->getLength();

                if (node_profit[j] > 0 && available_time > graph->timeArc(dist[j], 333.3))
                {
                    next_node = j;
                    stack.clear();
                    break;
                }
            }
        }
    }

    if (next_node != -1)
    {
        s = next_node;
        while (s != i)
        {
            dfs_arcs.push_back(make_pair(pred[s], s));
            used_arcs[pred[s]][s] = true;
            s = pred[s];
        }
        available_time -= graph->timeArc(dist[next_node], 333.3);
        reverse(dfs_arcs.begin(), dfs_arcs.end());
    }
    return next_node;
}

double WarmStart::compute_solution(Graph *graph, int max_time, float max_insecticide, vector<pair<int, int>> &x, vector<pair<int, int>> &y)
{
    // Initial infos
    int i, j, n = graph->getN();
    float available_time = max_time;
    float available_insecticide = max_insecticide;
    float best_profit, of = 0.0;
    Arc *aux_arc;
    x = vector<pair<int, int>>();
    y = vector<pair<int, int>>();
    vector<bool> serviced = vector<bool>(graph->getB(), false);

    // Get the blocks profit.
    vector<pair<int, float>> profits;
    WarmStart::blocks_profit(graph, profits);

    // Compute node profit
    vector<float> node_profit = vector<float>(n + 1, 0);
    vector<int> degree = vector<int>(n + 1, 0);
    vector<vector<bool>> used_arcs = vector<vector<bool>>(n + 1, vector<bool>(n + 1, false));

    for (auto node_pair : graph->nodes)
    {
        i = node_pair.first;
        for (auto block : node_pair.second)
            node_profit[i] += profits[block].second;
    }

    // Create greedy route
    i = n;
    while (available_time > 0 && available_insecticide > 0)
    {
        best_profit = 0;
        for (auto *arc : graph->arcs[i])
        {
            if (node_profit[arc->getD()] > best_profit)
            {
                best_profit = node_profit[arc->getD()], aux_arc = arc;
            }
        }

        if (best_profit <= 0)
        {

            vector<pair<int, int>> bfs_arcs;

            int next_node = WarmStart::bfs_first_profit(graph, i, available_time, used_arcs, node_profit, bfs_arcs);

            if (next_node != -1)
            {

                bool attend = WarmStart::attend_max_blocks(graph, next_node, available_time, max_time,
                                                           available_insecticide,
                                                           max_insecticide, profits,
                                                           of, node_profit, serviced, y);

                if (attend)
                    x.insert(x.end(), bfs_arcs.begin(), bfs_arcs.end());
                if (node_profit[next_node] > 0)
                    break;
                else
                    node_profit[next_node] = 0, i = next_node;
            }
            else
                break;
        }
        else
        {
            j = aux_arc->getD();

            // Check the use of the arc
            available_time -= graph->timeArc(aux_arc->getLength(), 333.3);

            if (available_time <= 0)
                break;

            bool attend = WarmStart::attend_max_blocks(graph, j, available_time, max_time,
                                                       available_insecticide,
                                                       max_insecticide, profits,
                                                       of, node_profit, serviced, y);
            if (attend)
            {
                x.push_back(make_pair(i, j));
                used_arcs[i][j] = true;
            }
            if (node_profit[j] > 0)
                break;
            else
                node_profit[j] = 0, i = j;
        }
    }

    // Add the return to the dummy depot
    auto last_arc = x.back();
    x.push_back(make_pair(last_arc.second, n));

    return of;
}