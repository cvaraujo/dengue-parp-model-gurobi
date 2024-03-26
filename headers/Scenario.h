#ifndef DPARP_SCENARIO_H
#define DPARP_SCENARIO_H

#include "Include.h"

class Scenario
{

public:
    float probability = 0;
    vector<int> cases_per_block;

    Scenario(float probability, vector<int> cases_per_block);

    Scenario();

    void SetCases(int block, int cases);
};

#endif