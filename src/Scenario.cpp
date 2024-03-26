#include "../headers/Scenario.h"

Scenario::Scenario(float probability, vector<int> cases_per_block)
{
    this->probability = probability;
    this->cases_per_block = cases_per_block;
}

Scenario::Scenario()
{
    this->probability = 0;
}

void Scenario::SetCases(int block, int cases)
{
    this->cases_per_block[block] = cases;
}