//
// Created by carlos on 06/03/19.
//

#include "../headers/Arc.h"

Arc::Arc(int o, int d, float length, int block, int cases) {
  this->o = o;
  this->d = d;
  this->length = length;
  this->block = block;
  this->cases = cases;
}

int Arc::getO() const {
    return o;
}

int Arc::getD() const {
    return d;
}

float Arc::getLength() const {
    return length;
}

int Arc::getBlock() const {
    return block;
}

int Arc::getCases() const {
    return cases;
}
