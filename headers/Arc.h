//
// Created by Carlos on 06/07/2021.
//

#ifndef DPARP_ARC_H
#define DPARP_ARC_H

class Arc
{

private:
  int o, d, block, cases;
  float length;

public:
  Arc(int o, int d, float length, int block, int cases);

  int getO() const;

  int getD() const;

  int getBlock() const;

  int getCases() const;

  float getLength() const;

  void setCases(int cases);

  void setBlock(int block);
};

#endif // MRP_ARC_H
