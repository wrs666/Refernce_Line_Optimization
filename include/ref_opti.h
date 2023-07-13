#ifndef REF_OPTI_H
#define REF_OPTI_H

#include "base.h"

class RefOptimization{
  private:
    Pose start;
    Pose end;
    Obstacle obs;
  public:
    RefOptimization(Point p1, double theta1, Point p2, double theta2, Obstacle o);
};

#endif