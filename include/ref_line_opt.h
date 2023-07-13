#pragma once
#include "base.h"
#include <cppad/ipopt/solve.hpp>

using CppAD::AD;

class RefLineOPT ;

struct KeyPoint{
  std::pair<double, double> left_most;
  std::pair<double, double> right_most;
  std::pair<double, double> fore_most;
  std::pair<double, double> rear_most;
};

class FG_eval {
  public:
    typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
    FG_eval(){}
    FG_eval(RefLineOPT *opt);
    void operator()(ADvector& fg, const ADvector& x);

  private:
    // std::vector<std::pair<double, double>> free_space_curve;
    // std::vector<double> way_points_heading;
    // std::vector<std::pair<double, double>> free_space_obs;
    RefLineOPT* opt_;
    std::vector<std::vector<KeyPoint>> obs_key_points;
};


class RefLineOPT {
  private:
    std::vector<double> points_x;
    std::vector<double> points_y;
    std::vector<double> init_heading;
    std::vector<Obstacle> obs_;
    FG_eval fg_eval;

  public:
    friend class FG_eval;
    RefLineOPT() {}
    RefLineOPT(std::vector<Point> points, const std::vector<Obstacle>& obs);
    RefLineOPT(std::vector<TrackPoint> track_points, const std::vector<Obstacle>& obs);
    void obstacleFilter(const std::vector<Obstacle>& obs);
    bool solve();
    void correctHeading();
    std::vector<Point> getPoints()
    {
      std::vector<Point> points;
      for(int i = 0; i < points_x.size(); i++){
        points.emplace_back(points_x.at(i), points_y.at(i));      
      }
      return points;
    }
};