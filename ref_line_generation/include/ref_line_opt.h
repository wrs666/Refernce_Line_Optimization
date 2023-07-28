#pragma once
#include "base.h"
#include <cppad/ipopt/solve.hpp>
#include <Eigen/Core>
#include <Eigen/Dense>

using CppAD::AD;

class RefLineOPT ;

struct KeyPoint{
  std::pair<double, double> left_most;
  std::pair<double, double> right_most;
  std::pair<double, double> fore_most;
  std::pair<double, double> rear_most;
};

class CubicSpline{
  private:
    std::vector<std::array<double, 4>> a;
    std::vector<std::array<double, 4>> b;
    std::vector<double> s;
  public:
    CubicSpline() {}
    CubicSpline(const std::vector<double>& x, const std::vector<double>& y);
    std::vector<double> getKinkHeading();
    double getKinkCurvature(int index);
    std::vector<TrackPoint> sampling(double point_margin);
};

class FG_eval {
  public:
    std::vector<double> key_x;
    std::vector<double> key_y;
    std::vector<double> initial_key_x;
    std::vector<double> initial_key_y;
    std::vector<double> key_theta;
    double start_heading;
    double end_heading;
    std::vector<std::vector<KeyPoint>> obs_key_points;
    std::vector<std::vector<int>> effective_obs_index;

    typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
    FG_eval(){}
    FG_eval(RefLineOPT *opt);
    void operator()(ADvector& fg, const ADvector& x);
    void refreshObsKey();
    void samplingSolve(double max_offset, double sampling_interval, CubicSpline& csp);
  private:
    RefLineOPT* opt_;
};


class RefLineOPT {
  private:
    std::vector<double> points_x;
    std::vector<double> points_y;
    std::vector<double> init_heading;

    std::vector<double> key_points_x;
    std::vector<double> key_points_y;

    std::vector<Obstacle> obs_;
    FG_eval fg_eval;

  public:
    friend class FG_eval;

    RefLineOPT() {}
    RefLineOPT(const std::vector<Point>& points, const std::vector<Obstacle>& obs);
    RefLineOPT(const std::vector<TrackPoint>& track_points, const std::vector<Obstacle>& obs);
    void obstacleFilter(const std::vector<Obstacle>& obs);
    bool solve(CubicSpline& csp);
    bool ssolve(CubicSpline& csp);
    void correctHeading();
    std::vector<Point> getPoints(){
      std::vector<Point> points;
      for(int i = 0; i < points_x.size(); i++){
        points.emplace_back(points_x.at(i), points_y.at(i));      
      }
      return points;
    }
    void refreshKeyPoints(double dist_thre, std::vector<double>& key_x, std::vector<double>& key_y, const std::vector<double>& x, const std::vector<double>& y);
    std::vector<Point> getKeyPoints(){
      std::vector<Point> points;
      for(int i = 0; i < this->fg_eval.key_x.size(); i++){
        points.emplace_back(this->fg_eval.key_x.at(i), this->fg_eval.key_y.at(i));      
      }
      return points;
    }
    // std::vector<TrackPoint> rawInterploration(double point_margin) const;
    // std::vector<TrackPoint> keyInterploration(double point_margin) const;
};