#pragma once
#include "base.h"
#include <memory>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <chrono>

class CubicSpline{
  private:
    std::vector<std::array<double, 4>> a;
    std::vector<std::array<double, 4>> b;
    std::vector<double> s;
  public:
    CubicSpline() {}
    CubicSpline(const std::vector<double>& x, const std::vector<double>& y);
    double getKinkHeading(int index){
      return atan2(b.at(index).at(1) + b.at(index).at(2) * s.at(index) * 2 + b.at(index).at(3) * s.at(index) * s.at(index) * 3, 
                   a.at(index).at(1) + a.at(index).at(2) * s.at(index) * 2 + a.at(index).at(3) * s.at(index) * s.at(index) * 3);
    }
    double getKinkCurvature(int index){
      double x1 = a.at(index).at(1) + a.at(index).at(2) * s.at(index) * 2 + a.at(index).at(3) * s.at(index) * s.at(index) * 3;
      double x2 = a.at(index).at(2) * 2 + a.at(index).at(3) * s.at(index) * 6;
      double y1 = b.at(index).at(1) + b.at(index).at(2) * s.at(index) * 2 + b.at(index).at(3) * s.at(index) * s.at(index) * 3;
      double y2 = b.at(index).at(2) * 2 + b.at(index).at(3) * s.at(index) * 6;
      return ((x1 * y2 - x2 * y1) / std::pow((x1 * x1 + y1 * y1), 1.5));
    }
    std::shared_ptr<RefLine> sampling(int start_index, int end_index, double point_margin);
    // std::vector<Point> samplingPoint(double point_margin);
};


class RefLineOPT {
  private:
    std::shared_ptr<RefLine> initial_path;
    std::vector<Obstacle> obstacles;
    std::vector<double> key_x;
    std::vector<double> key_y;
    std::vector<double> key_theta;
    std::vector<int> key_index;
    CubicSpline cubic_spline;

  public:
    RefLineOPT() {}
    RefLineOPT(Pose start, Pose end, const std::vector<Obstacle>& obs);
    void obstacleFilter(const std::vector<Obstacle>& obs);
    void extractKeyPoints();
    bool samplingSolve();
    std::shared_ptr<RefLine> getPoints();

    //for visualization when testing
    std::vector<Point> getKeyPoints(){
      std::vector<Point> key_points;
      int key_point_num = key_x.size();
      for(int i = 0; i < key_point_num; i++){
        key_points.emplace_back(key_x.at(i), key_y.at(i));
      }
      return key_points;
    }

    std::vector<double> init_key_x;
    std::vector<double> init_key_y;
};

std::shared_ptr<RefLine> quinticInterploration(Pose p1, Pose p2, double point_margin);