#pragma once
#include <vector>
#include <cmath>

typedef std::vector<std::pair<double, double>> Obstacle;

class Point{
  public:
    Point() {}
    Point(double x, double y): x_(x), y_(y) {}

    double x_;
    double y_;

    Point generate_Point(double heading, double offset) const
    {
      return Point(x_ + offset * cos(heading), y_ + offset * sin(heading));
    }
};

inline double distance(std::pair<double, double> p1, std::pair<double, double> p2){
  return sqrt(std::pow(p1.first - p2.first, 2) + std::pow(p1.second - p2.second, 2));
}

inline double distance(const Point& p1, const Point& p2){
  return distance(std::pair<double, double>(p1.x_, p1.y_), std::pair<double, double>(p2.x_, p2.y_));
}

class Pose : public Point{
  public:
    Pose() {}
    Pose(Point p): Point(p) {}
    Pose(Point p, double t): Point(p), theta_(t){}
    Pose(double x, double y, double t) {
      x_ = x;
      y_ = y;
      theta_ = t;
    }

    double theta_;

};

class TrackPoint : public Pose
{
  public:
    TrackPoint() {}
    TrackPoint(Pose p): Pose(p) {}
    TrackPoint(Pose p, double kappa): Pose(p), curvature_(kappa) {}
    TrackPoint(Point p, double t, double c){
      x_ = p.x_;
      y_ = p.y_;
      theta_ = t;
      curvature_ = c;
    }
    TrackPoint(double x, double y, double t, double c) {
      x_ = x;
      y_ = y;
      theta_ = t;
      curvature_ = c;
    }

    double curvature_;
};

typedef std::vector<TrackPoint> RefLine;

// std::vector<Pose> quinticInterploration(Pose p1, Pose p2, double point_margin);