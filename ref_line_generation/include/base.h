#ifndef BASE_H
#define BASE_H

#include <vector>
#include <cmath>
#include <geographic_msgs/GetGeographicMap.h>

typedef geographic_msgs::KeyValue prop;
typedef std::vector<prop> keyvalue_arr;
typedef std::vector<std::pair<double, double>> Obstacle;
using UUID = boost::array<uint8_t, 16>;

enum PointType{
  JUNTION,
  JUNCTION_POINT,
  ROAD_POINT,
};

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
    Pose(Point p, double t): Point(p), theta_(t){}
    Pose(double x, double y, double t) {
      x_ = x;
      y_ = y;
      theta_ = t;
    }
    double theta_;

};

class TrackPoint : public Point
{
  public:
    TrackPoint() {}
    TrackPoint(Point p, double t, double c): Point(p), theta_(t), curvature_(c) {}
    TrackPoint(double x, double y, double t, double c) {
      x_ = x;
      y_ = y;
      theta_ = t;
      curvature_ = c;
    }
    double theta_;
    double curvature_;
};

typedef std::vector<TrackPoint> ref_line;

class Line
{
  public:
    Point start;
    Point end;
    double heading;
    Line() {}
    Line(Point p1, Point p2) : start(p1), end(p2)
    {
      heading = atan2(end.y_ - start.y_, end.x_ - start.x_);
    }

    Line generate_line(double offset, int direction) const //direction 1 is left and -1 is right
    {
      Point p1 = start.generate_Point(this->heading + direction * M_PI / 2, offset);
      Point p2 = end.generate_Point(this->heading + direction * M_PI / 2, offset);
      return Line(p1, p2);
    }
};

typedef struct osm_point
{
  UUID id;
  Point p;
  keyvalue_arr props;
  PointType type;

  osm_point() {}
  osm_point(UUID id_p, Point p_p, keyvalue_arr props_p): id(id_p), p(p_p), props(props_p) {
  }

  //利用std::find，需要重载 ==
  bool operator== (const osm_point &p)
  {
    return (this -> id == p.id);
  }

}osm_point;

class Road
{
  public:    
    UUID id_;
    double width = 4;
    static int marker_id;
    Line geometry;
    std::vector<osm_point> way_points;
    bool in_junction = false;

    Road() {}
    Road(std::vector<osm_point> points, UUID id) : way_points(points), id_(id) {
      this->geometry = Line(way_points.front().p, way_points.back().p);
    }
    bool operator== (const Road &b)
    {
      //Road类添加ID；以ID相同为判断标准
      return (this -> id_ == b.id_);
    }
};

typedef std::pair<Road, osm_point> Interface;

class Junction
{    
  public:
    osm_point center_point;
    std::vector<osm_point> junction_points;
    std::vector<Interface> entries;
    std::vector<Interface> exits;
    std::vector<ref_line> ref_lines;
    std::vector<std::vector<Point>> optimized_ref_keys;
    std::vector<ref_line> optimized_ref_lines;

    std::vector<std::vector<Point>> initial_key_points;

    Junction() {}

    Junction(osm_point junction)
    {
      this -> center_point = junction;
    }

    bool operator== (const Junction &j)
    {
      //Road类添加ID；以ID相同为判断标准
      return (this -> center_point.id == j.center_point.id);
    }

    void ref_lines_generation(const std::vector<Obstacle>& obs);
};

#endif