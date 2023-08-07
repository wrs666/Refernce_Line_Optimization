#pragma once

#include "base.h"
#include <geographic_msgs/GetGeographicMap.h>
#include <array>
#include <string>

typedef geographic_msgs::KeyValue prop;
typedef std::vector<prop> keyvalue_arr;
using UUID = boost::array<uint8_t, 16>;

enum PointType{
  JUNTION,
  JUNCTION_POINT,
  ROAD_POINT,
};

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
    std::vector<RefLine> ref_lines;
    std::vector<std::vector<Point>> optimized_ref_keys;
    std::vector<RefLine> optimized_ref_lines;

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


class props_find_key
{
  public:
    props_find_key(const std::string& s): key(s) {}
    bool operator () (const prop &keyvalue) const
    {
      return(keyvalue.key == key);
    }

    std::string key;

};

class roads_find_id
{
  public:
    roads_find_id(UUID arr): id(arr) {}
    bool operator () (Road &road)
    {
      auto iter = find_if(road.way_points.begin(), road.way_points.end(), [this](osm_point &p) {return (p.id == this -> id);});
      return (iter != road.way_points.end());
    }
    UUID id;
};

class junction_find_point
{
  public:
    junction_find_point(osm_point junction_point): point(junction_point) {}
    
    bool operator () (Junction &junction)
    {
      auto iter = find(junction.junction_points.begin(), junction.junction_points.end(), point);//避免在程序逻辑出错时，反复push同一个元素，导致占用内存很大
      return (iter != junction.junction_points.end());
    }
    
    osm_point point;
};


class points_find_keyvalue
{
  public:
    points_find_keyvalue(std::string k, std::string v): key(k), value(v) {}
    bool operator () (osm_point &point)
    {
      auto iter = find_if(point.props.begin(), point.props.end(), [this](prop &p) {return (p.key == this -> key && p.value == this -> value);});
      return (iter != point.props.end());
    }
    std::string key;
    std::string value;
};