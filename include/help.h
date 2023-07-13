#ifndef HELP_H
#define HELP_H

#include <string>
#include "base.h"

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

std::vector<TrackPoint> spline_interploration(TrackPoint p1, TrackPoint p2, double point_margin);

#endif


