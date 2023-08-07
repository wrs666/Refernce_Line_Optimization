#pragma once

#include "map_base.h"
#include <string>
#include <ros/ros.h>

class VecMap
{
  public:
    geographic_msgs::GeographicMap raw_map;
    std::string osm_url;
    std::vector<osm_point> global_points;
    std::vector<Road> roads;
    std::vector<Junction> junctions;
    std::vector<Obstacle> obs; 

    VecMap() {}
    VecMap(std::string url, ros::NodeHandle nh);

    //assign values to global points and mission points,
    //construct junctions;
    void resolve_points();
    //assign values to roads
    void resolve_roads();
    void resolve_junctions();
    void resolve_map(ros::NodeHandle nh);
};