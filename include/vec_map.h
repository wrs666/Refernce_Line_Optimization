#ifndef VEC_MAP_H
#define VEC_MAP_H

#include "base.h"
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

    //vector<map_module::curvepoint> Section_Interploration(TrackPoint p1, TrackPoint p2, int kdxdy[], double &layback);
    std::vector<Point> spline_interploration(TrackPoint p1, TrackPoint p2);
    std::vector<Point> refLineOptimization(Point p1, double theta1, Point p2, double theta2);
};

#endif