#ifndef MAP_VIS_H
#define MAP_VIS_H

#include "vec_map.h"
#include <visualization_msgs/MarkerArray.h>
#include <nav_msgs/Path.h>

class MapVis {
  private:
    visualization_msgs::MarkerArray roads_vis;
    visualization_msgs::MarkerArray obs_vis;
    std::vector<nav_msgs::Path> ref_lines_vis;
    visualization_msgs::MarkerArray optimized_refs_vis;

    ros::Publisher roads_vis_pub;
    ros::Publisher obs_vis_pub;
    std::vector<ros::Publisher> ref_lines_vis_pubs;
    ros::Publisher optimized_refs_vis_pub;

  public:
    MapVis(const VecMap& map, ros::NodeHandle nh);
    void show();
};

#endif