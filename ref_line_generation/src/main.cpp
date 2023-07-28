#include "map_vis.h"

const std::string map_url = "package://ref_line_generation/map_data/junction.osm";

int main(int argc, char** argv){
  ros::init(argc, argv, "curve_optimization_node");
  ros::NodeHandle nh;

  VecMap vec_map(map_url, nh);
  
  MapVis map_vis(vec_map, nh);

  ros::Rate loop_rate(1);
    
  while(ros::ok())
  {
    map_vis.show();
    loop_rate.sleep();
    ros::spinOnce();
  }

  return 0;
}