#include "map_vis.h"
#include <type_cast.h>
#include <visualization_msgs/MarkerArray.h>
#include <string>


MapVis::MapVis(const VecMap& vec_map, ros::NodeHandle nh){
  // ros::Publisher roads_vis_pub = nh.advertise<visualization_msgs::MarkerArray>("/roads_vis", 1);
  // ros::Publisher obs_vis_pub = nh.advertise<visualization_msgs::MarkerArray>("/obs_vis", 1);
  // visualization_msgs::MarkerArray roads_vis;
  // visualization_msgs::MarkerArray obs_vis;
  
  //fill road_vis
  int road_number = 0;
  roads_vis.markers.resize(0);
  for(const Road& road : vec_map.roads){
    if(road.in_junction)
      continue;
    std::cout<<"fill road vis"<<std::endl;
    road_number++;
    visualization_msgs::Marker single_road_vis;
    single_road_vis.header.frame_id = "world";
    single_road_vis.header.stamp = ros::Time::now();
    single_road_vis.header.seq = road_number;
    single_road_vis.id = road_number;
    single_road_vis.ns = "roads";
    // single_road_vis.lifetime = ros::Duration();
    single_road_vis.action = visualization_msgs::Marker::ADD;
    single_road_vis.pose.orientation.w = 1.0;
    single_road_vis.type = visualization_msgs::Marker::LINE_LIST;
    single_road_vis.scale.x = 0.15;

    single_road_vis.id = road_number;
    single_road_vis.color.r = 1.0f;
    single_road_vis.color.g = 1.0f;
    single_road_vis.color.b = 1.0f;
    single_road_vis.color.a = 1.0;

    single_road_vis.points.push_back(getGeoPoint(road.geometry.generate_line(road.width / 2, 1).start));
    single_road_vis.points.push_back(getGeoPoint(road.geometry.generate_line(road.width / 2, 1).end));
    single_road_vis.points.push_back(getGeoPoint(road.geometry.generate_line(road.width / 2, -1).start));
    single_road_vis.points.push_back(getGeoPoint(road.geometry.generate_line(road.width / 2, -1).end));

    roads_vis.markers.push_back(single_road_vis);
  }

  //fill obs_vis
  int obs_number = 0;
  for(const Obstacle& ob : vec_map.obs){
    std::cout<<"fill obs vis"<<std::endl;
    obs_number++;
    visualization_msgs::Marker ob_vis;
    ob_vis.header.frame_id = "world";
    ob_vis.header.seq = obs_number;
    ob_vis.header.stamp = ros::Time::now();
    ob_vis.ns = "obs";
    ob_vis.id = obs_number;
    ob_vis.action = visualization_msgs::Marker::ADD;
    ob_vis.pose.orientation.w = 1.0;
    ob_vis.type = visualization_msgs::Marker::LINE_STRIP;
    ob_vis.scale.x = 0.15;

    ob_vis.color.r = 1.0f;
    ob_vis.color.a = 1.0;

    for(auto pos : ob)
      ob_vis.points.push_back(getGeoPoint(pos));
    
    obs_vis.markers.push_back(ob_vis);
  }

  
  
  for(const Junction& junc : vec_map.junctions){
    //fill ref_lines_vis
    int ref_lines_num = 0;
    for(const ref_line& r : junc.ref_lines){
      ref_lines_num++;
      nav_msgs::Path ref_path;
      ref_path.header.frame_id = "world";
      ref_path.header.seq = ref_lines_num;
      ref_path.header.stamp = ros::Time::now();
      ref_path.poses.resize(r.size());
      for(int n = 0; n < r.size(); n++){
        ref_path.poses.at(n).header.frame_id = "world";
        ref_path.poses.at(n).pose = track_to_pose(r.at(n));
      }
      ref_lines_vis.push_back(ref_path);
    }

    // fill initial key points
    int init_key_sets_num = 0;
    for(const std::vector<Point>& init_keys : junc.initial_key_points){
      init_key_sets_num++;
      std::cout<<"fill initial key points vis"<<std::endl;
      visualization_msgs::Marker init_keys_vis;
      init_keys_vis.header.frame_id = "world";
      init_keys_vis.header.seq = ref_lines_num;
      init_keys_vis.header.stamp = ros::Time::now();
      init_keys_vis.ns = "init_keys";
      init_keys_vis.id = ref_lines_num;
      init_keys_vis.action = visualization_msgs::Marker::ADD;
      init_keys_vis.pose.orientation.w = 1.0;
      init_keys_vis.type = visualization_msgs::Marker::POINTS;
      init_keys_vis.scale.x = 0.5;
      init_keys_vis.scale.y = 0.5;

      init_keys_vis.color.r = 0.0f;
      init_keys_vis.color.g = 0.0f;
      init_keys_vis.color.b = 0.0f;
      init_keys_vis.color.a = 1.0;

      for(const Point& p : init_keys)
        init_keys_vis.points.push_back(getGeoPoint(p));
      
      this->initial_keys_vis.markers.push_back(init_keys_vis);
    }

    // fill optimized key points
    int opti_key_sets_num = 0;
    for(const std::vector<Point>& opti_keys : junc.optimized_ref_keys){
      opti_key_sets_num++;
      std::cout<<"fill optimized key points vis"<<std::endl;
      visualization_msgs::Marker opti_keys_vis;
      opti_keys_vis.header.frame_id = "world";
      opti_keys_vis.header.seq = ref_lines_num;
      opti_keys_vis.header.stamp = ros::Time::now();
      opti_keys_vis.ns = "opti_keys";
      opti_keys_vis.id = ref_lines_num;
      opti_keys_vis.action = visualization_msgs::Marker::ADD;
      opti_keys_vis.pose.orientation.w = 1.0;
      opti_keys_vis.type = visualization_msgs::Marker::POINTS;
      opti_keys_vis.scale.x = 0.5;
      opti_keys_vis.scale.y = 0.5;

      opti_keys_vis.color.r = 1.0f;
      opti_keys_vis.color.g = 1.0f;
      opti_keys_vis.color.b = 0.0f;
      opti_keys_vis.color.a = 1.0;

      for(const Point& p : opti_keys)
        opti_keys_vis.points.push_back(getGeoPoint(p));
      
      this->optimized_keys_vis.markers.push_back(opti_keys_vis);
    }

    int opti_ref_lines_num = 0;
    for(const std::vector<TrackPoint>& opti_ref : junc.optimized_ref_lines){
      opti_ref_lines_num++;
      std::cout<<"fill optimized ref lines vis"<<std::endl;
      visualization_msgs::Marker opti_ref_vis;
      opti_ref_vis.header.frame_id = "world";
      opti_ref_vis.header.seq = ref_lines_num;
      opti_ref_vis.header.stamp = ros::Time::now();
      opti_ref_vis.ns = "opti_refs";
      opti_ref_vis.id = ref_lines_num;
      opti_ref_vis.action = visualization_msgs::Marker::ADD;
      opti_ref_vis.pose.orientation.w = 1.0;
      opti_ref_vis.type = visualization_msgs::Marker::POINTS;
      opti_ref_vis.scale.x = 0.2;
      opti_ref_vis.scale.y = 0.2;

      opti_ref_vis.color.r = 0.0f;
      opti_ref_vis.color.g = 1.0f;
      opti_ref_vis.color.b = 0.0f;
      opti_ref_vis.color.a = 1.0;

      for(const TrackPoint& p : opti_ref)
        opti_ref_vis.points.push_back(getGeoPoint(p));
      
      this->optimized_refs_vis.markers.push_back(opti_ref_vis);
    }
  }


  //init publisher
  roads_vis_pub = nh.advertise<visualization_msgs::MarkerArray>("/roads_vis", 1);
  obs_vis_pub = nh.advertise<visualization_msgs::MarkerArray>("/obs_vis", 1);
  this->optimized_refs_vis_pub = nh.advertise<visualization_msgs::MarkerArray>("/optimized_reference_lines", 1);
  this->optimized_keys_vis_pub = nh.advertise<visualization_msgs::MarkerArray>("/optimized_key_points", 1);
  this->initial_keys_vis_pub = nh.advertise<visualization_msgs::MarkerArray>("/initial_key_points", 1);
  ref_lines_vis_pubs.resize(ref_lines_vis.size());
  for(int i = 0; i < ref_lines_vis_pubs.size(); i++)
    ref_lines_vis_pubs.at(i) = nh.advertise<nav_msgs::Path>("/ref_lines_vis" + std::to_string(i), 1);
}

void MapVis::show(){
  roads_vis_pub.publish(roads_vis);
  obs_vis_pub.publish(obs_vis);
  optimized_refs_vis_pub.publish(optimized_refs_vis);
  optimized_keys_vis_pub.publish(optimized_keys_vis);
  initial_keys_vis_pub.publish(initial_keys_vis);
  for(int i = 0; i < ref_lines_vis.size(); i++)
    ref_lines_vis_pubs.at(i).publish(ref_lines_vis.at(i));
}