#include "type_cast.h"

geometry_msgs::Point getGeoPoint(const Point& p){
    geometry_msgs::Point gp;
    gp.x = p.x_;
    gp.y = p.y_;
    gp.z = 0;
    return gp;
}

geometry_msgs::Point getGeoPoint(const std::pair<double, double>& pos){
    geometry_msgs::Point gp;
    gp.x = pos.first;
    gp.y = pos.second;
    gp.z = 0;
    return gp;
}

geometry_msgs::Pose track_to_pose(const TrackPoint& tp){
  geometry_msgs::Pose pose;
  pose.position.x = tp.x_;
  pose.position.y = tp.y_;
  pose.position.z = 0;
  pose.orientation.x = 0;
  pose.orientation.y = 0;
  pose.orientation.z = sin(tp.theta_/2);
  pose.orientation.w = cos(tp.theta_/2);
  return pose;
}