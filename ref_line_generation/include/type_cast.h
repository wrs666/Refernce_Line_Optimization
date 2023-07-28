#ifndef TYPE_CAST_H
#define TYPE_CAST_H

#include <base.h>
#include <geometry_msgs/Point.h>
#include <nav_msgs/Path.h>

geometry_msgs::Point getGeoPoint(const Point& p);
geometry_msgs::Point getGeoPoint(const std::pair<double, double>& pos);
geometry_msgs::Pose track_to_pose(const TrackPoint& tp);

#endif