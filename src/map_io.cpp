#include "vec_map.h"
#include "help.h"
#include <geodesy/utm.h>

void VecMap::resolve_points()
{
  global_points.resize(0);
  int points_num = raw_map.points.size();
  for(int i = 0; i < points_num; i++)
  {
    geographic_msgs::WayPoint way_point = raw_map.points[i];
    geodesy::UTMPoint utm_point;
    geodesy::fromMsg(way_point.position, utm_point);

    Point p(utm_point.easting - 605010, utm_point.northing - 4085620);

    global_points.emplace_back(way_point.id.uuid, p, way_point.props);
  
    auto iter = find_if(way_point.props.begin(), way_point.props.end(), props_find_key("junction"));
    if(iter != way_point.props.end())
    {
      Junction junction(osm_point(way_point.id.uuid, p, way_point.props));
      junctions.push_back(junction);
    }
  }

  std::cout<<"*****global point size: "<<global_points.size()<<std::endl;

  ROS_INFO("FINISH point resolved.");
}

void VecMap::resolve_roads()
{
  roads.resize(0);
  //根据feature构建道路
  //todo 将任务点，路口点读进Road
  int feature_number = raw_map.features.size();
  for(int j = 0; j < feature_number; j++)
  {
    //找到相关的点
    std::vector<osm_point> road_point;
    road_point.reserve(raw_map.features[j].components.size());
    road_point.resize(0);

    for(int l = 0; l < raw_map.features[j].components.size(); l++)
    {
      UUID id = raw_map.features[j].components[l].uuid;
      auto iter = find_if(global_points.begin(), global_points.end(), [id](osm_point &p) { return (p.id == id);});
      if(iter != global_points.end())
        road_point.push_back(*iter);
      else
        ROS_INFO("WRONG! Found point in edge but not in global map");
    
    }

    auto prop_iter = find_if(raw_map.features[j].props.begin(), raw_map.features[j].props.end(), props_find_key("landuse"));
    if(prop_iter != raw_map.features[j].props.end()){
      Obstacle ob;
      for(osm_point pt : road_point)
        ob.emplace_back(pt.p.x_, pt.p.y_);
      this->obs.push_back(ob);
      continue;
    }

    prop_iter = find_if(raw_map.features[j].props.begin(), raw_map.features[j].props.end(), props_find_key("in_junction"));
  
    if(raw_map.features[j].props.size() > 0)//Not labeled edge will not be resolved
    {
      Road road(road_point, raw_map.features[j].id.uuid);
      if(prop_iter != raw_map.features[j].props.end())
          road.in_junction = true;
      roads.push_back(road);
    }
  } 

  ROS_INFO("FINISH roads resolved.");
}

void VecMap::resolve_junctions()
{
  //构建junctuon
  //对每个junction，添加与其相连的junction_point
  std::vector<Junction>::iterator iter_junction;//查找路口的迭代器
  std::vector<Road> roads_contains_junctions;
  osm_point junction_point;
  std::vector<Road>::iterator iter_road;//查找道路的迭代器
  int junction_number = 0;
  int road_number = 0;
  for(iter_junction = junctions.begin(); iter_junction != junctions.end(); iter_junction = next(iter_junction)){//for each junction
    iter_road = this->roads.begin();
    while(iter_road != this->roads.end()){//for each road in junction
      iter_road = find_if(iter_road, this->roads.end(), roads_find_id(iter_junction->center_point.id));
      if(iter_road == roads.end())
        continue;
      else if(!iter_road->in_junction)
        continue;
      for(auto wp : iter_road->way_points){
        auto prop_iter = find_if(wp.props.begin(), wp.props.end(), props_find_key("highway"));
        if(prop_iter != wp.props.end()){
          if(prop_iter->value == "junction_point"){//find relative junction point
            //find junction point
            iter_junction->junction_points.push_back((wp));
            auto iter_highway_road = this->roads.begin();
            while(iter_highway_road != this->roads.end()){
              iter_highway_road = find_if(iter_highway_road, roads.end(), roads_find_id(wp.id));
              if(iter_highway_road == this->roads.end()){
                ROS_WARN("invalid junction point which is not linked to any road");
              }
              else if(iter_highway_road->in_junction)
                iter_highway_road = next(iter_highway_road);
              else
                break;
            }
            if(iter_highway_road == this->roads.end()){
              ROS_WARN("invalid junction point which is not linked to any highway road");
              continue;
            }
            Interface interface((*iter_highway_road), wp);
            if(iter_highway_road->way_points.front() == wp){
              iter_junction->exits.push_back(interface);
              std::cout<<"one exit solved, junction ("<<iter_junction->center_point.p.x_<<", "<<iter_junction->center_point.p.y_<<"), junction_point ("<<wp.p.x_<<", "<<wp.p.y_<<")"<<std::endl;
            }
              
            else{
              iter_junction->entries.push_back(interface);
              std::cout<<"one entry solved, junction ("<<iter_junction->center_point.p.x_<<", "<<iter_junction->center_point.p.y_<<"), junction_point ("<<wp.p.x_<<", "<<wp.p.y_<<")"<<std::endl;
            }
              
          }
        }
      }
      iter_road = next(iter_road);
    }
    iter_junction->ref_lines_generation(this->obs);
  }

  ROS_INFO("FINISH junctions resolved.");
}

void VecMap::resolve_map(ros::NodeHandle nh)
{
  ros::service::waitForService("get_geographic_map");
  ros::ServiceClient osm_client;
  osm_client = nh.serviceClient<geographic_msgs::GetGeographicMap>("get_geographic_map");
  geographic_msgs::GetGeographicMap osm_srv;
  osm_srv.request.url = this -> osm_url;
  if (!osm_client.call(osm_srv)) {
    ROS_ERROR_STREAM("error in map file " << this -> osm_url << " reading.");
  }

  raw_map = osm_srv.response.map;

  //地图数据元素
  //获取全局点列表
  resolve_points();

  //根据feature构建道路
  //todo 将任务点，路口点读进Road
  resolve_roads();
  
  //构建junctuon
  //对每个junction，添加与其相连的junction_point
  resolve_junctions();
}

VecMap::VecMap(std::string url, ros::NodeHandle nh) : osm_url(url) {
  resolve_map(nh);
}