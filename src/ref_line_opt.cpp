#include "ref_line_opt.h"
#include <cassert>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <float.h>

//TODO: extract params

//sampling param
const double initial_sampling_margin = 1.0;
const double result_sampling_margin = 0.1;

//distance for avoiding collision params
const double lateral_threshold = 2;
const double longitude_threshold = 5;

//sampling optimization params
//  cost weight param
const double smooth_weight = 1e+3; //cost value about 10
const double boundary_orientation_weight = 5e+3;
const double obs_weight = 1;
const double length_weight = 3e+2;
const double exponential_scale = 10;
const double distance_threshold = 2.5;
//  sampling params
double sampling_solve_offset = 30.0; 
double sampling_solve_interval = 2.0;
const double max_key_points_interval = 25.0;

//extract polyline key points param
const double threshold_for_extraction = 0.5;


CubicSpline::CubicSpline(const std::vector<double>& x, const std::vector<double>& y){//toImprove : combine loop
  int n = x.size();
  this->a.resize(n - 1);
  this->b.resize(n - 1);

  int N = 8 * (n - 1); 

  this->s = std::vector<double>(n, 0);
  
  Eigen::MatrixXd A(N, N);
  A.setZero();
  Eigen::MatrixXd X(N, 1);
  Eigen::MatrixXd B(N, 1);
  B.setZero();

  //fill B and s
  B(0, 0) = x.front();
  B(1, 0) = y.front();
  A(0, 0) = 1;
  A(1, 4) = 1;

  for(int i = 1; i < n - 1; i++){
    //fill B
    B(4 * i - 2, 0) = x.at(i);
    B(4 * i - 1, 0) = y.at(i);
    B(4 * i, 0) = x.at(i);
    B(4 * i + 1, 0) = y.at(i);

    //fill s
    s.at(i) = s.at(i - 1) + sqrt(std::pow(x.at(i) - x.at(i - 1), 2) + std::pow(y.at(i) - y.at(i - 1), 2));

    //fill A
    double x_start = 4 * i - 2;
    double y_start = 8 * (i - 1);
    for(int j = 0; j < 4; j++){
      A(x_start + j, y_start + j * 4) = 1;
      A(x_start + j, y_start + j * 4 + 1) = s.at(i);
      A(x_start + j, y_start + j * 4 + 2) = s.at(i) * s.at(i);
      A(x_start + j, y_start + j * 4 + 3) = s.at(i) * s.at(i) * s.at(i);
    }
    A(4 * n - 8 + 4 * i, y_start + 1) = -1;
    A(4 * n - 8 + 4 * i, y_start + 2) = -2 * s.at(i);
    A(4 * n - 8 + 4 * i, y_start + 3) = -3 * s.at(i) * s.at(i);
    A(4 * n - 8 + 4 * i, y_start + 9) = 1;
    A(4 * n - 8 + 4 * i, y_start + 10) = 2 * s.at(i);
    A(4 * n - 8 + 4 * i, y_start + 11) = 3 * s.at(i) * s.at(i);
    A(4 * n - 7 + 4 * i, y_start + 2) = -2;
    A(4 * n - 7 + 4 * i, y_start + 3) = -6 * s.at(i);
    A(4 * n - 7 + 4 * i, y_start + 10) = 2;
    A(4 * n - 7 + 4 * i, y_start + 11) = 6 * s.at(i);
    A(4 * n - 6 + 4 * i, y_start + 5) = -1;
    A(4 * n - 6 + 4 * i, y_start + 6) = -2 * s.at(i);
    A(4 * n - 6 + 4 * i, y_start + 7) = -3 * s.at(i) * s.at(i);
    A(4 * n - 6 + 4 * i, y_start + 13) = 1;
    A(4 * n - 6 + 4 * i, y_start + 14) = 2 * s.at(i);
    A(4 * n - 6 + 4 * i, y_start + 15) = 3 * s.at(i) * s.at(i);
    A(4 * n - 5 + 4 * i, y_start + 6) = -2;
    A(4 * n - 5 + 4 * i, y_start + 7) = -6 * s.at(i);
    A(4 * n - 5 + 4 * i, y_start + 14) = 2;
    A(4 * n - 5 + 4 * i, y_start + 15) = 6 * s.at(i);


  }

  B(4 * n - 6, 0) = x.back();
  B(4 * n - 5, 0) = y.back();

  s.back() = s.at(n - 2) + sqrt(std::pow(x.back() - x.at(n - 2), 2) + std::pow(y.back() - y.at(n - 2), 2));

  A(4 * n - 6, 8 * (n - 2)) = 1;
  A(4 * n - 6, 8 * (n - 2) + 1) = s.back();
  A(4 * n - 6, 8 * (n - 2) + 2) = s.back() * s.back();
  A(4 * n - 6, 8 * (n - 2) + 3) = s.back() * s.back() * s.back();
  A(4 * n - 5, 8 * (n - 2) + 4) = 1;
  A(4 * n - 5, 8 * (n - 2) + 5) = s.back();
  A(4 * n - 5, 8 * (n - 2) + 6) = s.back() * s.back();
  A(4 * n - 5, 8 * (n - 2) + 7) = s.back() * s.back() * s.back();

  A(8 * n - 12, 2) = 2;
  A(8 * n - 11, 6) = 2;
  A(8 * n - 10, 8 * n - 14) = 2;
  A(8 * n - 10, 8 * n - 13) = 6 * s.back();
  A(8 * n - 9, 8 * n - 10) = 2;
  A(8 * n - 9, 8 * n - 9) = 6 * s.back();

  X = A.colPivHouseholderQr().solve(B);

  for(int i = 0; i < n - 1; i++){
    for(int j = 0; j < 4; j++){
      this->a.at(i).at(j) = X(8 * i + j, 0);
      this->b.at(i).at(j) = X(8 * i + 4 + j, 0);
    }
  }
}

std::shared_ptr<RefLine> quinticInterploration(Pose p1, Pose p2, double point_margin){
  std::shared_ptr<RefLine> curve_points = std::make_shared<RefLine>();
  curve_points->resize(0);

  double x1, y1, theta1, kappa1, x2, y2, theta2, kappa2;
  x1 = p1.x_;
  y1 = p1.y_;
  theta1 = p1.theta_;
  kappa1 = 0;
  x2 = p2.x_;
  y2 = p2.y_;
  theta2 = p2.theta_;
  kappa2 = 0;
  
  double s_distance = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));

  curve_points->reserve((int)ceil(s_distance / point_margin) + 1);

  double vx1, vx2, vy1, vy2, ax1, ax2, ay1,ay2;
  vx1 = cos(theta1);
  vx2 = cos(theta2);
  vy1 = sin(theta1);
  vy2 = sin(theta2);
  ax1 = -sin(theta1) * kappa1;
  ax2 = -sin(theta2) * kappa2;
  ay1 = cos(theta1) * kappa1;
  ay2 = cos(theta2) * kappa2;

  double a[6];
  a[0] = x1;
  a[1] = vx1;
  a[2] = ax1 / 2;
  a[3] = (-3 * ax1 + ax2) / (2 * s_distance) + (-6 * vx1 - 4 * vx2) / pow(s_distance, 2) - 10 * (x1 -x2) / pow(s_distance, 3);
  a[4] = (1.5 * ax1 - ax2) / pow(s_distance, 2) + (8 * vx1 + 7 * vx2) / pow(s_distance, 3) + 15 * (x1 - x2) / pow(s_distance, 4);
  a[5] = (-ax1 + ax2) / (2 * pow(s_distance, 3)) - 3 * (vx1 + vx2) / pow(s_distance, 4) + (-6 * x1 + 6 * x2) / pow(s_distance, 5);

  double b[6];
  b[0] = y1;
  b[1] = vy1;
  b[2] = ay1 / 2;
  b[3] = (-3 * ay1 + ay2) / (2 * s_distance) + (-6 * vy1 - 4 * vy2) / pow(s_distance, 2) - 10 * (y1 -y2) / pow(s_distance, 3);
  b[4] = (1.5 * ay1 - ay2) / pow(s_distance, 2) + (8 * vy1 + 7 * vy2) / pow(s_distance, 3) + 15 * (y1 - y2) / pow(s_distance, 4);
  b[5] = (-ay1 + ay2) / (2 * pow(s_distance, 3)) - 3 * (vy1 + vy2) / pow(s_distance, 4) + (-6 * y1 + 6 * y2) / pow(s_distance, 5); 

  curve_points->push_back(p1);

  //sampling
  double s = point_margin;
  TrackPoint p;
  double vx, vy, ax, ay, theta;

  Eigen::VectorXd coef_x(6);
  Eigen::VectorXd coef_y(6);
  Eigen::VectorXd coef_vx(6);
  Eigen::VectorXd coef_vy(6);
  Eigen::VectorXd coef_ax(6);
  Eigen::VectorXd coef_ay(6);
    
  coef_x << a[0], a[1], a[2], a[3], a[4], a[5];
  coef_y << b[0], b[1], b[2], b[3], b[4], b[5];
  coef_vx << a[1], 2 * a[2], 3 * a[3], 4 * a[4], 5 * a[5], 0;
  coef_vy << b[1], 2 * b[2], 3 * b[3], 4 * b[4], 5 * b[5], 0;
  coef_ax << 2 * a[2], 6 * a[3], 12 * a[4], 20 * a[5], 0, 0;
  coef_ay << 2 * b[2], 6 * b[3], 12 * b[4], 20 * b[5], 0, 0;

  while(s < s_distance)
  {
    Eigen::VectorXd s_power(6);
    s_power << 1, s, pow(s, 2), pow(s, 3), pow(s, 4), pow(s, 5);
    
    p.x_ = coef_x.dot(s_power);
    p.y_ = coef_y.dot(s_power);
    vx = coef_vx.dot(s_power);
    vy = coef_vy.dot(s_power);
    ax = coef_ax.dot(s_power);
    ay = coef_ay.dot(s_power);

    p.theta_ = std::atan2(vy, vx);
    p.curvature_ = (vx * ay - ax * vy) / pow(vx * vx + vy * vy, 1.5);

    curve_points->push_back(p);
    s = s + point_margin;
  }

  if(s_distance - s > point_margin / 2)
      curve_points->push_back(p2);
  // std::cout<<"Finish quintic spline interplorartion, points size is "<<(int)curve_points->size()<<std::endl;
  return curve_points;
}

std::shared_ptr<RefLine> CubicSpline::sampling(int start_index, int end_index, double point_margin){//no curvature calculation
  //sampling
  int n = this->s.size();
  std::shared_ptr<RefLine> sampling_points = std::make_shared<RefLine>();
  double station = this->s.at(start_index) + point_margin;
  int segment_index = start_index;
  while(station < s.at(end_index)){
    while(station > s.at(segment_index + 1)){      
      segment_index++;
      if(segment_index >= end_index)
          break;
    }
    sampling_points->emplace_back(a.at(segment_index).at(0) * 1 + a.at(segment_index).at(1) * station + a.at(segment_index).at(2) * station * station + a.at(segment_index).at(3) * station * station * station, 
                                  b.at(segment_index).at(0) * 1 + b.at(segment_index).at(1) * station + b.at(segment_index).at(2) * station * station + b.at(segment_index).at(3) * station * station * station, 
                                  atan2(b.at(segment_index).at(1) + b.at(segment_index).at(2) * station * 2 + b.at(segment_index).at(3) * station * station * 3, a.at(segment_index).at(1) + a.at(segment_index).at(2) * station * 2 + a.at(segment_index).at(3) * station * station * 3), 0);
    station += point_margin;
    // std::cout<<"add one point: x "<<sampling_points.back().x_<<", y "<<sampling_points.back().y_<<std::endl;
  }
  // std::cout<<"Finish cubic interploration and sampling, points size "<<(int)sampling_points->size()<<std::endl;
  return sampling_points;
}


void extractPolyKeyPoints(const std::vector<std::pair<double, double>>& polyline, double min_distance, int start, int end, std::vector<bool>& black_list ){
  if(end - start < 2)
    return;
  double max_distance = 0;
  int max_offset_index;
  double start_x = polyline.at(start).first;
  double start_y = polyline.at(start).second;
  double end_x = polyline.at(end).first;
  double end_y = polyline.at(end).second;
  for(int i = start + 1; i < end; i++){
    double current_distance = fabs((start_y - end_y) * polyline.at(i).first + (end_x - start_x) * polyline.at(i).second + start_x * end_y - end_x * start_y) / sqrt((start_y - end_y) * (start_y - end_y) + (start_x - end_x) * (start_x - end_x));
    if(current_distance > max_distance){
      max_distance = current_distance;
      max_offset_index = i;
    }
  }

  if(max_distance < min_distance){
    for(int i = start + 1; i < end; i++)
    black_list.at(i) = true;
  }
  else{
    extractPolyKeyPoints(polyline, min_distance, start, max_offset_index, black_list);
    extractPolyKeyPoints(polyline, min_distance, max_offset_index, end, black_list);
  }
}

Obstacle singleObstacleFilter(const Obstacle& ob){
  Obstacle filtered_ob;
  std::vector<bool> black_list(ob.size(), false);
  extractPolyKeyPoints(ob, threshold_for_extraction, 0, (int)ob.size() - 1, black_list);
  for(int i = 0; i < ob.size(); i++){
    if(!black_list.at(i))
      filtered_ob.push_back(ob.at(i));
  }
  return filtered_ob;
}


void RefLineOPT::obstacleFilter(const std::vector<Obstacle>& obs){    //toFix : extract key points of obstacle polyline
    double start_to_end = distance(this->initial_path->front(), this->initial_path->back());
    std::pair<double, double> start(this->initial_path->front().x_, this->initial_path->front().y_);
    std::pair<double, double> end(this->initial_path->back().x_, this->initial_path->back().y_);
    this->obstacles.resize(0);
    for(const Obstacle& ob : obs){
        for(std::pair<double, double> p : ob){
            if(distance(p, start) < start_to_end && distance(p, end) < start_to_end){
              Obstacle filtered_ob = singleObstacleFilter(ob);
              this->obstacles.push_back(filtered_ob);
              // std::cout<<"filtering obstacle, choose obstacle whose first point is ("<<filtered_ob.front().first<<", "<<filtered_ob.front().second<<")"<<std::endl;
              break;
            }
        }
    }
}

void RefLineOPT::extractKeyPoints(){
  int n = this->initial_path->size();
  this->key_y.resize(0);
  this->key_x.resize(0);
  this->key_theta.resize(0);
  key_x.push_back(this->initial_path->front().x_);
  key_y.push_back(this->initial_path->front().y_);
  key_theta.push_back(this->initial_path->front().theta_);
  key_index.resize(0);

  bool selected = false, last_selected = false;
  bool refresh = true;
  bool increasing = true;
  double distance_measurement, last_distance_measurement = 1000;

  int obs_size = this->obstacles.size();
  for(int i = 1; i < n - 1; i++){//for each way points
    distance_measurement = 1000;
    selected = false;
    refresh = false;
    // double theta;
    // std::cout<<"For this point : "<<std::endl;
    double current_theta = initial_path->at(i).theta_;
    double current_x = initial_path->at(i).x_;
    double current_y = initial_path->at(i).y_;
    for(int j = 0; j < obs_size; j++){//for each obstacle in this freesapce
      double MIN_X = 100;//TODO: extract params 
      double MAX_X = -100;
      double MIN_Y = 100;
      double MAX_Y = -100;
      // std::cout<<"  For this obstacle : "<<std::endl;
      Obstacle current_ob = this->obstacles.at(j);
      for(std::pair<double, double> p : current_ob){//toFix : the complexity can be reduced further
        double relative_x = cos(current_theta) * (p.first - current_x) + sin(current_theta) * (p.second - current_y);
        double relative_y = -sin(current_theta) * (p.first - current_x) + cos(current_theta) * (p.second - current_y);
        
        if(relative_x > MAX_X)
            MAX_X = relative_x;

        if(relative_x < MIN_X)
            MIN_X = relative_x;

        if(relative_y > MAX_Y)
            MAX_Y = relative_y;
        
        if(relative_y < MIN_Y)
            MIN_Y = relative_y;
      }

      double lateral_distance_measurement = (MIN_X + MAX_X) * (MIN_X + MAX_X) / ((MAX_X - MIN_X + 2 * longitude_threshold) * (MAX_X - MIN_X + 2 * longitude_threshold));
      double longitude_distance_measurement = (MIN_Y + MAX_Y) * (MIN_Y + MAX_Y) / ((MAX_Y - MIN_Y + 2 * lateral_threshold) * (MAX_Y - MIN_Y + 2 * lateral_threshold));
      
      if(lateral_distance_measurement < 1 && longitude_distance_measurement < 1){//toFix : check whether the uniform standard is better
          selected = true;
          refresh = !last_selected;
      }

      if(lateral_distance_measurement + longitude_distance_measurement < distance_measurement)
        distance_measurement = lateral_distance_measurement + longitude_distance_measurement;

    }

    if(selected){//no need to insert key points for uniform
      // std::cout<<i<<"th point selected, distance measurement "<<distance_measurement<<std::endl;
      if(refresh){//start
        key_index.push_back(i);
        // std::cout<<i<<"th point added to key points "<<std::endl;
      }
      else{
        if(increasing ^ (distance_measurement > last_distance_measurement)){//insert extremum
          int insert_key_num = floor(initial_sampling_margin * (i - 1 - key_index.back()) / max_key_points_interval);
          if(insert_key_num > 0){
            // std::cout<<"!!insert key points, number "<<insert_key_num<<", after "<<(int)key_index.size()<<" key points."<<std::endl;
            // std::cout<<"distance to last key point : "<<initial_sampling_margin * (i - 1 - key_index.back())<<std::endl;
            int index_increment = (i - 1 - key_index.back()) / (insert_key_num + 1);
            for(int j = 1; j <= insert_key_num; j++){
                key_index.push_back(key_index.back() + index_increment);
                // std::cout<<"  "<<key_index.back()<<"th point added to key points "<<std::endl;
            }
          }
          key_index.push_back(i - 1);
          // std::cout<<"  "<<i - 1<<"th point added to key points "<<std::endl;
        }
      }
    }
    else if(last_selected){//insert end in one section
      int insert_key_num = ceil(initial_sampling_margin * (i - 1 - key_index.back()) / max_key_points_interval) - 1;
      if(insert_key_num > 0){
        // std::cout<<"!!insert key points, number "<<insert_key_num<<", after "<<(int)key_index.size()<<"th key points."<<std::endl;
        // std::cout<<"distance to last key point : "<<initial_sampling_margin * (i - 1 - key_index.back())<<std::endl;
        int index_increment = (i - 1 - key_index.back()) / (insert_key_num + 1);
        for(int j = 1; j <= insert_key_num; j++){
          key_index.push_back(key_index.back() + index_increment);
          // std::cout<<"  "<<key_index.back()<<"th point added to key points "<<std::endl;
        }
      }

      key_index.push_back(i - 1);
      // std::cout<<"  "<<i - 1<<"th point added to key points "<<std::endl;
    }

    increasing = (distance_measurement > last_distance_measurement);

    last_distance_measurement = distance_measurement;
    last_selected = selected;
  }

  for(int index : key_index){
    key_x.push_back(initial_path->at(index).x_);
    key_y.push_back(initial_path->at(index).y_);
    key_theta.push_back(initial_path->at(index).theta_);
  }

  key_x.push_back(initial_path->back().x_);
  key_y.push_back(initial_path->back().y_);
  key_theta.push_back(initial_path->back().theta_);

  //for test
  init_key_x = key_x;
  init_key_y = key_y;

  // std::cout<<"Number of key points selected is "<<(int)this->key_x.size()<<". Number of key index is "<<(int)this->key_index.size()<<std::endl;
}

inline double pointToLineDistance(double px, double py, double p1_x, double p1_y, double p2_x, double p2_y){
  // std::cout<<"       point ("<<px<<", "<<py<<"),  line start ("<<p1_x<<", "<<p1_y<<"), end ("<<p2_x<<", "<<p2_y<<").    ";
  //i is p1, j is p2
  double ip_x = px - p1_x;
  double ip_y = py - p1_y;
  double ij_x = p2_x - p1_x;
  double ij_y = p2_y - p1_y;
  double pj_x = p2_x - px;
  double pj_y = p2_y - py;
  double r = (ip_x * ij_x + ip_y * ij_y) / (ij_x * ij_x + ij_y * ij_y);
  double current_d;
  if(r < 0)
    current_d = sqrt(ip_x * ip_x + ip_y * ip_y);
  else if(r > 1){
    current_d = sqrt(pj_x* pj_x + pj_y * pj_y);
  }
  else{
    current_d = fabs(ij_y * px - ij_x * py + p2_x * p1_y - p1_x * p2_y) / sqrt(ij_x * ij_x + ij_y * ij_y);
  }
  // std::cout<<"      distance "<<current_d<<std::endl;
  return current_d;
}

inline double lineToObstacleDistance(double p1_x, double p1_y, double p2_x, double p2_y, const Obstacle& obs){
  int n = obs.size();
  double dist = 100;

  double x_diff = p2_x - p1_x;
  double y_diff = p2_y - p1_y;

  for(int i = 0; i < n - 1; i++){
    // whether is cross
    double ix = obs.at(i).first;
    double iy = obs.at(i).second;
    double jx = obs.at(i + 1).first;
    double jy = obs.at(i + 1).second;
    double deno = x_diff * (jy - iy) - y_diff * (jx - ix);
    double mem1 = (p1_y - iy) * (jx - ix) - (p1_x - ix) * (jy - iy);
    double mem2 = (p1_y - iy) * x_diff - (p1_x - ix) * y_diff;
    double r = mem1 / deno;
    double s = mem2 / deno;
    if(r > 0 && r < 1 && s > 0 && s < 1)
      return 0;

    //distance from obstacle point to line
    double d = pointToLineDistance(ix, iy, p1_x, p1_y, p2_x, p2_y);
    if(d < dist)
      dist = d;
    d = pointToLineDistance(p1_x, p1_y, ix, iy, jx, jy);
    if(d < dist)
      dist = d;
    d = pointToLineDistance(p2_x, p2_y, ix, iy, jx, jy);
    if(d < dist)
      dist = d;
  }
  double d = pointToLineDistance(obs.back().first, obs.back().second, p1_x, p1_y, p2_x, p2_y);
  if(d < dist)
    dist = d;

  return dist;
}

inline bool pointInPolygon(std::pair<double, double> p, std::vector<std::pair<double, double>> polygon){//points of polygon arranged counterclockwise
  int n = polygon.size();
  for(int i = 0, j = n - 1; i < n - 1; j = i++){
    if((polygon.at(i).first - polygon.at(j).first) * (p.second * polygon.at(j).second) - (p.first - polygon.at(j).first) * (polygon.at(i).second - polygon.at(j).second) < 0)
      return false;
  }
  return true;
}

inline bool crossWithPolygon(double p1_x, double p1_y, double p2_x, double p2_y, std::vector<std::pair<double, double>> polygon){
  int n = polygon.size();

  double x_diff = p2_x - p1_x;
  double y_diff = p2_y - p1_y;

  for(int i = 0, j = n - 1; i < n - 1; j = i++){
    // whether is cross
    double deno = x_diff * (polygon.at(i).second - polygon.at(j).second) - y_diff * (polygon.at(i).first - polygon.at(j).first);
    double mem1 = (p1_y - polygon.at(j).second) * (polygon.at(i).first - polygon.at(j).first) - (p1_x - polygon.at(j).first) * (polygon.at(i).second - polygon.at(j).second);
    double mem2 = (p1_y - polygon.at(j).second) * x_diff - (p1_x - polygon.at(j).first) * y_diff;

    double r = mem1 / deno;
    double s = mem2 / deno;
    if(r > 0 && r < 1 && s > 0 && s < 1)
      return true;
  }
  return false;
}


bool RefLineOPT::samplingSolve(){
  std::cout<<"************Implementing Sampling Solving***********"<<std::endl;
  auto t1 = std::chrono::steady_clock::now();

  int m = 2 * (sampling_solve_offset / sampling_solve_interval) + 1;
  int n = this->key_x.size();
  // std::cout<<"Number of sampling states per layer : "<<m<<". Key points size : "<<n<<std::endl;
  std::vector<std::vector<double>> dp(n, std::vector<double>(m, 0));
  std::vector<std::unordered_map<double, double>> link(n);
  std::vector<double> last_heading(m), current_heading(m);
  std::vector<std::vector<double>> cumulative_length(n, std::vector<double>(m, DBL_MAX));

  // std::cout<<"For the end point"<<std::endl;
  //filter obstacle for this level
  //whether the first point is in the box:
      //in : choose this obstacle
      //out : each line is not cross with the box
  // std::cout<<"  Start filtering obstacles for the n - 2 level."<<std::endl; 
  std::vector<std::pair<double, double>> path_polygon;
  path_polygon.emplace_back(key_x.back(), key_y.back());
  double max_offset = sampling_solve_offset + 2.0;
  double origin_x = key_x.at(n - 2);
  double origin_y = key_y.at(n - 2);
  double tangent_theta = key_theta.at(n - 2) - M_PI / 2;
  path_polygon.emplace_back(origin_x - max_offset * cos(tangent_theta), origin_y - max_offset * sin(tangent_theta));
  path_polygon.emplace_back(origin_x + max_offset * cos(tangent_theta), origin_y + max_offset * sin(tangent_theta));
  std::vector<int> effective_obs_index {};
  for(int ob_index = 0; ob_index < this->obstacles.size(); ob_index++){
    Obstacle ob = this->obstacles.at(ob_index);
    if(pointInPolygon(ob.front(), path_polygon)){
      effective_obs_index.push_back(ob_index);
      // std::cout<<"    the first point is in the path polygon, the first point is ("<<ob.front().first<<", "<<ob.front().second<<")"<<std::endl;
    }
      
    else{
      int n = ob.size();
      for(int i = 0; i < n - 1; i++){
        if(crossWithPolygon(ob.at(i).first, ob.at(i).second, ob.at(i + 1).first, ob.at(i + 1).second, path_polygon)){
          effective_obs_index.push_back(ob_index);
          // std::cout<<"    line cross with the path polygon, ("<<ob.at(i).first<<", "<<ob.at(i).second<<")  ("<<ob.at(i + 1).first<<", "<<ob.at(i + 1).second<<")"<<std::endl;
          break;
        }
      }
    }
  }

  for(int i = 0; i < m; i++){//toFix : improve the efficiency
    double d = -sampling_solve_offset + i * sampling_solve_interval;
    // std::cout<<"  offset "<<d<<std::endl;
    double moved_x = origin_x + d * cos(tangent_theta);
    double moved_y = origin_y + d * sin(tangent_theta);
    double section_theta = atan2(key_y.back() - moved_y, key_x.back() - moved_x);
    last_heading.at(i) = section_theta;
    double ob_distance = 100;
    //filter obstacle for this level
    for(int ob_index : effective_obs_index){
      double tmp_ob_distance = lineToObstacleDistance(moved_x, moved_y, key_x.back(), key_y.back(), this->obstacles.at(ob_index));
      if(tmp_ob_distance < ob_distance)
        ob_distance = tmp_ob_distance;
    }

    double smooth_cost = boundary_orientation_weight * (section_theta - key_theta.back()) * (section_theta - key_theta.back());
    // std::cout<<"smoooth cost "<<smooth_cost<<std::endl;
    double obs_cost = obs_weight * exp((distance_threshold - ob_distance) * exponential_scale);
    // std::cout<<"obs cost "<<obs_cost<<std::endl;
    double length_cost = length_weight * sqrt((key_x.back() - moved_x) * (key_x.back() - moved_x) + (key_y.back() - moved_y) * (key_y.back() - moved_y));
    // std::cout<<"length cost "<<length_cost<<std::endl;
    dp.at(n - 2).at(i) = smooth_cost + obs_cost + length_cost;
  }


  //index : n-3 to 1
  for(int i = n - 3; i > 0; i--){//toFix : improve the efficiency
    // std::cout<<"Sampling "<<i<<"th point"<<std::endl;

    double current_origin_x = key_x.at(i);
    double current_origin_y = key_y.at(i);
    double current_tangent_theta = key_theta.at(i) - M_PI / 2;
    double last_origin_x = key_x.at(i + 1);
    double last_origin_y = key_y.at(i + 1);
    double last_tangent_theta = key_theta.at(i + 1) - M_PI / 2;

    //filter obstacles for this level
    // std::cout<<"  Start filtering obstacles for the "<<i<<" level."<<std::endl; 
    std::vector<std::pair<double, double>> current_path_polygon;
    current_path_polygon.emplace_back(last_origin_x - max_offset * cos(last_tangent_theta), last_origin_y - max_offset * sin(last_tangent_theta));
    current_path_polygon.emplace_back(current_origin_x - max_offset * cos(current_tangent_theta), current_origin_y - max_offset * sin(current_tangent_theta));
    current_path_polygon.emplace_back(current_origin_x + max_offset * cos(current_tangent_theta), current_origin_y + max_offset * sin(current_tangent_theta));
    current_path_polygon.emplace_back(last_origin_x + max_offset * cos(last_tangent_theta), last_origin_y + max_offset * sin(last_tangent_theta));
    std::vector<int> current_effective_obs_index {};
    for(int ob_index = 0; ob_index < this->obstacles.size(); ob_index++){
      Obstacle ob = this->obstacles.at(ob_index);
      int n = ob.size();
      if(pointInPolygon(ob.front(), current_path_polygon)){
        current_effective_obs_index.push_back(ob_index);
        // std::cout<<"    the first point is in the path polygon, the first point is ("<<ob.front().first<<", "<<ob.front().second<<")"<<std::endl;
      }
        
      else{
        for(int i = 0; i < n - 1; i++){
          if(crossWithPolygon(ob.at(i).first, ob.at(i).second, ob.at(i + 1).first, ob.at(i + 1).second, current_path_polygon)){
            current_effective_obs_index.push_back(ob_index);
            // std::cout<<"    line cross with the path polygon, ("<<ob.at(i).first<<", "<<ob.at(i).second<<")  ("<<ob.at(i + 1).first<<", "<<ob.at(i + 1).second<<")"<<std::endl;
            break;
          }
        }
      }
    }

    for(int j = 0; j < m; j++){//for the current level
      double d = -sampling_solve_offset + j * sampling_solve_interval;
      double moved_x = current_origin_x + d * cos(current_tangent_theta);
      double moved_y = current_origin_y + d * sin(current_tangent_theta);
      // std::cout<<"  Offset "<<d<<std::endl;

      //for each states in the last level
      double min_cost = DBL_MAX;
      for(int last_j = 0; last_j < m; last_j++){
        double last_d = -sampling_solve_offset + last_j * sampling_solve_interval;
        // std::cout<<"    last offset "<<last_d<<std::endl;
        double last_x =  last_origin_x + last_d * cos(last_tangent_theta);
        double last_y =  last_origin_y + last_d * sin(last_tangent_theta);

        double current_ob_distance = 100;
        for(int ob_index : current_effective_obs_index){
          double tmp_ob_distance = lineToObstacleDistance(moved_x, moved_y, last_x, last_y, this->obstacles.at(ob_index));

          if(tmp_ob_distance < current_ob_distance)
            current_ob_distance = tmp_ob_distance;
        }
        // std::cout<<"      ob distance "<<current_ob_distance<<std::endl;
        //get obs cost
        double current_obs_cost = obs_weight * exp((distance_threshold - current_ob_distance) * exponential_scale);
        // std::cout<<"      current obs cost "<<current_obs_cost<<std::endl;
        //get length cost
        double current_length_cost = length_weight * sqrt((last_x - moved_x) * (last_x - moved_x) + (last_y - moved_y) * (last_y - moved_y));
        // std::cout<<"      current length cost "<<current_length_cost<<std::endl;
        //get smooth cost
        double section_theta = atan2(last_y - moved_y, last_x - moved_x);
        double current_smooth_cost = smooth_weight * (last_heading.at(last_j) - section_theta) * (last_heading.at(last_j) - section_theta);
        double current_cost = current_length_cost + current_obs_cost + current_smooth_cost;
        if(current_cost + dp.at(i + 1).at(last_j) < min_cost){
          current_heading.at(j) = section_theta;
          link.at(i)[d] = last_d;
          min_cost = current_cost + dp.at(i + 1).at(last_j);
        }
      }
      dp.at(i).at(j) = min_cost;
    }
    last_heading = current_heading;
  }

  //index 1
  // std::cout<<"For the start point. "<<std::endl;    
  //filter obstacles for this level
  // std::cout<<"  Start filtering obstacles for the start point."<<std::endl; 
  std::vector<std::pair<double, double>> current_path_polygon;
  current_path_polygon.emplace_back(key_x.at(1) - (sampling_solve_offset + 2) * cos(key_theta.at(1) - M_PI / 2), key_y.at(1) - (sampling_solve_offset + 2) * sin(key_theta.at(1) - M_PI / 2));
  current_path_polygon.emplace_back(key_x.front(), key_y.front());
  current_path_polygon.emplace_back(key_x.at(1) + (sampling_solve_offset + 2) * cos(key_theta.at(1) - M_PI / 2), key_y.at(1) + (sampling_solve_offset + 2) * sin(key_theta.at(1) - M_PI / 2));
  std::vector<int> current_effective_obs_index {};
  for(int ob_index = 0; ob_index < this->obstacles.size(); ob_index++){
    Obstacle ob = this->obstacles.at(ob_index);
    int n = ob.size();
    if(pointInPolygon(ob.front(), current_path_polygon)){
      current_effective_obs_index.push_back(ob_index);
      // std::cout<<"    the first point is in the path polygon, the first point is ("<<ob.front().first<<", "<<ob.front().second<<")"<<std::endl;
    }
      
    else{
      for(int i = 0; i < n - 1; i++){
        if(crossWithPolygon(ob.at(i).first, ob.at(i).second, ob.at(i + 1).first, ob.at(i + 1).second, current_path_polygon)){
          current_effective_obs_index.push_back(ob_index);
          // std::cout<<"    line cross with the path polygon, ("<<ob.at(i).first<<", "<<ob.at(i).second<<")  ("<<ob.at(i + 1).first<<", "<<ob.at(i + 1).second<<")"<<std::endl;
          break;
        }
      }
    }
  }

  double current_origin_x = key_x.front();
  double current_origin_y = key_y.front();
  double current_tangent_theta = key_theta.front() - M_PI / 2;
  double last_origin_x = key_x.at(1);
  double last_origin_y = key_y.at(1);
  double last_tangent_theta = key_theta.at(1) - M_PI / 2;

  double min_cost = DBL_MAX;
  for(int j = 0; j < m; j++){//toFix : improve the efficiency
    double d = -sampling_solve_offset + j * sampling_solve_interval;
    double last_x = last_origin_x + d * cos(last_tangent_theta);
    double last_y = last_origin_y + d * sin(last_tangent_theta);

    double current_ob_distance = 100;
    for(int ob_index : current_effective_obs_index){
      double tmp_ob_distance = lineToObstacleDistance(current_origin_x, current_origin_y, last_x, last_y, this->obstacles.at(ob_index));

      if(tmp_ob_distance < current_ob_distance)
        current_ob_distance = tmp_ob_distance;
    }
    //get obs cost
    double current_obs_cost = obs_weight * exp((distance_threshold - current_ob_distance) * exponential_scale);

    //get length cost
    double current_length_cost = length_weight * sqrt((last_x - current_origin_x) * (last_x - current_origin_x) + (last_y - current_origin_y) * (last_y - current_origin_y));
    // std::cout<<"current length cost "<<current_length_cost<<std::endl;
    //get smooth cost
    double section_theta = atan2(last_y - current_origin_y, last_x - current_origin_x);
    double current_smooth_cost = smooth_weight * (last_heading.at(j) - section_theta) * (last_heading.at(j) - section_theta);
    double current_cost = current_length_cost + current_obs_cost + current_smooth_cost;
    
    current_cost += boundary_orientation_weight * (section_theta - key_theta.front()) * (section_theta - key_theta.front());
    if(current_cost + dp.at(1).at(j) < min_cost){
      link.at(0)[0] = d;
      min_cost = current_cost + dp.at(1).at(j);
    }
  }

  
  // std::cout<<"Final cost "<<min_cost<<std::endl;

  auto t2 = std::chrono::steady_clock::now();
  double dr_ms=std::chrono::duration<double,std::milli>(t2-t1).count();
  std::cout<<"Time consumption of sampling optimization: "<<dr_ms<<"ms."<<std::endl;

  //backtrack to fiil key_x and key_y
  double last_d = 0;
  bool success = true;
  for(int i = 1; i < n - 1; i++){
    if(link.at(i-1).count(last_d)){
      double d = link.at(i - 1).at(last_d);//i-1 is the index of previous point
      key_x.at(i) = key_x.at(i) + d * cos(key_theta.at(i) - M_PI / 2);
      key_y.at(i) = key_y.at(i) + d * sin(key_theta.at(i) - M_PI / 2);
      last_d = d;      
    }
    else{
      success = false;
    }
  }

  return success;
}


RefLineOPT::RefLineOPT(Pose start, Pose end, const std::vector<Obstacle>& obs){
  //Connect the start and the terminal to generate a quintic spline curve
  this->initial_path = quinticInterploration(start, end, initial_sampling_margin);

  //first filter out the unnecessary obstacles
  obstacleFilter(obs);

  //extract key points
  extractKeyPoints();

  samplingSolve();

  sampling_solve_offset = 2.0;
  sampling_solve_interval = 0.4;

  //solve
  if(samplingSolve())
    std::cout<<"@@@@ Solve SUCCESSFULLY!  ~(￣▽￣)~* (￣y▽￣)╭ Ohohoho..... @@@@"<<std::endl;
  else
    std::cout<<"@@@@ mission FAILED!  (T__T) (T__T) (T__T) (T__T) (T__T)  @@@@@"<<std::endl;
}


std::shared_ptr<RefLine> RefLineOPT::getPoints(){
  int n = key_x.size();
  CubicSpline csp(this->key_x, this->key_y);
  
  std::shared_ptr<RefLine> ref_line = quinticInterploration(Pose(key_x.front(), key_y.front(), key_theta.front()), Pose(key_x.at(1), key_y.at(1), csp.getKinkHeading(1)), result_sampling_margin);
  std::shared_ptr<RefLine> middle_segments = csp.sampling(1, n - 2, result_sampling_margin);
  ref_line->insert(ref_line->end(), middle_segments->begin(), middle_segments->end());
  std::shared_ptr<RefLine> end_segment = quinticInterploration(Pose(key_x.at(n - 2), key_y.at(n - 2), csp.getKinkHeading(n - 2)), Pose(key_x.back(), key_y.back(), key_theta.back()), result_sampling_margin);
  ref_line->insert(ref_line->end(), end_segment->begin(), end_segment->end());

  return ref_line;
}