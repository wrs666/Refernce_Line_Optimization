#include "ref_line_opt.h"
#include <cassert>
#include <chrono>
#include <array>
#include <unordered_map>
#include <unordered_set>

const double lateral_threshold = 2;
const double longitude_threshold = 4;//TODO: extract params

const double smooth_weight = 1e+3; //cost value about 10
const double boundary_orientation_weight = 5e+3;
const double obs_weight = 1;
const double length_weight = 3e+2;
const double deviation_weight = 10;
const double sigmoid_scale = 10;

const double threshold_for_extraction = 0.5;

const double sampling_solve_offset = 40.0; 
const double sampling_solve_interval = 0.5;
const double sampling_distance_threshold = 1.8;
const double distance_threshold_for_crossing = 1.5;

static const std::string EnumStrings[] = {
  "not_defined",
  "success",
  "maxiter_exceeded",
  "stop_at_tiny_step",
  "stop_at_acceptable_point",
  "local_infeasibility",
  "user_requested_stop",
  "feasible_point_found",
  "diverging_iterates",
  "restoration_failure",
  "error_in_step_computation",
  "invalid_number_detected",
  "too_few_degrees_of_freedom",
  "internal_error",
  "unknown"
};

std::vector<int> combineIndex(std::vector<int> index1, std::vector<int> index2){
  std::vector<int> combined_index;
  combined_index.reserve(index1.size() + index2.size());
  int p1 = 0, p2 = 0;
  while(p1 < index1.size() && p2 < index2.size()){
    if(index1.at(p1) > index2.at(p2)){
      combined_index.push_back(index2.at(p2++));
    }
    else if(index1.at(p1) < index2.at(p2)){
      combined_index.push_back(index1.at(p1++));
    }
    else{
      combined_index.push_back(index1.at(p1++));
      p2++;
    }    
  }
  if(p2 < (int)index2.size()){
    combined_index.insert(combined_index.end(), index2.begin() + p2, index2.end());
  }
  else if(p1 < (int)index1.size()){
    combined_index.insert(combined_index.end(), index1.begin() + p1, index1.end());
  }
  return combined_index;
}

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
  for(int i = 1; i < n - 1; i++){
    B(4 * i - 2, 0) = x.at(i);
    B(4 * i - 1, 0) = y.at(i);
    B(4 * i, 0) = x.at(i);
    B(4 * i + 1, 0) = y.at(i);
    s.at(i) = s.at(i - 1) + sqrt(std::pow(x.at(i) - x.at(i - 1), 2) + std::pow(y.at(i) - y.at(i - 1), 2));


  }
  B(4 * n - 6, 0) = x.back();
  B(4 * n - 5, 0) = y.back();
  s.back() = s.at(n - 2) + sqrt(std::pow(x.back() - x.at(n - 2), 2) + std::pow(y.back() - y.at(n - 2), 2));



  //fill A
  A(0, 0) = 1;
  A(1, 4) = 1;
  double x_start;
  double y_start;
  for(int i = 1; i < n - 1; i++){
    x_start = 4 * i - 2;
    y_start = 8 * (i - 1);
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

  // for(int i = 0; i < N; i++){
  //   for(int j = 0; j < N; j++){
  //     std::cout<<A(i, j)<<" ";
  //   }
  //   std::cout<<std::endl;
  // }

  X = A.colPivHouseholderQr().solve(B);

  for(int i = 0; i < n - 1; i++){
    for(int j = 0; j < 4; j++){
      this->a.at(i).at(j) = X(8 * i + j, 0);
      this->b.at(i).at(j) = X(8 * i + 4 + j, 0);
    }
  }

  // for(int k = 0; k < N; k++){
  //   std::cout<<X(k, 0)<<" ";
  // }
  // std::cout<<std::endl;

  //sampling
  // std::vector<TrackPoint> sampling_points = {};
  // double station = 0;
  // int segment_index = 0;
  // while(station < s.back()){
  //   while(station > s.at(segment_index + 1)){      
  //     segment_index++;
  //     if(segment_index == n - 1)
  //         break;
  //   }
  //   if(segment_index < n - 1){
  //     sampling_points.emplace_back(X(8 * segment_index, 0) * 1 + X(8 * segment_index + 1, 0) * station + X(8 * segment_index + 2, 0) * station * station + X(8 * segment_index + 3, 0) * station * station * station, 
  // X(8 * segment_index + 4, 0) * 1 + X(8 * segment_index + 5, 0) * station + X(8 * segment_index + 6, 0) * station * station + X(8 * segment_index + 7, 0) * station * station * station, 
  // atan2(X(8 * segment_index + 5, 0) + X(8 * segment_index + 6, 0) * station * 2 + X(8 * segment_index + 7, 0) * station * station * 3, X(8 * segment_index + 1, 0) + X(8 * segment_index + 2, 0) * station * 2 + X(8 * segment_index + 3, 0) * station * station * 3), 0);
  //     station += point_margin;
  //     std::cout<<"add one point: x "<<sampling_points.back().x_<<", y "<<sampling_points.back().y_<<std::endl;
  //   }
  // }
  // std::cout<<"Finish interploration and sampling, points size "<<(int)sampling_points.size()<<std::endl;
  // return sampling_points;
}

std::vector<TrackPoint> CubicSpline::sampling(double point_margin){
  //sampling
  int n = this->s.size();
  std::vector<TrackPoint> sampling_points = {};
  double station = this->s.at(1);
  int segment_index = 1;
  while(station < s.at(n - 2)){
    while(station > s.at(segment_index + 1)){      
      segment_index++;
      if(segment_index >= n - 2)
          break;
    }
    if(segment_index < n - 1){
      sampling_points.emplace_back(a.at(segment_index).at(0) * 1 + a.at(segment_index).at(1) * station + a.at(segment_index).at(2) * station * station + a.at(segment_index).at(3) * station * station * station, 
                                   b.at(segment_index).at(0) * 1 + b.at(segment_index).at(1) * station + b.at(segment_index).at(2) * station * station + b.at(segment_index).at(3) * station * station * station, 
                                   atan2(b.at(segment_index).at(1) + b.at(segment_index).at(2) * station * 2 + b.at(segment_index).at(3) * station * station * 3, a.at(segment_index).at(1) + a.at(segment_index).at(2) * station * 2 + a.at(segment_index).at(3) * station * station * 3), 0);
      station += point_margin;
      // std::cout<<"add one point: x "<<sampling_points.back().x_<<", y "<<sampling_points.back().y_<<std::endl;
    }
  }
  std::cout<<"Finish interploration and sampling, points size "<<(int)sampling_points.size()<<std::endl;
  return sampling_points;
}

std::vector<double> CubicSpline::getKinkHeading(){
  std::vector<double> heading(this->s.size());
  for(int i = 0; i < this->s.size() - 1; i++){
    heading.at(i) = atan2(b.at(i).at(1) + b.at(i).at(2) * s.at(i) * 2 + b.at(i).at(3) * s.at(i) * s.at(i) * 3, a.at(i).at(1) + a.at(i).at(2) * s.at(i) * 2 + a.at(i).at(3) * s.at(i) * s.at(i) * 3);
  }
  heading.back() = atan2(b.back().at(1) + b.back().at(2) * s.back() * 2 + b.back().at(3) * s.back() * s.back() * 3, a.back().at(1) + a.back().at(2) * s.back() * 2 + a.back().at(3) * s.back() * s.back() * 3);
  return heading;
}

double CubicSpline::getKinkCurvature(int i){
  double x1 = a.at(i).at(1) + a.at(i).at(2) * s.at(i) * 2 + a.at(i).at(3) * s.at(i) * s.at(i) * 3;
  double x2 = a.at(i).at(2) * 2 + a.at(i).at(3) * s.at(i) * 6;
  double y1 = b.at(i).at(1) + b.at(i).at(2) * s.at(i) * 2 + b.at(i).at(3) * s.at(i) * s.at(i) * 3;
  double y2 = b.at(i).at(2) * 2 + b.at(i).at(3) * s.at(i) * 6;
  return ((x1 * y2 - x2 * y1) / std::pow((x1 * x1 + y1 * y1), 1.5));
}


FG_eval::FG_eval(RefLineOPT *opt){
  this->opt_ = opt;
  int n = opt->points_x.size();
  this->obs_key_points.reserve(n);
  this->key_y.reserve(n);
  this->key_x.reserve(n);
  this->key_theta.reserve(n);
  this->key_y.resize(0);
  this->key_x.resize(0);
  this->key_theta.resize(0);
  this->obs_key_points.resize(0);
  key_x.push_back(opt->points_x.front());
  key_y.push_back(opt->points_y.front());
  key_theta.push_back(opt->init_heading.front());//key points is 2 more than obs key points

  int obs_size = opt->obs_.size();

  bool selected = false, last_selected = false;
  bool refresh = true;
  bool increasing = false;
  std::vector<KeyPoint> kps, last_kps;
  double distance, last_distance = 1000;
  std::vector<int> last_active_obs_index;

  for(int i = 1; i < n - 1; i++){//for each way points
    double relative_x;
    double relative_y;
    KeyPoint kp;
    double MIN_X = 100;//TODO: extract params 
    double MAX_X = -100;
    double MIN_Y = 100;
    double MAX_Y = -100;
    distance = 1000;
    selected = false;
    refresh = false;
    kps.clear();
    std::vector<int> active_obs_index = {};
    // double theta;
    // std::cout<<"For this point : "<<std::endl;
    for(int j = 0; j < obs_size; j++){//for each obstacle in this freesapce
      MIN_X = 100;//TODO: extract params 
      MAX_X = -100;
      MIN_Y = 100;
      MAX_Y = -100;
      // std::cout<<"  For this obstacle : "<<std::endl;
      for(const std::pair<double, double>& p : opt->obs_.at(j)){
        relative_x = cos(opt->init_heading.at(i)) * (p.first - opt->points_x.at(i)) + sin(opt->init_heading.at(i)) * (p.second - opt->points_y.at(i));
        relative_y = -sin(opt->init_heading.at(i)) * (p.first - opt->points_x.at(i)) + cos(opt->init_heading.at(i)) * (p.second - opt->points_y.at(i));
        if(relative_x > MAX_X){
            kp.fore_most = p;
            MAX_X = relative_x;
        }
        if(relative_x < MIN_X){
            kp.rear_most = p;
            MIN_X = relative_x;
        }
        if(relative_y > MAX_Y){
            kp.left_most = p;
            MAX_Y = relative_y;
        }
        if(relative_y < MIN_Y){
            kp.right_most = p;
            MIN_Y = relative_y;
        }
      }

      distance = std::min(distance, std::pow((MIN_X + MAX_X) / (MAX_X - MIN_X + 3) , 2) + std::pow((MIN_Y + MAX_Y) / (MAX_Y - MIN_Y + 10), 2));

      if(std::pow((MIN_X + MAX_X) / (MAX_X - MIN_X + 4) , 2) < 1 && std::pow((MIN_Y + MAX_Y) / (MAX_Y - MIN_Y + 2), 2) < 1){
          selected = true;
          refresh = !last_selected;
      }

      if(std::pow((MIN_X + MAX_X) / (MAX_X - MIN_X + 6) , 2) < 1 && std::pow((MIN_Y + MAX_Y) / (MAX_Y - MIN_Y + 6), 2) < 1){
        kps.push_back(kp);
        // active_obs_index.push_back(j);
      }
      
      if(std::pow((MIN_X + MAX_X) / (MAX_X - MIN_X + 80) , 2) < 1 && std::pow((MIN_Y + MAX_Y) / (MAX_Y - MIN_Y + 80), 2) < 1){
        // kps.push_back(kp);
        active_obs_index.push_back(j);
      }
        
    }

    if(selected){
      std::cout<<i<<"th point, distance "<<distance<<std::endl;
      if(refresh){//start
        key_x.push_back(opt->points_x.at(i));
        key_y.push_back(opt->points_y.at(i));
        key_theta.push_back(opt->init_heading.at(i));
        this->obs_key_points.push_back(kps);
        this->effective_obs_index.push_back(active_obs_index);
        std::cout<<i<<"th point added to key points "<<std::endl;
      }
      else{
        if(increasing && distance < last_distance){//insert maximum
          key_x.push_back(opt->points_x.at(i - 1));
          key_y.push_back(opt->points_y.at(i - 1));
          key_theta.push_back(opt->init_heading.at(i - 1));
          this->obs_key_points.push_back(last_kps);
          this->effective_obs_index.push_back(last_active_obs_index);
          std::cout<<i - 1<<"th point added to key points "<<std::endl;
        }
        else if(!increasing && distance > last_distance){//insert minimum
          key_x.push_back(opt->points_x.at(i - 1));
          key_y.push_back(opt->points_y.at(i - 1));
          key_theta.push_back(opt->init_heading.at(i - 1));
          this->obs_key_points.push_back(last_kps);
          this->effective_obs_index.push_back(last_active_obs_index);
          std::cout<<i - 1<<"th point added to key points "<<std::endl;
        }
      }
    }
    else if(last_selected){//insert end in one section
      key_x.push_back(opt->points_x.at(i - 1));
      key_y.push_back(opt->points_y.at(i - 1));
      key_theta.push_back(opt->init_heading.at(i - 1));
      this->obs_key_points.push_back(last_kps);
      this->effective_obs_index.push_back(last_active_obs_index);
      std::cout<<i - 1<<"th point added to key points "<<std::endl;
    }

    increasing = (distance > last_distance);

    last_distance = distance;
    last_kps = kps;
    last_selected = selected;
    last_active_obs_index = active_obs_index;
  }
  key_x.push_back(opt->points_x.back());
  key_y.push_back(opt->points_y.back());
  key_theta.push_back(opt->init_heading.back());
  this->start_heading = opt->init_heading.front();
  this->end_heading = opt->init_heading.back();
  this->initial_key_x = key_x;
  this->initial_key_y = key_y;
  std::cout<<"Number of key points selected is "<<(int)this->key_x.size()<<". Number of obs vector recorded is "<<(int)this->obs_key_points.size()<<std::endl;
}


void FG_eval::operator()(ADvector& fg, const ADvector& x)
{
  // assert(fg.size() == 5);
  assert(x.size() == this->key_x.size() * 2);
    // variables
  int N = x.size() / 2;
  std::vector<AD<double>> theta(N);
  std::vector<AD<double>> adjacent_theta(N - 1);
  std::vector<std::array<AD<double>, 2>> adjacent_vec(N - 1);
  std::vector<AD<double>> distance(N - 1);
  for(int index = 1; index < N - 1; index++){
    theta.at(index) = CppAD::atan2(x[N + index + 1] - x[N + index - 1], x[index + 1] - x[index - 1]);
    // std::cout<<"theta "<<index<<" "<<theta.at(index)<<std::endl;
    adjacent_theta.at(index) = CppAD::atan2(x[N + index + 1] - x[N + index], x[index + 1] - x[index]);
    // std::cout<<"adj theta "<<index<<" "<<adjacent_theta.at(index)<<std::endl;
    distance.at(index) = CppAD::sqrt((x[index + 1] - x[index]) * (x[index + 1] - x[index]) + (x[N + index + 1] - x[N + index]) * (x[N + index + 1] - x[N + index]));
    adjacent_vec.at(index).at(0) = x[index + 1] - x[index];
    adjacent_vec.at(index).at(1) = x[index + 1 + N] - x[index + N];
  }
  adjacent_theta.front() = CppAD::atan2(x[N + 1] - x[N], x[1] - x[0]);
  distance.front() = CppAD::sqrt((x[1] - x[0]) * (x[1] - x[0]) + (x[N + 1] - x[N]) * (x[N + 1] - x[N]));
  adjacent_vec.front().at(0) = x[1] - x[0];
  adjacent_vec.front().at(1) = x[N + 1] - x[N];

        
  AD<double> smooth_cost = 0, obs_cost = 0, deviation_cost = 0;
  for(int index = 1; index < N - 1; index++){
    deviation_cost += CppAD::log10( 200 * ( (x[index] - this->initial_key_x.at(index)) * (x[index] - this->initial_key_x.at(index)) + (x[index + N] - this->initial_key_y.at(index)) * (x[index + N] - this->initial_key_y.at(index)) + 1 ));
    // AD<double> kappa_smooth_cost = (adjacent_theta.at(index) - adjacent_theta.at(index - 1)) * (adjacent_theta.at(index) - adjacent_theta.at(index - 1)) / ( (distance.at(index - 1) + distance.at(index)) * (distance.at(index - 1) + distance.at(index)));
    // std::cout<<"kappa smooth cost "<<kappa_smooth_cost<<std::endl;
    // smooth_cost += kappa_smooth_cost;

    // AD<double> current_smooth_cost = (adjacent_vec.at(index).at(0) - adjacent_vec.at(index-1).at(0)) * (adjacent_vec.at(index).at(0) - adjacent_vec.at(index-1).at(0)) + (adjacent_vec.at(index).at(1) - adjacent_vec.at(index-1).at(1)) * (adjacent_vec.at(index).at(1) - adjacent_vec.at(index-1).at(1));
    AD<double> current_smooth_cost = (adjacent_theta.at(index) - adjacent_theta.at(index - 1)) * (adjacent_theta.at(index) - adjacent_theta.at(index - 1));
    std::cout<<"current smooth cost "<<current_smooth_cost<<std::endl;
    smooth_cost += current_smooth_cost;         
    std::cout<<"For key point "<<index<<", theta "<<adjacent_theta.at(index)<<std::endl;
    for(const KeyPoint& obs_p : this->obs_key_points.at(index - 1)){
      // std::cout<<"  Obs information: "<<std::endl;
      // std::cout<<"    left "<<obs_p.left_most.first<<", "<<obs_p.left_most.second;
      // std::cout<<"  right "<<obs_p.right_most.first<<", "<<obs_p.right_most.second;
      // std::cout<<"  front "<<obs_p.fore_most.first<<", "<<obs_p.fore_most.second;
      // std::cout<<"  back "<<obs_p.rear_most.first<<", "<<obs_p.rear_most.second<<std::endl;
      // AD<double> left_y = -CppAD::sin(adjacent_theta.at(index)) * (obs_p.left_most.first - x[index]) + CppAD::cos(adjacent_theta.at(index)) * (obs_p.left_most.second - x[index + N]);
      // AD<double> right_y = -CppAD::sin(adjacent_theta.at(index)) * (obs_p.right_most.first - x[index]) + CppAD::cos(adjacent_theta.at(index)) * (obs_p.right_most.second - x[index + N]);
      // AD<double> fore_x = CppAD::cos(adjacent_theta.at(index)) * (obs_p.fore_most.first - x[index]) + CppAD::sin(adjacent_theta.at(index)) * (obs_p.fore_most.second - x[index + N]);
      // AD<double> rear_x = CppAD::cos(adjacent_theta.at(index)) * (obs_p.rear_most.first - x[index]) + CppAD::sin(adjacent_theta.at(index)) * (obs_p.rear_most.second - x[index + N]);
      
      AD<double> left_y = -CppAD::sin(theta.at(index)) * (obs_p.left_most.first - x[index]) + CppAD::cos(theta.at(index)) * (obs_p.left_most.second - x[index + N]);
      AD<double> right_y = -CppAD::sin(theta.at(index)) * (obs_p.right_most.first - x[index]) + CppAD::cos(theta.at(index)) * (obs_p.right_most.second - x[index + N]);
      AD<double> fore_x = CppAD::cos(theta.at(index)) * (obs_p.fore_most.first - x[index]) + CppAD::sin(theta.at(index)) * (obs_p.fore_most.second - x[index + N]);
      AD<double> rear_x = CppAD::cos(theta.at(index)) * (obs_p.rear_most.first - x[index]) + CppAD::sin(theta.at(index)) * (obs_p.rear_most.second - x[index + N]);
      
      std::cout<<"  left_y "<<left_y<<", right_y "<<right_y<<", fore_x "<<fore_x<<", rear_x "<<rear_x<<std::endl;
      AD<double> ellipse_distance = (left_y + right_y) * (left_y + right_y) / ((left_y - right_y + 2 * lateral_threshold) * (left_y - right_y + 2 * lateral_threshold)) + (fore_x + rear_x) * (fore_x + rear_x) / ((fore_x - rear_x + 2 * longitude_threshold) * (fore_x - rear_x + 2 * longitude_threshold)) - 1.5;
      std::cout<<"  current ellipse_distance "<<ellipse_distance<<std::endl;
        // obs_cost += 1 / (1 + CppAD::exp(ellipse_distance * sigmoid_scale));
      obs_cost += CppAD::exp(-ellipse_distance * sigmoid_scale);
        // obs_cost = ellipse_distance;
    }
  }
  AD<double> current_smooth_cost = (adjacent_theta.front() - this->start_heading) * (adjacent_theta.front() - this->start_heading);
  smooth_cost += current_smooth_cost * 10;
  std::cout<<"current smooth cost "<<current_smooth_cost<<std::endl;
  current_smooth_cost = (adjacent_theta.back() - this->end_heading) * (adjacent_theta.back() - this->end_heading);
  smooth_cost += current_smooth_cost * 10;
  // std::cout<<"last theta "<<adjacent_theta.back()<<", last distance "<<distance.back()<<std::endl;
  std::cout<<"current smooth cost "<<current_smooth_cost<<std::endl;

  std::cout<<"smooth cost "<<smooth_cost<<", obs cost "<<obs_cost<<std::endl;


    fg[0] = smooth_weight * smooth_cost + obs_weight * obs_cost;
    // fg[0] = smooth_cost;
    // fg[0] = obs_cost;
    // fg[1] = x[0];
    // fg[2] = x[N];
    // fg[3] = x[N-1];
    // fg[4] = x[2*N-1];
    return;
}

inline double distToObstacle(double px, double py, Obstacle obs){
  int n = obs.size();
  double d = 100;
  // std::cout<<"^-^ in functoin distToObstacle(), sampling point x "<<px<<" y "<<py<<std::endl;

  for(int i = 0, j = 1; i < n - 1; i = j++){
    std::cout<<"           obs 1 x "<<obs.at(i).first<<", y "<<obs.at(i).second<<" ; 2 x "<<obs.at(i + 1).first<<", y "<<obs.at(i + 1).second<<std::endl;
    std::pair<double, double> ip(px - obs.at(i).first, py - obs.at(i).second);
    std::pair<double, double> ij(obs.at(j).first - obs.at(i).first, obs.at(j).second - obs.at(i).second);
    std::pair<double, double> pj(obs.at(j).first - px, obs.at(j).second - py);
    double r = (ip.first * ij.first + ip.second * ij.second) / (ij.first * ij.first + ij.second * ij.second);//point to line
    std::cout<<"           r "<<r<<std::endl;
    double current_d;
    if(r < 0)
      current_d = sqrt(ip.first * ip.first + ip.second * ip.second);
    else if(r > 1){
      current_d = sqrt(pj.first * pj.first + pj.second * pj.second);
    }
    else{
      current_d = fabs((obs.at(j).second - obs.at(i).second) * px + (obs.at(i).first - obs.at(j).first) * py + obs.at(j).first * obs.at(i).second - obs.at(i).first * obs.at(j).second) / 
          sqrt((obs.at(j).second - obs.at(i).second) * (obs.at(j).second - obs.at(i).second) + (obs.at(i).first - obs.at(j).first) * (obs.at(i).first - obs.at(j).first));
    }
    std::cout<<"           current distance "<<current_d<<std::endl;
    d = std::min(d, current_d);
  }
      
  // std::cout<<"    distance "<<d<<std::endl;
  return d;
}

// inline bool isRepelling(double max_x, double min_x, double max_y, double min_y, std::pair<double, double> p1, std::pair<double, double> p2){
//   bool on_left = p1.first < min_x && p2.first < min_x;
//   bool on_right = p1.first > max_x && p2.first > max_x;
//   bool above = p1.second > max_y && p2.second > max_y;
//   bool below = p1.second < min_y && p2.second < min_y;
//   return on_left || on_right || above || below;
// }

inline bool isCross(double p1_x, double p1_y, double p2_x, double p2_y, Obstacle obs){
  int n = obs.size();
  bool cross = false;

  double max_x = std::max(p1_x, p2_x) + 1.5;
  double min_x = std::min(p1_x, p2_x) - 1.5;
  double max_y = std::max(p1_y, p2_y) + 1.5;
  double min_y = std::min(p1_y, p2_y) - 1.5;

  std::pair<double, double> v(p2_x - p1_x, p2_y - p1_y);

  for(int i = 0; i < n - 1; i++){

    bool on_left = obs.at(i).first < min_x && obs.at(i + 1).first < min_x;
    bool on_right = obs.at(i).first > max_x && obs.at(i + 1).first > max_x;
    bool above = obs.at(i).second > max_y && obs.at(i + 1).second > max_y;
    bool below = obs.at(i).second < min_y && obs.at(i + 1).second < min_y;
    if(on_left || on_right || above || below)
      continue;

    std::pair<double, double> v1i(obs.at(i).first - p1_x, obs.at(i).second - p1_y);
    std::pair<double, double> v1j(obs.at(i + 1).first - p1_x, obs.at(i + 1).second - p1_y);
    std::pair<double, double> vij(obs.at(i + 1).first - obs.at(i).first, obs.at(i + 1).second - obs.at(i).second);
    std::pair<double, double> vi1(p1_x - obs.at(i).first, p1_y - obs.at(i).second);
    std::pair<double, double> vi2(p2_x - obs.at(i).first, p2_y - obs.at(i).second);
    double m1iv = v1i.first * v.second - v.first * v1i.second;
    double m1jv = v1j.first * v.second - v.first * v1j.second;
    double mi1ij = vi1.first * vij.second - vij.first * vi1.second;
    double mi2ij = vi2.first * vij.second - vij.first * vi2.second;
    if(m1iv * m1jv <= 0 && mi1ij * mi2ij <= 0){
      cross = true;
      break;
    }
    else{
      double p2line = fabs((p2_y - p1_y) * obs.at(i).first + (p1_x - p2_x) * obs.at(i).second + p2_x * p1_y - p1_x * p2_y) / sqrt((p2_y - p1_y) * (p2_y - p1_y) + (p1_x - p2_x) * (p1_x - p2_x));
      // if(fabs((p2_y - p1_y) * obs.at(i).first + (p1_x - p2_x) * obs.at(i).second + p2_x * p1_y - p1_x * p2_y) / sqrt((p2_y - p1_y) * (p2_y - p1_y) + (p1_x - p2_x) * (p1_x - p2_x)) < distance_threshold_for_crossing){
      // std::cout<<"    obs point x "<<obs.at(i).first<<", y "<<obs.at(i).second<<".  distance from point to line "<<p2line<<std::endl;
      if(p2line < distance_threshold_for_crossing){
        cross = true;
        break;
      }
    }
  }
  if(!cross){
    // cross = fabs((p2_y - p1_y) * obs.back().first + (p1_x - p2_x) * obs.back().second + p2_x * p1_y - p1_x * p2_y) / sqrt((p2_y - p1_y) * (p2_y - p1_y) + (p1_x - p2_x) * (p1_x - p2_x)) < distance_threshold_for_crossing;
    double p2line = fabs((p2_y - p1_y) * obs.back().first + (p1_x - p2_x) * obs.back().second + p2_x * p1_y - p1_x * p2_y) / sqrt((p2_y - p1_y) * (p2_y - p1_y) + (p1_x - p2_x) * (p1_x - p2_x));
      // if(fabs((p2_y - p1_y) * obs.at(i).first + (p1_x - p2_x) * obs.at(i).second + p2_x * p1_y - p1_x * p2_y) / sqrt((p2_y - p1_y) * (p2_y - p1_y) + (p1_x - p2_x) * (p1_x - p2_x)) < distance_threshold_for_crossing){
    // std::cout<<"    obs back point x "<<obs.back().first<<", y "<<obs.back().second<<".  distance from point to line "<<p2line<<std::endl;
      if(p2line < distance_threshold_for_crossing){
        cross = true;
      }
  }
  return cross;
}

void FG_eval::samplingSolve(double max_offset, double sampling_interval, CubicSpline& csp){
  std::cout<<"************Implementing Sampling Solving***********"<<std::endl;
  auto t1 = std::chrono::steady_clock::now();

  int m = 2 * (max_offset / sampling_interval) + 1;
  int n = this->key_x.size();
  // std::cout<<"Number of sampling states per layer "<<m<<". Key points size "<<n<<std::endl;
  std::vector<std::vector<double>> dp(n, std::vector<double>(m, 0));
  std::vector<std::unordered_map<double, double>> link(n);
  std::vector<double> last_heading(m), current_heading(m);
  std::vector<std::unordered_set<int>> black_lists(n);
  std::vector<std::vector<double>> cumulative_length(n, std::vector<double>(m, DBL_MAX));

  //state trandsistion from n-2 to 1, leave alone n-1 and 0
  //index : n - 2
  // double d;//lateral movement distance of key point in current level
  // double last_d;//lateral movement distance of key point in last level
  // double moved_x, moved_y, last_x, last_y;//coordinates of points in this level and last level
  // double section_theta;//heading connecting current level to the last level
  for(int i = 0; i < m; i++){//toFix : improve the efficiency
    double d = -max_offset + i * sampling_interval;
    double moved_x = key_x.at(n - 2) + d * cos(key_theta.at(n - 2) - M_PI / 2);
    double moved_y = key_y.at(n - 2) + d * sin(key_theta.at(n - 2) - M_PI / 2);
    double section_theta = atan2(key_y.back() - moved_y, key_x.back() - moved_x);
    last_heading.at(i) = section_theta;
    double min_ob_distance = 100;
    for(int ob_index : effective_obs_index.back()){
      // std::cout<<"For the end point, offset "<<d<<std::endl;
      if(isCross(moved_x, moved_y, key_x.back(), key_y.back(), this->opt_->obs_.at(ob_index))){
        min_ob_distance = 0;
        black_lists.at(n - 2).insert(i);
        // std::cout<<"CROSS! offset "<<d<<std::endl;
        break;
      }
      else{
        // double dist = distToObstacle(moved_x, moved_y, this->opt_->obs_.at(ob_index));
        // std::cout<<"For the end point, offset "<<d<<", distance to obstacle "<<dist<<std::endl;
        min_ob_distance = std::min(min_ob_distance, distToObstacle(moved_x, moved_y, this->opt_->obs_.at(ob_index)));
      }

    }
    double smooth_cost = boundary_orientation_weight * (section_theta - end_heading) * (section_theta - end_heading);
    // std::cout<<"smoooth cost "<<smooth_cost<<std::endl;
    double obs_cost = obs_weight * exp((sampling_distance_threshold - min_ob_distance) * sigmoid_scale);
    // std::cout<<"obs cost "<<obs_cost<<std::endl;
    double length_cost = length_weight * sqrt((key_x.back() - moved_x) * (key_x.back() - moved_x) + (key_y.back() - moved_y) * (key_y.back() - moved_y));
    // std::cout<<"length cost "<<length_cost<<std::endl;
    dp.at(n - 2).at(i) = smooth_cost + obs_cost + length_cost;
  }

  //index : n-3 to 1
  for(int i = n - 3; i > 0; i--){//toFix : improve the efficiency
    std::cout<<"Sampling "<<i<<"th point"<<std::endl;
    const std::vector<int>& active_obs_index = effective_obs_index.at(i - 1);
    std::vector<int> section_obs_index = combineIndex(effective_obs_index.at(i - 1), effective_obs_index.at(i));
    for(int j = 0; j < m; j++){//for the current level
      double d = -max_offset + j * sampling_interval;
      double moved_x = key_x.at(i) + d * cos(key_theta.at(i) - M_PI / 2);
      double moved_y = key_y.at(i) + d * sin(key_theta.at(i) - M_PI / 2);
      std::cout<<"  offset "<<d<<std::endl;

      //get obs cost
      double min_ob_distance = 100, ob_distacne;
      for(int ob_index : active_obs_index){
        min_ob_distance = std::min(min_ob_distance, distToObstacle(moved_x, moved_y, this->opt_->obs_.at(ob_index)));
        if(min_ob_distance == 0)
          break;
      }
      std::cout<<"    min distance to ob "<<min_ob_distance<<std::endl;

      //for each states in the last level
      double min_cost = DBL_MAX;
      for(int last_j = 0; last_j < m; last_j++){
        if(black_lists.at(i + 1).count(last_j))
          continue;
        double last_d = -max_offset + last_j * sampling_interval;
        double last_x =  key_x.at(i + 1) + last_d * cos(key_theta.at(i + 1) - M_PI / 2);
        double last_y =  key_y.at(i + 1) + last_d * sin(key_theta.at(i + 1) - M_PI / 2);
        // std::cout<<"    last level offset "<<last_d<<std::endl;
        bool cross = false;
        for(int soi : section_obs_index){
          if(isCross(moved_x, moved_y, last_x, last_y, this->opt_->obs_.at(soi))){//Collision occurs between two Point Connection and Obstacle
            cross = true;
            break;
            // std::cout<<"    CROSS! last level offset "<<last_d<<std::endl;
          }
        }
        if(cross)
          continue;

        //get length cost
        double current_length_cost = length_weight * sqrt((last_x - moved_x) * (last_x - moved_x) + (last_y - moved_y) * (last_y - moved_y));
        // std::cout<<"    current length cost "<<current_length_cost<<std::endl;
        //get smooth cost
        double section_theta = atan2(last_y - moved_y, last_x - moved_x);
        double current_cost = current_length_cost + smooth_weight * (last_heading.at(last_j) - section_theta) * (last_heading.at(last_j) - section_theta);
        if(current_cost + dp.at(i + 1).at(last_j) < min_cost){
          current_heading.at(j) = section_theta;
          link.at(i)[d] = last_d;
          min_cost = current_cost + dp.at(i + 1).at(last_j);
        }
      }
      // if(link.at(i).count(d))
        // std::cout<<"    optimal offset of last level "<<link.at(i).at(d)<<std::endl;
      if(link.at(i).count(d) == 0){
        // std::cout<<"    No feasible path "<<std::endl;
        black_lists.at(i).insert(j);
        // std::cout<<i<<" level, "<<j<<" state is put into black list"<<std::endl;
      }
        

      //refresh dp array
      dp.at(i).at(j) = min_cost + obs_weight * exp((sampling_distance_threshold - min_ob_distance) * sigmoid_scale);
    }
    last_heading = current_heading;
  }

  //index 1
  double min_cost = DBL_MAX;
  for(int j = 0; j < m; j++){//toFix : improve the efficiency
    if(black_lists.at(1).count(j))
      continue;
    double d = -max_offset + j * sampling_interval;
    double moved_x = key_x.at(1) + d * cos(key_theta.at(1) - M_PI / 2);
    double moved_y = key_y.at(1) + d * sin(key_theta.at(1) - M_PI / 2);

    bool cross = false;
    for(int soi : effective_obs_index.front()){
      if(isCross(key_x.front(), key_y.front(), moved_x, moved_y, this->opt_->obs_.at(soi))){//Collision occurs between two Point Connection and Obstacle
        cross = true;
        break;
        // std::cout<<"CROSS! last level offset "<<d<<std::endl;
      }

    }
    if(cross)
      continue;

    //get length cost
    double current_length_cost = length_weight * sqrt((moved_x - key_x.front()) * (moved_x - key_x.front()) + (moved_y - key_y.front()) * (moved_y - key_y.front()));
    // std::cout<<"current length cost "<<current_length_cost<<std::endl;
    //get smooth cost
    double section_theta = atan2(moved_y - key_y.front(), moved_x - key_x.front());
    double current_cost = current_length_cost + smooth_weight * (last_heading.at(j) - section_theta) * (last_heading.at(j) - section_theta);
    current_cost += boundary_orientation_weight * (section_theta - start_heading) * (section_theta - start_heading);
    if(current_cost + dp.at(1).at(j) < min_cost){
      link.at(0)[0] = d;
      min_cost = current_cost + dp.at(1).at(j);
    }
  }

  
  std::cout<<"Final cost "<<min_cost<<std::endl;

  auto t2 = std::chrono::steady_clock::now();
  double dr_ms=std::chrono::duration<double,std::milli>(t2-t1).count();
  std::cout<<"Time consumption of sampling optimization: "<<dr_ms<<"ms."<<std::endl;

  //backtrack to fiil key_x and key_y
  double last_d = 0;
  for(int i = 1; i < n - 1; i++){
    if(link.at(i-1).count(last_d)){
      double d = link.at(i - 1).at(last_d);//i-1 is the index of previous point
      key_x.at(i) = key_x.at(i) + d * cos(key_theta.at(i) - M_PI / 2);
      key_y.at(i) = key_y.at(i) + d * sin(key_theta.at(i) - M_PI / 2);
      last_d = d;      
    }
    else{
      std::cout<<"Solutioan Failed, level "<<(i - 1)<<std::endl;
    }

  }

  //cubic interploration
  auto t3 = std::chrono::steady_clock::now();
  csp = CubicSpline(key_x, key_y);
  auto t4 = std::chrono::steady_clock::now();
  dr_ms=std::chrono::duration<double,std::milli>(t4-t3).count();
  std::cout<<"Time consumption of cubic interploration: "<<dr_ms<<"ms."<<std::endl;
}

void FG_eval::refreshObsKey(){//when x, y, theta are refreshed, refresh obs key points
  int n = this->key_x.size();
  this->obs_key_points.clear();
  this->obs_key_points.resize(n - 2);
  for(int i = 1; i < n - 1; i++){//for each way points
    double relative_x;
    double relative_y;
    double MIN_X = 50;//TODO: extract params 
    double MAX_X = -50;
    double MIN_Y = 50;
    double MAX_Y = -50;
    int obs_size = this->opt_->obs_.size();
    
    std::cout<<"*********"<<i<< "th point, theta "<<this->key_theta.at(i)<<std::endl;
    for(int j = 0; j < obs_size; j++){//for each obstacle in this freesapce
    std::cout<<"for this obs"<<std::endl;
      MIN_X = 50;//TODO: extract params 
      MAX_X = -50;
      MIN_Y = 50;
      MAX_Y = -50;
      // std::cout<<"  For this obstacle : "<<std::endl;
      KeyPoint kp;
      for(const std::pair<double, double>& p : this->opt_->obs_.at(j)){
        std::cout<<"  obs point in world "<<p.first<<", "<<p.second;
        relative_x = cos(this->key_theta.at(i)) * (p.first - this->key_x.at(i)) + sin(this->key_theta.at(i)) * (p.second - this->key_y.at(i));
        relative_y = -sin(this->key_theta.at(i)) * (p.first - this->key_x.at(i)) + cos(this->key_theta.at(i)) * (p.second - this->key_y.at(i));
        std::cout<<"  in relative frame "<<relative_x<<", "<<relative_y<<std::endl;
        if(relative_x > MAX_X){
            kp.fore_most = p;
            MAX_X = relative_x;
        }
        if(relative_x < MIN_X){
            kp.rear_most = p;
            MIN_X = relative_x;
        }
        if(relative_y > MAX_Y){
            kp.left_most = p;
            MAX_Y = relative_y;
        }
        if(relative_y < MIN_Y){
            kp.right_most = p;
            MIN_Y = relative_y;
        }
      }

      if(std::pow((MIN_X + MAX_X) / (MAX_X - MIN_X + 6) , 2) < 1 && std::pow((MIN_Y + MAX_Y) / (MAX_Y - MIN_Y + 6), 2) < 1){
        this->obs_key_points.at(i - 1).push_back(kp);
      }

    }
  }
}

void RefLineOPT::obstacleFilter(const std::vector<Obstacle>& obs){    
    std::pair<double, double> start(points_x.front(), points_y.front());
    std::pair<double, double> end(points_x.back(), points_y.back());
    double start_to_end = distance(start, end);
    this->obs_.resize(0);
    for(Obstacle ob : obs){
        for(std::pair<double, double> p : ob){
            if(distance(p, start) < start_to_end && distance(p, end) < start_to_end){
                if(ob.back() == ob.front())
                    ob.pop_back();
                this->obs_.push_back(ob);
                break;
            }
        }
    }
}

RefLineOPT::RefLineOPT(const std::vector<Point>& points, const std::vector<Obstacle>& obs){
  int n = points.size();
  points_x.resize(n);
  points_y.resize(n);
  init_heading.resize(n);
  points_x.front() = points.front().x_;
  points_y.front() = points.front().y_;
  points_x.back() = points.back().x_;
  points_y.back() = points.back().y_;
  init_heading.front() = atan2(points.at(1).y_ - points.front().y_, points.at(1).x_ - points.front().x_);
  init_heading.back() = atan2(points.back().y_ - points.at(n-2).y_, points.back().x_ - points.at(n-2).x_);

  for(int i = 1; i < n - 1; i++){
      points_x.at(i) = points.at(i).x_;
      points_y.at(i) = points.at(i).y_;
      init_heading.at(i) = atan2(points.at(i+1).y_ - points.at(i-1).y_, points.at(i+1).x_ - points.at(i-1).x_);
  }

  obstacleFilter(obs);
  fg_eval = FG_eval(this);

}

RefLineOPT::RefLineOPT(const std::vector<TrackPoint>& track_points, const std::vector<Obstacle>& obs){
  int n = track_points.size();
  points_x.resize(n);
  points_y.resize(n);
  init_heading.resize(n);

  for(int i = 0; i < n; i++){
      points_x.at(i) = track_points.at(i).x_;
      points_y.at(i) = track_points.at(i).y_;
      init_heading.at(i) = track_points.at(i).theta_;
  }

  obstacleFilter(obs);
  fg_eval = FG_eval(this); 
}

bool RefLineOPT::ssolve(CubicSpline& csp){
  this->fg_eval.samplingSolve(sampling_solve_offset, sampling_solve_interval, csp);
  return true;
}

bool RefLineOPT::solve(CubicSpline& csp){
  typedef CPPAD_TESTVECTOR(double) Dvector;
  auto t1 = std::chrono::steady_clock::now();

  int max_iter = 5;
  int min_iter = 2;
  int iter_time = 0;

  bool ok = false;
  

  // object that computes objective and constraints
  // options
  std::string options;
  // turn off any printing
  //options += "Retape  true\n";
  options += "Integer print_level  0\n";
  options += "String sb            yes\n";
  // maximum iterations
  options += "Integer max_iter     100\n";
  //approximate accuracy in first order  conditions;
  // see Mathematical Programming, Volume 106, Number 1,
  // Pages 25-57, Equation (6)
  options += "Numeric tol          0.1\n";
  //derivative tesing
  options += "String derivative_test   second-order\n";
  // maximum amount of random pertubation; e.g.,
  // when evaluation finite diff
  options += "Numeric point_perturbation_radius   0.\n";

  double final_F;

  while(iter_time <= min_iter){
    iter_time++;

    int N = this->fg_eval.key_x.size();
    if(N != this->fg_eval.key_y.size() || N != this->fg_eval.key_theta.size())
        throw("Unequal number of coordinates or heading before optimization");

    size_t nx = N * 2; // number of varibles
    size_t ng = 0; // number of constraints
    Dvector x0(nx), xl(nx), xu(nx); // initial condition of varibles
    Dvector gl(ng), gu(ng);// lower and upper bounds for constraints

    x0[0] = this->points_x.front();
    x0[N] = this->points_y.front();
    x0[N - 1] = this->points_x.back();
    x0[2 * N - 1] = this->points_y.back();

    xl[0] = this->points_x.front();
    xu[0] = this->points_x.front();
    xl[N-1] = this->points_x.back();
    xu[N-1] = this->points_x.back();
    xl[N] = this->points_y.front();
    xu[N] = this->points_y.front();
    xl[2 * N - 1] = this->points_y.back();
    xu[2 * N - 1] = this->points_y.back();    

    CppAD::ipopt::solve_result<Dvector> solution;

    for(int i = 1; i < N - 1; i++){//no need to transfer by fg_eval.key_point???
      x0[i] = this -> fg_eval.key_x.at(i);
      x0[i + N] = this -> fg_eval.key_y.at(i);
      xl[i] = this -> fg_eval.key_x.at(i) - 3.0;
      xl[i + N] = this -> fg_eval.key_y.at(i) - 3.0;
      xu[i] = this -> fg_eval.key_x.at(i) + 3.0;
      xu[i + N] = this -> fg_eval.key_y.at(i) + 3.0;
    }


    auto t2=std::chrono::steady_clock::now();
    CppAD::ipopt::solve<Dvector, FG_eval>(options, x0, xl, xu, gl, gu, this->fg_eval, solution); // solve the problem
    auto t3=std::chrono::steady_clock::now();
    double dr_ms=std::chrono::duration<double,std::milli>(t3-t2).count();
    std::cout<<"Time consumption of optimization: "<<dr_ms<<"ms."<<std::endl;

    if(solution.status == CppAD::ipopt::solve_result<Dvector>::status_type::success)
      ok = true;
    else{
      std::cout<<"\033[31;1mThe Solution was NOT Successfull\033[0m"<<std::endl;
      std::cout<<"Status "<<EnumStrings[solution.status]<<std::endl;
    }

    //refresh x and y
    for(int i = 1; i < N - 1; i++){
      this->fg_eval.key_x.at(i) = solution.x[i];
      this->fg_eval.key_y.at(i) = solution.x[i + N];
    }
    
    // std::vector<double> extract_x, extract_y;
    // refreshKeyPoints(threshold_for_extraction, extract_x, extract_y, this->fg_eval.key_x, this->fg_eval.key_y);

    //cubic interploration
    csp = CubicSpline(this->fg_eval.key_x, this->fg_eval.key_y);
    if(iter_time == min_iter){
      break;
    }
    //refresh heading
    // this->fg_eval.key_theta = csp.getKinkHeading();
    correctHeading();
    this->fg_eval.key_theta.front() = this->fg_eval.start_heading;
    this->fg_eval.key_theta.back() = this->fg_eval.end_heading;
    //refresh obs key points
    this->fg_eval.refreshObsKey();

    final_F = solution.obj_value;
  }

  csp = CubicSpline(this->fg_eval.key_x, this->fg_eval.key_y);

  auto t4=std::chrono::steady_clock::now();
  //毫秒级
  double whole_ms=std::chrono::duration<double,std::milli>(t4-t1).count();

  std::cout<<"the whole process takes "<<whole_ms<<"ms."<<std::endl;

  std::cout<<"Final result F "<<final_F<<std::endl;
  return ok;
}

void RefLineOPT::correctHeading(){
  int n = this->fg_eval.key_x.size();
    for(int i = 1; i < n - 1; i++)
      this->fg_eval.key_theta.at(i) = atan2(this->fg_eval.key_y.at(i + 1) - this->fg_eval.key_y.at(i), this->fg_eval.key_x.at(i + 1) - this->fg_eval.key_x.at(i));
}
// void RefLineOPT::correctHeading(){
//     this->init_heading.resize(points_x.size());
//     for(int i = 1; i < points_x.size() - 1; i++)
//         init_heading.at(i) = atan2(points_y.at(i + 1) - points_y.at(i), points_x.at(i) - points_x.at(i));
// }

void extarctKeyPoints(const std::vector<double>& x, const std::vector<double>& y, double min_distance, int start, int end, std::vector<int>& black_list ){
  if(end - start < 2)
    return;
  double max_distance = 0;
  double current_distance;
  int max_offset_index;
  for(int i = start + 1; i < end; i++){
    current_distance = fabs((y.at(start) - y.at(end)) * x.at(i) + (x.at(end) - x.at(start)) * y.at(i) + x.at(start) * y.at(end) - x.at(end) * y.at(start)) / sqrt(std::pow(y.at(start) - y.at(end), 2) + std::pow(x.at(start) - x.at(end), 2));
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
    extarctKeyPoints(x, y, min_distance, start, max_offset_index, black_list);
    extarctKeyPoints(x, y, min_distance, max_offset_index, end, black_list);
  }
}

void RefLineOPT::refreshKeyPoints(double dist_thre, std::vector<double>& key_x, std::vector<double>& key_y, const std::vector<double>& x, const std::vector<double>& y){
    key_x.resize(0);
    key_y.resize(0);
    std::vector<int> black_list(x.size(), false);
    extarctKeyPoints(x, y, dist_thre, 0, x.size() - 1, black_list);
    for(int i = 0; i < x.size(); i++){
        if(!black_list.at(i)){
            key_x.push_back(x.at(i));
            key_y.push_back(y.at(i));
        }
    }
}

// CubicSpline::CubicSpline(const std::vector<double>& x, const std::vector<double>& y){//toImprove : combine loop
//   int n = x.size();
//   int N = 8 * (n - 1); 

//   this->s = std::vector<double>(n, 0);
  
//   Eigen::MatrixXd A(N, N);
//   A.setZero();
//   Eigen::MatrixXd X(N, 1);
//   Eigen::MatrixXd B(N, 1);
//   B.setZero();

//   //fill B and s
//   B(0, 0) = x.front();
//   B(1, 0) = y.front();
//   for(int i = 1; i < n - 1; i++){
//     B(4 * i - 2, 0) = x.at(i);
//     B(4 * i - 1, 0) = y.at(i);
//     B(4 * i, 0) = x.at(i);
//     B(4 * i + 1, 0) = y.at(i);
//     s.at(i) = s.at(i - 1) + sqrt(std::pow(x.at(i) - x.at(i - 1), 2) + std::pow(y.at(i) - y.at(i - 1), 2));


//   }
//   B(4 * n - 6, 0) = x.back();
//   B(4 * n - 5, 0) = y.back();
//   s.back() = s.at(n - 2) + sqrt(std::pow(x.back() - x.at(n - 2), 2) + std::pow(y.back() - y.at(n - 2), 2));



//   //fill A
//   A(0, 0) = 1;
//   A(1, 4) = 1;
//   double x_start;
//   double y_start;
//   for(int i = 1; i < n - 1; i++){
//     x_start = 4 * i - 2;
//     y_start = 8 * (i - 1);
//     for(int j = 0; j < 4; j++){
//       A(x_start + j, y_start + j * 4) = 1;
//       A(x_start + j, y_start + j * 4 + 1) = s.at(i);
//       A(x_start + j, y_start + j * 4 + 2) = s.at(i) * s.at(i);
//       A(x_start + j, y_start + j * 4 + 3) = s.at(i) * s.at(i) * s.at(i);
//     }
    
//     A(4 * n - 8 + 4 * i, y_start + 1) = -1;
//     A(4 * n - 8 + 4 * i, y_start + 2) = -2 * s.at(i);
//     A(4 * n - 8 + 4 * i, y_start + 3) = -3 * s.at(i) * s.at(i);
//     A(4 * n - 8 + 4 * i, y_start + 9) = 1;
//     A(4 * n - 8 + 4 * i, y_start + 10) = 2 * s.at(i);
//     A(4 * n - 8 + 4 * i, y_start + 11) = 3 * s.at(i) * s.at(i);
//     A(4 * n - 7 + 4 * i, y_start + 2) = -2;
//     A(4 * n - 7 + 4 * i, y_start + 3) = -6 * s.at(i);
//     A(4 * n - 7 + 4 * i, y_start + 10) = 2;
//     A(4 * n - 7 + 4 * i, y_start + 11) = 6 * s.at(i);
//     A(4 * n - 6 + 4 * i, y_start + 5) = -1;
//     A(4 * n - 6 + 4 * i, y_start + 6) = -2 * s.at(i);
//     A(4 * n - 6 + 4 * i, y_start + 7) = -3 * s.at(i) * s.at(i);
//     A(4 * n - 6 + 4 * i, y_start + 13) = 1;
//     A(4 * n - 6 + 4 * i, y_start + 14) = 2 * s.at(i);
//     A(4 * n - 6 + 4 * i, y_start + 15) = 3 * s.at(i) * s.at(i);
//     A(4 * n - 5 + 4 * i, y_start + 6) = -2;
//     A(4 * n - 5 + 4 * i, y_start + 7) = -6 * s.at(i);
//     A(4 * n - 5 + 4 * i, y_start + 14) = 2;
//     A(4 * n - 5 + 4 * i, y_start + 15) = 6 * s.at(i);
//   }

//   A(4 * n - 6, 8 * (n - 2)) = 1;
//   A(4 * n - 6, 8 * (n - 2) + 1) = s.back();
//   A(4 * n - 6, 8 * (n - 2) + 2) = s.back() * s.back();
//   A(4 * n - 6, 8 * (n - 2) + 3) = s.back() * s.back() * s.back();
//   A(4 * n - 5, 8 * (n - 2) + 4) = 1;
//   A(4 * n - 5, 8 * (n - 2) + 5) = s.back();
//   A(4 * n - 5, 8 * (n - 2) + 6) = s.back() * s.back();
//   A(4 * n - 5, 8 * (n - 2) + 7) = s.back() * s.back() * s.back();

//   A(8 * n - 12, 2) = 2;
//   A(8 * n - 11, 6) = 2;
//   A(8 * n - 10, 8 * n - 14) = 2;
//   A(8 * n - 10, 8 * n - 13) = 6 * s.back();
//   A(8 * n - 9, 8 * n - 10) = 2;
//   A(8 * n - 9, 8 * n - 9) = 6 * s.back();

//   // for(int i = 0; i < N; i++){
//   //   for(int j = 0; j < N; j++){
//   //     std::cout<<A(i, j)<<" ";
//   //   }
//   //   std::cout<<std::endl;
//   // }

//   X = A.colPivHouseholderQr().solve(B);

//   for(int i = 0; i < n - 1; i++){
//     for(int j = 0; j < 4; j++){
//       this->a.at(i).at(j) = X(8 * i + j, 0);
//       this->b.at(i).at(j) = X(8 * i + 4 + j, 0);
//     }
//   }

//   // for(int k = 0; k < N; k++){
//   //   std::cout<<X(k, 0)<<" ";
//   // }
//   // std::cout<<std::endl;

//   //sampling
//   // std::vector<TrackPoint> sampling_points = {};
//   // double station = 0;
//   // int segment_index = 0;
//   // while(station < s.back()){
//   //   while(station > s.at(segment_index + 1)){      
//   //     segment_index++;
//   //     if(segment_index == n - 1)
//   //         break;
//   //   }
//   //   if(segment_index < n - 1){
//   //     sampling_points.emplace_back(X(8 * segment_index, 0) * 1 + X(8 * segment_index + 1, 0) * station + X(8 * segment_index + 2, 0) * station * station + X(8 * segment_index + 3, 0) * station * station * station, 
//   // X(8 * segment_index + 4, 0) * 1 + X(8 * segment_index + 5, 0) * station + X(8 * segment_index + 6, 0) * station * station + X(8 * segment_index + 7, 0) * station * station * station, 
//   // atan2(X(8 * segment_index + 5, 0) + X(8 * segment_index + 6, 0) * station * 2 + X(8 * segment_index + 7, 0) * station * station * 3, X(8 * segment_index + 1, 0) + X(8 * segment_index + 2, 0) * station * 2 + X(8 * segment_index + 3, 0) * station * station * 3), 0);
//   //     station += point_margin;
//   //     std::cout<<"add one point: x "<<sampling_points.back().x_<<", y "<<sampling_points.back().y_<<std::endl;
//   //   }
//   // }
//   // std::cout<<"Finish interploration and sampling, points size "<<(int)sampling_points.size()<<std::endl;
//   // return sampling_points;
// }

// std::vector<TrackPoint> CubicSpline::sampling(double point_margin){
//   //sampling
//   int n = this->s.size();
//   std::vector<TrackPoint> sampling_points = {};
//   double station = 0;
//   int segment_index = 0;
//   while(station < s.back()){
//     while(station > s.at(segment_index + 1)){      
//       segment_index++;
//       if(segment_index == n - 1)
//           break;
//     }
//     if(segment_index < n - 1){
//       sampling_points.emplace_back(a.at(segment_index).at(0) * 1 + a.at(segment_index).at(1) * station + a.at(segment_index).at(2) * station * station + a.at(segment_index).at(3) * station * station * station, 
//                                    b.at(segment_index).at(0) * 1 + b.at(segment_index).at(1) * station + b.at(segment_index).at(2) * station * station + b.at(segment_index).at(3) * station * station * station, 
//                                    atan2(b.at(segment_index).at(1) + b.at(segment_index).at(2) * station * 2 + b.at(segment_index).at(3) * station * station * 3, a.at(segment_index).at(1) + a.at(segment_index).at(2) * station * 2 + a.at(segment_index).at(3) * station * station * 3), 0);
//       station += point_margin;
//       // std::cout<<"add one point: x "<<sampling_points.back().x_<<", y "<<sampling_points.back().y_<<std::endl;
//     }
//   }
//   std::cout<<"Finish interploration and sampling, points size "<<(int)sampling_points.size()<<std::endl;
//   return sampling_points;
// }

// std::vector<double> CubicSpline::getKinkHeading(){
//   std::vector<double> heading(this->s.size());
//   for(int i = 0; i < this->s.size() - 1; i++){
//     heading.at(i) = atan2(b.at(i).at(1) + b.at(i).at(2) * s.at(0) * 2 + b.at(i).at(3) * s.at(i) * s.at(i) * 3, a.at(i).at(1) + a.at(i).at(2) * s.at(i) * 2 + a.at(i).at(3) * s.at(i) * s.at(i) * 3);
//   }
//   heading.back() = atan2(b.back().at(1) + b.back().at(2) * s.back() * 2 + b.back().at(3) * s.back() * s.back() * 3, a.back().at(1) + a.back().at(2) * s.back() * 2 + a.back().at(3) * s.back() * s.back() * 3);
//   return heading;
// }

// std::vector<TrackPoint> RefLineOPT::keyInterploration(double point_margin) const {
//     // this->refreshKeyPoints(threshold_for_extraction);
//     return cubic_interploration(this->key_points_x, this->key_points_y, point_margin);
// }

// std::vector<TrackPoint> RefLineOPT::rawInterploration(double point_margin) const {
//     return cubic_interploration(this->points_x, this->points_y, point_margin);
// }