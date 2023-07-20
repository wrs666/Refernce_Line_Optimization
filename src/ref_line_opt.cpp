#include "ref_line_opt.h"
#include <cassert>
#include <chrono>
#include <array>

const double lateral_threshold = 2;
const double longitude_threshold = 4;//TODO: extract params

const double smooth_weight = 1e+5; //cost value about 10
const double obs_weight = 10;
const double deviation_weight = 10;
const double sigmoid_scale = 20;

double threshold_for_extraction = 0.5;

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
  double station = 0;
  int segment_index = 0;
  while(station < s.back()){
    while(station > s.at(segment_index + 1)){      
      segment_index++;
      if(segment_index == n - 1)
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
  for(int i = 1; i < n - 1; i++){//for each way points
    double relative_x;
    double relative_y;
    double MIN_X = 100;//TODO: extract params 
    double MAX_X = -100;
    double MIN_Y = 100;
    double MAX_Y = -100;
    bool selected = false;
    std::vector<KeyPoint> kps;
    KeyPoint kp;
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

      if(std::pow((MIN_X + MAX_X) / (MAX_X - MIN_X + 4) , 2) < 1 && std::pow((MIN_Y + MAX_Y) / (MAX_Y - MIN_Y + 4), 2) < 1){
        selected = true;
      }

       if(std::pow((MIN_X + MAX_X) / (MAX_X - MIN_X + 6) , 2) < 1 && std::pow((MIN_Y + MAX_Y) / (MAX_Y - MIN_Y + 6), 2) < 1){
        kps.push_back(kp);
      }
        
        
    }
    if(selected){
      key_x.push_back(opt->points_x.at(i));
      key_y.push_back(opt->points_y.at(i));
      key_theta.push_back(opt->init_heading.at(i));
      this->obs_key_points.push_back(kps);
    }
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
    AD<double> current_smooth_cost = (adjacent_vec.at(index).at(0) - adjacent_vec.at(index-1).at(0)) * (adjacent_vec.at(index).at(0) - adjacent_vec.at(index-1).at(0)) + (adjacent_vec.at(index).at(1) - adjacent_vec.at(index-1).at(1)) * (adjacent_vec.at(index).at(1) - adjacent_vec.at(index-1).at(1));
    std::cout<<"current smooth cost "<<current_smooth_cost<<std::endl;
    smooth_cost += current_smooth_cost;         
    std::cout<<"For key point "<<index<<", theta "<<adjacent_theta.at(index)<<std::endl;
    for(const KeyPoint& obs_p : this->obs_key_points.at(index - 1)){
      std::cout<<"  Obs information: "<<std::endl;
      std::cout<<"    left "<<obs_p.left_most.first<<", "<<obs_p.left_most.second;
      std::cout<<"  right "<<obs_p.right_most.first<<", "<<obs_p.right_most.second;
      std::cout<<"  front "<<obs_p.fore_most.first<<", "<<obs_p.fore_most.second;
      std::cout<<"  back "<<obs_p.rear_most.first<<", "<<obs_p.rear_most.second<<std::endl;
      AD<double> left_y = -CppAD::sin(adjacent_theta.at(index)) * (obs_p.left_most.first - x[index]) + CppAD::cos(adjacent_theta.at(index)) * (obs_p.left_most.second - x[index + N]);
      AD<double> right_y = -CppAD::sin(adjacent_theta.at(index)) * (obs_p.right_most.first - x[index]) + CppAD::cos(adjacent_theta.at(index)) * (obs_p.right_most.second - x[index + N]);
      AD<double> fore_x = CppAD::cos(adjacent_theta.at(index)) * (obs_p.fore_most.first - x[index]) + CppAD::sin(adjacent_theta.at(index)) * (obs_p.fore_most.second - x[index + N]);
      AD<double> rear_x = CppAD::cos(adjacent_theta.at(index)) * (obs_p.rear_most.first - x[index]) + CppAD::sin(adjacent_theta.at(index)) * (obs_p.rear_most.second - x[index + N]);
      std::cout<<"  left_y "<<left_y<<", right_y "<<right_y<<", fore_x "<<fore_x<<", rear_x "<<rear_x<<std::endl;
      AD<double> ellipse_distance = (left_y + right_y) * (left_y + right_y) / ((left_y - right_y + 2 * lateral_threshold) * (left_y - right_y + 2 * lateral_threshold)) + (fore_x + rear_x) * (fore_x + rear_x) / ((fore_x - rear_x + 2 * longitude_threshold) * (fore_x - rear_x + 2 * longitude_threshold)) - 2.0;
      std::cout<<"  current ellipse_distance "<<ellipse_distance<<std::endl;
        // obs_cost += 1 / (1 + CppAD::exp(ellipse_distance * sigmoid_scale));
      obs_cost += CppAD::exp(-ellipse_distance * sigmoid_scale);
        // obs_cost = ellipse_distance;
    }
  }
  AD<double> current_smooth_cost = (adjacent_theta.front() - this->start_heading) * (adjacent_theta.front() - this->start_heading) / (distance.front() * distance.front());
  smooth_cost += current_smooth_cost * 1e+6;
  std::cout<<"current smooth cost "<<current_smooth_cost<<std::endl;
  current_smooth_cost = (adjacent_theta.back() - this->end_heading) * (adjacent_theta.back() - this->end_heading) / (distance.back() * distance.back());
  smooth_cost += current_smooth_cost * 1e+6;
  // std::cout<<"last theta "<<adjacent_theta.back()<<", last distance "<<distance.back()<<std::endl;
  std::cout<<"current smooth cost "<<current_smooth_cost<<std::endl;



    fg[0] = smooth_weight * smooth_cost + obs_weight * obs_cost;
    // fg[0] = smooth_cost;
    // fg[0] = obs_cost;
    // fg[1] = x[0];
    // fg[2] = x[N];
    // fg[3] = x[N-1];
    // fg[4] = x[2*N-1];
    return;
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

bool RefLineOPT::solve(CubicSpline& csp){
  typedef CPPAD_TESTVECTOR(double) Dvector;
  auto t1 = std::chrono::steady_clock::now();

  int max_iter = 5;
  int min_iter = 3;
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