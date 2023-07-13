#include "ref_line_opt.h"
#include <cassert>
#include <chrono>
#include <array>

const double lateral_threshold = 2;
const double longitude_threshold = 4;//TODO: extract params

const double smooth_weight = 100; //cost value about 10
const double obs_weight = 10;
const double deviation_weight = 10;
const double sigmoid_scale = 10;

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


FG_eval::FG_eval(RefLineOPT *opt){
    this->opt_ = opt;
    int n = opt_->points_x.size();
    this->obs_key_points.resize(n);
    double init_obs_cost = 0;
    std::vector<double> theta(n-1);
    std::vector<double> dist(n - 1);
    int obs_size = this->opt_->obs_.size();
    for(int i = 0; i < n; i++){//for each way points
        double relative_x;
        double relative_y;
        double MIN_X = 50;//TODO: extract params 
        double MAX_X = -50;
        double MIN_Y = 50;
        double MAX_Y = -50;
        this->obs_key_points.at(i).resize(obs_size);
        // double theta;
        // std::cout<<"For this point : "<<std::endl;
        for(int j = 0; j < obs_size; j++){//for each obstacle in this freesapce
            MIN_X = 50;//TODO: extract params 
            MAX_X = -50;
            MIN_Y = 50;
            MAX_Y = -50;
            // std::cout<<"  For this obstacle : "<<std::endl;
            KeyPoint& kp = this->obs_key_points.at(i).at(j);
            for(const std::pair<double, double>& p : this->opt_->obs_.at(j)){
                relative_x = cos(opt_->init_heading.at(i)) * (p.first - opt_->points_x.at(i)) + sin(opt_->init_heading.at(i)) * (p.second - opt_->points_y.at(i));
                relative_y = -sin(opt_->init_heading.at(i)) * (p.first - opt_->points_x.at(i)) + cos(opt_->init_heading.at(i)) * (p.second - opt_->points_y.at(i));
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
            // std::cout<<MAX_X<<" "<<MIN_X<<" "<<MAX_Y<<" "<<MIN_Y<<std::endl;
            // double current_obs_cost = std::pow((MIN_X + MAX_X) / (MAX_X - MIN_X + 15) , 2) + std::pow((MIN_Y + MAX_Y) / (MAX_Y - MIN_Y + 3), 2);
            // std::cout<<"current obs cost "<<current_obs_cost<<std::endl;
            // init_obs_cost += current_obs_cost;
            
        }
        if(i > 0){
            theta.at(i - 1) = atan2(opt_->points_y.at(i) - opt_->points_y.at(i - 1), opt_->points_x.at(i) - opt_->points_x.at(i - 1));
            dist.at(i - 1) = std::sqrt(std::pow(opt_->points_y.at(i) - opt_->points_y.at(i - 1), 2) + std::pow(opt_->points_x.at(i) - opt_->points_x.at(i - 1), 2));
        }
    }
    // std::cout<<"&&&&&&&&&&&&&initial obs cost "<<init_obs_cost<<std::endl;
    // init_obs_cost = 0;

    double init_smooth_cost = 0;
    for(int i = 1; i < n-1; i++){
        double current_smooth_cost = std::pow(theta.at(i) - theta.at(i - 1), 2) / std::pow(dist.at(i) + dist.at(i - 1), 2);
        std::cout<<"current smooth cost "<<current_smooth_cost<<std::endl;
        init_smooth_cost += current_smooth_cost;
        // for(const KeyPoint& obs_p : this->obs_key_points.at(i)){
        //     double left_y = -sin(opt_->init_heading.at(i)) * (obs_p.left_most.first - this->opt_->points_x.at(i)) + cos(opt_->init_heading.at(i)) * (obs_p.left_most.second - this->opt_->points_y.at(i));
        //     double right_y = -sin(opt_->init_heading.at(i)) * (obs_p.right_most.first - this->opt_->points_x.at(i)) + cos(opt_->init_heading.at(i)) * (obs_p.right_most.second - this->opt_->points_y.at(i));
        //     double fore_x = cos(opt_->init_heading.at(i)) * (obs_p.fore_most.first - this->opt_->points_x.at(i)) + sin(opt_->init_heading.at(i)) * (obs_p.fore_most.second - this->opt_->points_y.at(i));
        //     double rear_x = cos(opt_->init_heading.at(i)) * (obs_p.rear_most.first - this->opt_->points_x.at(i)) + sin(opt_->init_heading.at(i)) * (obs_p.rear_most.second - this->opt_->points_y.at(i));
        //     std::cout<<fore_x<<" "<<rear_x<<" "<<left_y<<" "<<right_y<<std::endl;
        //     double current_obs_cost = std::pow((fore_x + rear_x) / (fore_x - rear_x + 15) , 2) + std::pow((left_y + right_y) / (left_y - right_y + 3), 2);
        //     std::cout<<"current obs cost "<<current_obs_cost<<std::endl;
        //     init_obs_cost += current_obs_cost;     
        // }
    }
    double current_smooth_cost = std::pow(opt_->init_heading.front() - theta.front(), 2) / std::pow(dist.front(), 2);
    init_smooth_cost += current_smooth_cost;
    std::cout<<"current smooth cost "<<current_smooth_cost<<std::endl;
    current_smooth_cost = std::pow(opt_->init_heading.back() - theta.back(), 2) / std::pow(dist.back(), 2);
    std::cout<<"last theta "<<theta.back()<<", last distance "<<dist.back();
    std::cout<<"current smooth cost "<<current_smooth_cost<<std::endl;
    init_smooth_cost += current_smooth_cost;
    std::cout<<"&&&&&&&&&&&&&&&initial smooth cost "<<init_smooth_cost<<std::endl;
    // std::cout<<"&&&&&&&&&&&&&&&initial obs cost "<<init_obs_cost<<std::endl;
}


void FG_eval::operator()(ADvector& fg, const ADvector& x)
{
    // assert(fg.size() == 5);
    assert(x.size() == this->opt_->points_x.size() * 2);
    // variables
    // AD<double> x1 = x[0];
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
        deviation_cost += CppAD::log10( 200 * ( (x[index] - this->opt_->points_x.at(index)) * (x[index] - this->opt_->points_x.at(index)) + (x[index + N] - this->opt_->points_y.at(index)) * (x[index + N] - this->opt_->points_y.at(index)) + 1 ));
        //AD<double> current_smooth_cost = (adjacent_theta.at(index) - adjacent_theta.at(index - 1)) * (adjacent_theta.at(index) - adjacent_theta.at(index - 1)) / ( (distance.at(index - 1) + distance.at(index)) * (distance.at(index - 1) + distance.at(index)));
        
        AD<double> current_smooth_cost = (adjacent_vec.at(index).at(0) - adjacent_vec.at(index-1).at(0)) * (adjacent_vec.at(index).at(0) - adjacent_vec.at(index-1).at(0)) + (adjacent_vec.at(index).at(1) - adjacent_vec.at(index-1).at(1)) * (adjacent_vec.at(index).at(1) - adjacent_vec.at(index-1).at(1));
        std::cout<<"current smooth cost "<<current_smooth_cost<<std::endl;
        smooth_cost += current_smooth_cost;
        for(const KeyPoint& obs_p : this->obs_key_points.at(index)){
            AD<double> left_y = -CppAD::sin(adjacent_theta.at(index)) * (obs_p.left_most.first - x[index]) + CppAD::cos(adjacent_theta.at(index)) * (obs_p.left_most.second - x[index + N]);
            AD<double> right_y = -CppAD::sin(adjacent_theta.at(index)) * (obs_p.right_most.first - x[index]) + CppAD::cos(adjacent_theta.at(index)) * (obs_p.right_most.second - x[index + N]);
            AD<double> fore_x = CppAD::cos(adjacent_theta.at(index)) * (obs_p.fore_most.first - x[index]) + CppAD::sin(adjacent_theta.at(index)) * (obs_p.fore_most.second - x[index + N]);
            AD<double> rear_x = CppAD::cos(adjacent_theta.at(index)) * (obs_p.rear_most.first - x[index]) + CppAD::sin(adjacent_theta.at(index)) * (obs_p.rear_most.second - x[index + N]);
            AD<double> ellipse_distance = (left_y + right_y) * (left_y + right_y) / ((left_y - right_y + 2 * lateral_threshold) * (left_y - right_y + 2 * lateral_threshold)) + (fore_x + rear_x) * (fore_x + rear_x) / ((fore_x - rear_x + 2 * longitude_threshold) * (fore_x - rear_x + 2 * longitude_threshold)) - 1;
            std::cout<<"current ellipse_distance "<<ellipse_distance<<std::endl;
            // obs_cost += 1 / (1 + CppAD::exp(ellipse_distance * sigmoid_scale));
            obs_cost += CppAD::exp(-ellipse_distance * sigmoid_scale);
            // obs_cost = ellipse_distance;
        }
    }
    AD<double> current_smooth_cost = (adjacent_theta.front() - this->opt_->init_heading.front()) * (adjacent_theta.front() - this->opt_->init_heading.front()) / (distance.front() * distance.front());
    smooth_cost += current_smooth_cost * 1000;
    std::cout<<"current smooth cost "<<current_smooth_cost<<std::endl;
    current_smooth_cost = (adjacent_theta.back() - this->opt_->init_heading.back()) * (adjacent_theta.back() - this->opt_->init_heading.back()) / (distance.back() * distance.back());
    smooth_cost += current_smooth_cost * 1000;
    // std::cout<<"last theta "<<adjacent_theta.back()<<", last distance "<<distance.back()<<std::endl;
    std::cout<<"current smooth cost "<<current_smooth_cost<<std::endl;



    fg[0] = smooth_weight * smooth_cost + obs_weight * obs_cost + deviation_weight * deviation_cost;
    // fg[0] = smooth_cost;
    // fg[0] = obs_cost;
    // fg[1] = x[0];
    // fg[2] = x[N];
    // fg[3] = x[N-1];
    // fg[4] = x[2*N-1];
    return;
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

RefLineOPT::RefLineOPT(std::vector<Point> points, const std::vector<Obstacle>& obs){
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
}

RefLineOPT::RefLineOPT(std::vector<TrackPoint> track_points, const std::vector<Obstacle>& obs){
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
}

bool RefLineOPT::solve(){
    auto t1=std::chrono::steady_clock::now();
    bool ok = false;
    int N = points_x.size();
    if(N != points_y.size() || N != init_heading.size())
        throw("Unequal number of coordinates or heading before optimization");

    size_t i;
    typedef CPPAD_TESTVECTOR(double) Dvector;

    size_t nx = N * 2; // number of varibles
    size_t ng = 0; // number of constraints
    Dvector x0(nx), xl(nx), xu(nx); // initial condition of varibles
    Dvector gl(ng), gu(ng);// lower and upper bounds for constraints
    
    for(i = 0; i < N; i++){
        x0[i] = this -> points_x.at(i);
        x0[i + N] = this -> points_y.at(i);
        xl[i] = this -> points_x.at(i) - 5.0;
        xl[i + N] = this -> points_y.at(i) - 5.0;
        xu[i] = this -> points_x.at(i) + 5.0;
        xu[i + N] = this -> points_y.at(i) + 5.0;
    }
    xl[0] = this->points_x.front();
    xu[0] = this->points_x.front();
    xl[N-1] = this->points_x.back();
    xu[N-1] = this->points_x.back();
    xl[N] = this->points_y.front();
    xu[N] = this->points_y.front();
    xl[2 * N - 1] = this->points_y.back();
    xu[2 * N - 1] = this->points_y.back();

    // gl[0] = points_x.front();
    // gu[0] = points_x.front();
    // gl[1] = points_y.front();
    // gu[1] = points_y.front();
    // gl[2] = points_x.back();
    // gu[2] = points_x.back();
    // gl[3] = points_y.back();
    // gu[3] = points_y.back();

    // object that computes objective and constraints
    FG_eval fg_eval(this);
    // options
    std::string options;
    // turn off any printing
    //options += "Retape  true\n";
    options += "Integer print_level  5\n";
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


    CppAD::ipopt::solve_result<Dvector> solution; // solution
    CppAD::ipopt::solve<Dvector, FG_eval>(options, x0, xl, xu, gl, gu, fg_eval, solution); // solve the problem

    if(solution.status == CppAD::ipopt::solve_result<Dvector>::status_type::success)
      ok = true;
    else{
      std::cout<<"\033[31;1mThe Solution was NOT Successfull\033[0m"<<std::endl;
      std::cout<<"Status "<<EnumStrings[solution.status]<<std::endl;
    }

    auto t2=std::chrono::steady_clock::now();
    double dr_s=std::chrono::duration<double>(t2-t1).count();
    //毫秒级
    double dr_ms=std::chrono::duration<double,std::milli>(t2-t1).count();
    //微妙级
    double dr_us=std::chrono::duration<double,std::micro>(t2-t1).count();
    //纳秒级
    double dr_ns=std::chrono::duration<double,std::nano>(t2-t1).count(); 

    std::cout<<"ns "<<dr_ns<<", us "<<dr_us<<", ms "<<dr_ms<<", s"<<dr_s<<std::endl;

    
    for(i = 0; i < N; i++){
        this->points_x.at(i) = solution.x[i];
        this->points_y.at(i) = solution.x[i + N];
    }

    correctHeading();
      

    std::cout<<"Final result F "<<solution.obj_value<<std::endl;
    ok = true;
    return ok;
}

void RefLineOPT::correctHeading(){
    this->init_heading.resize(points_x.size());
    for(int i = 1; i < points_x.size() - 1; i++)
        init_heading.at(i) = atan2(points_y.at(i + 1) - points_y.at(i - 1), points_x.at(i + 1) - points_x.at(i - 1));
}

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

// std::vector<Point> RefLineOPT::getPoints(){
//     std::vector<int> black_list(this->points_x.size(), false);
//     double min_thre = 0.15;
//     extarctKeyPoints(this->points_x, this->points_y, min_thre, 0, this->points_x.size() - 1, black_list);
//     std::vector<Point> points = {};
//     for(int i = 0; i < points_x.size(); i++){
//         if(!black_list.at(i))
//             points.emplace_back(points_x.at(i), points_y.at(i));      
//     }
//     return points;

// }