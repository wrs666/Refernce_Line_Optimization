#include "help.h"
#include "ref_line_opt.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <cmath>

using namespace Eigen;

double sampling_interval = 2.0;

std::vector<TrackPoint> spline_interploration(TrackPoint p1, TrackPoint p2, double point_margin){
  std::vector<TrackPoint> curve_points;
  curve_points.resize(0);

  double x1, y1, theta1, kappa1, x2, y2, theta2, kappa2;
  x1 = p1.x_;
  y1 = p1.y_;
  theta1 = p1.theta_;
  kappa1 = p1.curvature_;
  x2 = p2.x_;
  y2 = p2.y_;
  theta2 = p2.theta_;
  kappa2 = p2.curvature_;

  std::cout<<"x1, y1, theta1 : "<< x1 << " " << y1 << " " << theta1 <<std::endl;
  std::cout<<"x2, y2, theta2 : "<< x2 << " " << y2 << " " << theta2 <<std::endl;
  
  double s_distance = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));

  curve_points.reserve((int)ceil(s_distance / point_margin) + 1);

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

  double kappa_start = (vx1 * ay1 - ax1 * vy1) / pow( (pow(vx1, 2) + pow(vy1, 2) ), 1.5);
  double kappa_end = (vx2 * ay2 - ax2 * vy2) / pow( (pow(vx2, 2) + pow(vy2, 2) ), 1.5);

  std::cout<<"kappa start : "<< kappa_start << std::endl;
  std::cout<<"kappa end : "<< kappa_end <<std::endl;

  curve_points.push_back(p1);

  double s = point_margin;
  TrackPoint p;
  double vx, vy, ax, ay, theta;

  while(s < s_distance)
  {
    VectorXd s_power(6);
    VectorXd coef_x(6);
    VectorXd coef_y(6);
    VectorXd coef_vx(6);
    VectorXd coef_vy(6);
    VectorXd coef_ax(6);
    VectorXd coef_ay(6);
      
    coef_x << a[0], a[1], a[2], a[3], a[4], a[5];
    coef_y << b[0], b[1], b[2], b[3], b[4], b[5];
    coef_vx << a[1], 2 * a[2], 3 * a[3], 4 * a[4], 5 * a[5], 0;
    coef_vy << b[1], 2 * b[2], 3 * b[3], 4 * b[4], 5 * b[5], 0;
    coef_ax << 2 * a[2], 6 * a[3], 12 * a[4], 20 * a[5], 0, 0;
    coef_ay << 2 * b[2], 6 * b[3], 12 * b[4], 20 * b[5], 0, 0;


    s_power << 1, s, pow(s, 2), pow(s, 3), pow(s, 4), pow(s, 5);
    
    p.x_ = coef_x.dot(s_power);
    p.y_ = coef_y.dot(s_power);
    vx = coef_vx.dot(s_power);
    vy = coef_vy.dot(s_power);
    ax = coef_ax.dot(s_power);
    ay = coef_ay.dot(s_power);

    p.theta_ = std::atan2(vy, vx);

    // if(vx < 0)
    //   if(vy < 0)
    //     p.theta_ = theta - PI;
    //   else
    //     p.theta_ = theta + PI;
    // else
    //   p.theta_ = theta; 

    p.curvature_ = (vx * ay - ax * vy) / pow( (pow(vx, 2) + pow(vy, 2) ), 1.5);

    curve_points.push_back(p);
    s = s + point_margin;//2*2^1/2 / PI = 0.9003.....
  }

  if(s_distance - s > point_margin / 2)
      curve_points.push_back(p2);
  std::cout<<"Finish one spline interplorartion, points size is "<<(int)curve_points.size()<<std::endl;
  return curve_points;
}

void Junction::ref_lines_generation(const std::vector<Obstacle>& obs){
  for(const Interface& in : this->entries){
    TrackPoint start(in.second.p, in.first.geometry.heading, 0);
    for(const Interface& out : this->exits){
      TrackPoint end(out.second.p, out.first.geometry.heading, 0);
      std::vector<TrackPoint> init_ref_line = spline_interploration(start, end, sampling_interval);
      // ref_line simple_ref_line;
      // simple_ref_line.resize(init_ref_line.size());
      // for(int i = 0; i < init_ref_line.size(); i++){
      //   simple_ref_line.at(i) = init_ref_line.at(i);
      // }
      RefLineOPT ref_line_opt(init_ref_line, obs);
      bool solved = ref_line_opt.solve();
      solved = ref_line_opt.solve();
      solved = ref_line_opt.solve();
      solved = ref_line_opt.solve();
      solved = ref_line_opt.solve();
      solved = ref_line_opt.solve();
      solved = ref_line_opt.solve();
      if(solved){
        this->optimized_ref_lines.push_back(ref_line_opt.getPoints());
      }
      this->ref_lines.push_back(init_ref_line);
    }
  }
}

