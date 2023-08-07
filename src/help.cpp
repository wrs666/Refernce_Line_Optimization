#include "map_base.h"
#include "ref_line_opt.h"
#include <cmath>


void Junction::ref_lines_generation(const std::vector<Obstacle>& obs){
  for(const Interface& in : this->entries){
    TrackPoint start(in.second.p, in.first.geometry.heading, 0);
    for(const Interface& out : this->exits){

      auto t1 = std::chrono::steady_clock::now();
      RefLineOPT ref_line_opt(Pose(in.second.p, in.first.geometry.heading), Pose(out.second.p, out.first.geometry.heading), obs);
      auto t2 = std::chrono::steady_clock::now();
      std::shared_ptr<RefLine> optimized_ref_line = ref_line_opt.getPoints();
 
      
      double dr_ms=std::chrono::duration<double,std::milli>(t2-t1).count();
      std::cout<<"Time consumption of Whole process: "<<dr_ms<<"ms."<<std::endl;


      std::shared_ptr<RefLine> init_ref_line = quinticInterploration(Pose(in.second.p,in.first.geometry.heading), Pose(out.second.p, out.first.geometry.heading), 0.1);
      
      this->ref_lines.push_back((*init_ref_line));
      this->optimized_ref_lines.push_back((*optimized_ref_line));

      std::vector<Point> optimized_key_points = ref_line_opt.getKeyPoints();
      
      std::vector<Point> init_key_points {};
      for(int i = 0; i < ref_line_opt.init_key_x.size(); i++){
        init_key_points.emplace_back(ref_line_opt.init_key_x.at(i), ref_line_opt.init_key_y.at(i));
      }
      this->initial_key_points.push_back(init_key_points);
      this->optimized_ref_keys.push_back(optimized_key_points);

    }
  }
}

