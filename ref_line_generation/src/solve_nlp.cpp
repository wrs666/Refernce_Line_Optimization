#include "solve_nlp.h"

using namespace std;
using CppAD::AD;
void FG_eval::operator()(ADvector& fg, const ADvector& x)
{
    assert(x.size() == this->origin_radius.size() - 1);
    assert(fg.size() == x.size() + 1);
    // variables
    // AD<double> x1 = x[0];

    double car_fake_diagonal = sqrt(pow(car_width, 2) + pow(2 * car_length_front, 2));
    double radius, distance, curvature;

    vector<AD<double>> kappa_segment(x.size());
    vector<AD<double>> kappa_slope_segment(x.size() - 1);

    //fill kappa segment
    kappa_segment.front() = heading_diff.front() / this->points_distance + (sign(road_curvature.at(1)) / pow(this->points_distance, 2)) * x[0];
    kappa_segment.at(1) = heading_diff.at(1) / this->points_distance + (sign(road_curvature.at(2)) / pow(this->points_distance, 2)) * x[1] - (2 * sign(road_curvature.at(1)) / pow(this->points_distance, 2)) * x[0];
    kappa_slope_segment.front() = (kappa_segment.at(1) - kappa_segment.front()) / this->points_distance;
    
    for(size_t k = 2; k < x.size(); k++){
        kappa_segment.at(k) += heading_diff.at(k) / this->points_distance + (sign(road_curvature.at(k + 1)) / pow(this->points_distance, 2)) * x[k] + (sign(road_curvature.at(k - 1)) / pow(this->points_distance, 2)) * x[k - 2] -  (2 * sign(road_curvature.at(k)) / pow(this->points_distance, 2)) * x[k - 1];
        kappa_slope_segment.at(k - 1) = (kappa_segment.at(k) - kappa_segment.at(k - 1)) / this->points_distance;
    }

    double K = this->pcl_.getSlopeWeight();

    fg[0] = 0.0;
    // fg[0] += heading_diff.front() / this->points_distance + (sign(road_curvature.at(1)) / pow(this->points_distance, 2)) * x[0];
    // fg[0] += heading_diff.at(1) / this->points_distance + (sign(road_curvature.at(2)) / pow(this->points_distance, 2)) * x[1] - (2 * sign(road_curvature.at(1)) / pow(this->points_distance, 2)) * x[0];
    for(size_t i = 0; i < x.size() - 1; i++){
        //objective function
        //fg[0] += K * kappa_slope_segment.at(i) * kappa_slope_segment.at(i);
        //fg[0] += kappa_segment.at(i) * kappa_segment.at(i);

        // constraints
        fg[i + 1] = kappa_segment.at(i);
    }

    for(size_t k = 1; k < x.size(); k++){
        //x[k]作用于第k+1个点
        radius = this->origin_radius.at(k + 1);
        curvature = road_curvature.at(k + 1);

        if(curvature < 0.01)//toFix : extract global prameters
            fg[0] += 4*x[k]*x[k];
        else{
            fg[0] += (radius*2 +(x[k] - radius)*cos((sign(road_curvature.at(k + 1)) * x[k] - sign(road_curvature.at(k)) * x[k-1])/this->points_distance) + car_width/2 - sqrt(pow(car_fake_diagonal, 2)/4 + (radius - x[k])*(radius - x[k]) + car_fake_diagonal*(radius-x[k])*sin(alpha + (sign(road_curvature.at(k + 1)) * x[k] - sign(road_curvature.at(k)) * x[k-1])/this->points_distance)))
                    * (radius*2 +(x[k] - radius)*cos((sign(road_curvature.at(k + 1)) * x[k] - sign(road_curvature.at(k)) * x[k-1])/this->points_distance) + car_width/2 - sqrt(pow(car_fake_diagonal, 2)/4 + (radius - x[k])*(radius - x[k]) + car_fake_diagonal*(radius-x[k])*sin(alpha + (sign(road_curvature.at(k + 1)) * x[k] - sign(road_curvature.at(k)) * x[k-1])/this->points_distance)));
        }
    }
    //fg[0] += kappa_segment.back() * kappa_segment.back();
    fg[x.size()] = kappa_segment.back();
    return;
}