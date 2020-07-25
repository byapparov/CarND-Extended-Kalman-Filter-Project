#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * Calculates the RMSE.
   */
  
   VectorXd rmse(4);
   rmse << 0,0,0,0;
   
   // check the validity of the following inputs:
   //  * the estimation vector size should not be zero
   //  * the estimation vector size should equal ground truth vector size
   if (estimations.size() == 0 || 
         ground_truth.size() == 0 || 
         ground_truth.size() != estimations.size()) {
      return rmse;
   }
   
   // accumulate squared residuals
   
   // int observations = estimations.size();
   VectorXd err(4);
   VectorXd err_square(4);
   
   for (int i=0; i < estimations.size(); ++i) {
      
      err = estimations[i] - ground_truth[i];
      err_square = err.array() * err.array();
      rmse += err_square;
      
   }
   // calculating the mean
   rmse = rmse / estimations.size();
   
   // calculating the squared root
   rmse = rmse.array().sqrt();

   return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * Calculates a Jacobian Matrix.
   */
  
   MatrixXd Hj(3, 4);
   // recover state parameters
   float px = x_state(0);
   float py = x_state(1);
   float vx = x_state(2);
   float vy = x_state(3);
   
    // check division by zero
   
   if(px == 0 && py == 0) {
      cout << "CalculateJacobian() - Error - Division by Zero";
      return Hj;
   }
   
   // compute the Jacobian matrix
   float l_2 = px * px + py * py;
   float l = sqrt(l_2);
   
   // p
   Hj(0,0) = px / l;
   Hj(0,1) = py / l;
   Hj(0,2) = 0;
   Hj(0,3) = 0;
   
   // phi 
   Hj(1,0) = -py / l_2;
   Hj(1,1) = px / l_2;
   Hj(1,2) = 0;
   Hj(1,3) = 0;
   
   // p_speed
   Hj(2,0) = py * (vx * py - vy * px) / pow(l_2, 3/2);
   Hj(2,1) = px * (vy * px - vx * py) / pow(l_2, 3/2);
   Hj(2,2) = px / l;
   Hj(2,3) = py / l;
   
   return Hj;
}