#include "kalman_filter.h"
#include "tools.h"
#include <math.h>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;


/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}
KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
  
}

void KalmanFilter::Predict() {
  x_ = F_ * x_; // + u_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}


/**
 * Updates the state by using Kalman Filter equations
 * 
 * @param z cortesian position from the latests measurment
 */
void KalmanFilter::Update(const VectorXd &z) {

  VectorXd y = z - H_ * x_;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;
  
  // new state
  MatrixXd I = MatrixXd::Identity(4, 4);
  x_ = x_ + (K * y);
  P_ = (I - K * H_) * P_;
} 

/**
 * Updates the state by using Extended Kalman Filter equations
 * 
 * @param z polar position from the latest measrument
 */
void KalmanFilter::UpdateEKF(const VectorXd &z) {

  Tools tools;
  
  MatrixXd Hj = MatrixXd(3, 4);
  Hj = tools.CalculateJacobian(x_);
  
  float l = sqrt(pow(x_[0], 2) + pow(x_[1], 2));
 
  if (l < 0.001) {
    cout << "Devision by zero in the Radar update" << endl;
    return;
  }

  VectorXd y = z - tools.PolarCoordinates(x_);
  y = tools.ConformPolarPosition(y);

  MatrixXd Hjt = Hj.transpose();
  MatrixXd S = Hj * P_ * Hjt + R_;
  MatrixXd K =  P_ * Hjt * S.inverse();
  
  // new state
  MatrixXd I = MatrixXd::Identity(4, 4);
  x_ = x_ + (K * y);
  P_ = (I - K * Hj) * P_; 
}
