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

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * Updates the state by using Kalman Filter equations
   */
  
  VectorXd y = z - H_ * x_;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;
  
  // new state
  MatrixXd I = MatrixXd::Identity(4, 4);
  // x_ = x_ + (K * y);
  // P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * Updates the state by using Extended Kalman Filter equations
   */
  Tools tools;

  VectorXd z_ = VectorXd(3);
  z_ << z[0], z[1], z[2];

  MatrixXd Hj = MatrixXd(3, 4);
  Hj = tools.CalculateJacobian(x_);
  
  /** 
  cout << "z_ : " << z_ << endl;
  cout << "Hj is calculated OK" << endl;
  cout << "H Jacobian: " << endl << Hj << endl;
  cout << "x_ " << endl << x_ << endl;
  */
  
  VectorXd z_dash = VectorXd(3);
  
  float l = sqrt(pow(x_[0], 2) + pow(x_[1], 2));
 
  if (l < 0.001) {
    cout << "Devision by zero in the Radar update" << endl;
    return;
  }
  
  /**
  l = z[0];
  float angle = z[1];
  float l_velocity = z[2];
  x_ << cos(angle) * l,
             sin(angle) * l,
             cos(angle) * l_velocity,
             sin(angle) * l_velocity; 
  return; 
  */
  
  z_dash << l,
            atan2(x_[1], x_[0]), 
            (x_[0] * x_[2] + x_[1] * x_[3]) / l;
  
  cout << "z_ : " << z_ << endl;
  cout << "z_dash : " << z_dash << endl;
  
 
  VectorXd y = z_ - z_dash;

  // adjust the angle to the (-pi, pi) range expected 
  // by the calculation in the Kalman Filter
  if (y[1] > M_PI || y[1] < -M_PI) {
    y[1] = atan2(sin(y[1]), cos(y[1]));
  }

  MatrixXd Hjt = Hj.transpose();
  MatrixXd S = Hj * P_ * Hjt + R_radar_;
  MatrixXd K =  P_ * Hjt * S.inverse();
  
  // new state
  MatrixXd I = MatrixXd::Identity(4, 4);
  x_ = x_ + (K * y);
  P_ = (I - K * Hj) * P_; 
  
}
