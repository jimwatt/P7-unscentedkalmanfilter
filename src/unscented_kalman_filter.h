#pragma once

#include "Eigen/Dense"

namespace unscented_kalman_filter {

// Straight forwad Extended Kalman filter predict step
void Predict(Eigen::VectorXd &x, Eigen::MatrixXd &P, const Eigen::MatrixXd &F, const Eigen::MatrixXd &Q);

// Straight forward Extended Kalman filter update step
void Update(Eigen::VectorXd &x, Eigen::MatrixXd &P, 
  const Eigen::MatrixXd &H, const Eigen::MatrixXd &R, 
  const Eigen::VectorXd& yk);


}   // namespace unscented_kalman_filter


