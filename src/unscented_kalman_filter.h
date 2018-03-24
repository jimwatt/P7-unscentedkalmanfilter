#pragma once

#include "Eigen/Dense"

// Uses the UKF to update the given state.  These functions use no global variables, and do not maintain state.

namespace unscented_kalman_filter {

void UpdateState(const double delta_t, const Eigen::VectorXd& process_noise_std, 
	const Eigen::VectorXd& meas_noise_std, const Eigen::VectorXd& z, 
	Eigen::VectorXd measurement_function(const Eigen::VectorXd&),
	Eigen::VectorXd& x, Eigen::MatrixXd& P);

}   // namespace unscented_kalman_filter


